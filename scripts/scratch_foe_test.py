#!/usr/bin/env python3
"""Inspect clusters in scene13_01 and check epipolar consistency with camera velocity."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path("python/scripts").resolve()))
sys.path.insert(0, str(Path("python").resolve()))

import numpy as np
import json
import argparse

import evimo2_common as ec
import rotcomp as rc
import object_motion as om
import imo_labels as il
from survey_sequences import pose_series
from evimo2.detector import detect_imo
from evimo2.geometry import qrot, qconj

SEQ_PATH = Path("data/npz_left_camera_imo/left_camera/imo/train/scene13_dyn_test_01_000000")

def run():
    print(f"Loading data from {SEQ_PATH}...")
    with open(SEQ_PATH / "dataset_intrinsics.json") as f:
        intr = json.load(f)
    info = ec.load_info(SEQ_PATH)
    t_ev, xy_ev, p_ev = ec.load_events(SEQ_PATH)
    
    # IMU calibration
    raw_imu = np.load(SEQ_PATH / "dataset_imu.npy")
    cal_ts = raw_imu[:, 0]
    cal_w = rc.calibrated_gyro(raw_imu[:, 1:4])
    
    # Camera trajectory from Vicon (GT)
    traj = info["meta"].get("full_trajectory") or []
    cam_ts, cam_q, cam_p = pose_series(traj, "cam")
    
    # Kinematics
    _, tracks = om.load_tracks(info)
    kin = il.kinematics(tracks, "R*S", 0.025)
    
    # GT masks
    z_mask = np.load(SEQ_PATH / "dataset_mask.npz")
    
    # Frames
    series_dict = ec.find_series(info)
    frames_list = series_dict.get("meta/frames") or series_dict.get("frames") or []
    fts = {int(d.get("fid", d.get("frame_id", i))): float(d["ts"]) for i, d in enumerate(frames_list) if "ts" in d}
    
    W = 0.040 # 40 ms
    d_lag = -0.0309 # calibrated lag for scene13
    
    # Pick frame 60 where camera is translating at 0.165 m/s
    fid = 60
    tsf = fts[fid]
    win = rc.prep_window(tsf, W, t_ev, xy_ev, intr)
    xs, ys = rc.warp(win, cal_ts, cal_w, d_lag)
    te = win["te"]
    
    # Camera velocity at tsf
    idx = np.searchsorted(cam_ts, tsf)
    dt = cam_ts[idx+2] - cam_ts[idx-2]
    v_world = (cam_p[idx+2] - cam_p[idx-2]) / dt
    q_wc = cam_q[idx:idx+1]
    q_cw = qconj(q_wc)
    v_cam = qrot(q_cw, v_world[None, :])[0]
    speed = np.linalg.norm(v_cam)
    
    print(f"\n================ Frame {fid} at t={tsf:.3f}s ================")
    print(f"Camera linear velocity v_cam: [{v_cam[0]:.3f}, {v_cam[1]:.3f}, {v_cam[2]:.3f}] m/s (|v|={speed:.3f} m/s)")
    
    # Run v1 detector (r=3, tau=0.3)
    pred_imo, n_clusters = detect_imo(xs, ys, te, tsf, blur=3, tau=0.3, support=1.0, min_cells=8)
    print(f"Detector v1: found {n_clusters} clusters, {pred_imo.sum()} IMO events.")
    
    # GT regions
    from imo_detect import regions_for_frame
    a_opts = argparse.Namespace(slow_v=0.05, fast_v=0.3, w_thr=1.0, dilate=6)
    excl, fmask = regions_for_frame(z_mask, fid, tsf, kin, a_opts)
    
    pix_x = np.clip(win["pix"][:, 0].astype(int), 0, 639)
    pix_y = np.clip(win["pix"][:, 1].astype(int), 0, 479)
    posm = fmask[pix_y, pix_x]
    negm = (~excl[pix_y, pix_x])
    
    # Early vs Late masks
    m_early = te < tsf
    m_late = te >= tsf
    
    # Analyze connected components / clusters
    # Let's map events to clusters
    # In python detector: CCL is done on 4x4 cells
    H, W_px = 480, 640
    # Let's inspect each connected component
    from evimo2.detector import label_cc
    # Recalculate cell grid
    cell_grid = np.zeros((H // 4, W_px // 4), bool)
    for i in range(len(xs)):
        if pred_imo[i]:
            cx = int(xs[i]) // 4
            cy = int(ys[i]) // 4
            if 0 <= cx < (W_px // 4) and 0 <= cy < (H // 4):
                cell_grid[cy, cx] = True
                
    labeled, areas = label_cc(cell_grid)
    num_comp = len(areas) - 1
    print(f"Connected components on cell grid: {num_comp}")
    
    print("\n" + "="*85)
    print(f"{'Klaster':<8s} {'Rozmiar':<9s} {'Środek (x, y)':<16s} {'Typ GT':<14s} {'Wektor przemieszczenia':<22s} {'Kąt z linią epipolarną':<20s}")
    print("="*85)
    
    for c_id in range(1, num_comp + 1):
        # Find events belonging to this component
        mask_cells = (labeled == c_id)
        # Bounding box of cells
        cy_idx, cx_idx = np.where(mask_cells)
        x_min, x_max = cx_idx.min() * 4, (cx_idx.max() + 1) * 4
        y_min, y_max = cy_idx.min() * 4, (cy_idx.max() + 1) * 4
        
        in_bbox = (xs >= x_min) & (xs <= x_max) & (ys >= y_min) & (ys <= y_max) & pred_imo
        n_ev = int(in_bbox.sum())
        if n_ev < 100:
            continue
            
        # Ground truth status
        n_pos = int((in_bbox & posm).sum())
        n_neg = int((in_bbox & negm).sum())
        if n_pos > n_neg:
            gt_type = "PRAWDZIWE IMO"
        else:
            gt_type = "PARALAKSA (TŁO)"
            
        # Calculate cluster centroid
        cx_mean = float(xs[in_bbox].mean())
        cy_mean = float(ys[in_bbox].mean())
        
        # Calculate Early centroid vs Late centroid
        ev_e = in_bbox & m_early
        ev_l = in_bbox & m_late
        if ev_e.sum() < 20 or ev_l.sum() < 20:
            continue
            
        xe_mean, ye_mean = float(xs[ev_e].mean()), float(ys[ev_e].mean())
        xl_mean, yl_mean = float(xs[ev_l].mean()), float(ys[ev_l].mean())
        
        dx = xl_mean - xe_mean
        dy = yl_mean - ye_mean
        shift_len = np.hypot(dx, dy)
        
        # Theoretical epipolar direction at (cx_mean, cy_mean):
        # In normalized pinhole coordinates:
        x_norm = (cx_mean - intr["cx"]) / intr["fx"]
        y_norm = (cy_mean - intr["cy"]) / intr["fy"]
        # Background apparent velocity:
        # dx_norm/dt = (-v_cam[0] + x_norm * v_cam[2])
        # dy_norm/dt = (-v_cam[1] + y_norm * v_cam[2])
        # In pixels:
        u_dot = intr["fx"] * (-v_cam[0] + x_norm * v_cam[2])
        v_dot = intr["fy"] * (-v_cam[1] + y_norm * v_cam[2])
        epi_len = np.hypot(u_dot, v_dot)
        
        if epi_len > 1e-6 and shift_len > 0.2:
            u_epi_norm = np.array([u_dot, v_dot]) / epi_len
            u_shift_norm = np.array([dx, dy]) / shift_len
            # Dot product (cosine of angle)
            cos_theta = np.clip(np.dot(u_shift_norm, u_epi_norm), -1.0, 1.0)
            angle_deg = np.degrees(np.arccos(cos_theta))
            # Or perpendicular deviation (sine of angle)
            sin_theta = abs(np.cross(u_shift_norm, u_epi_norm))
            
            vec_str = f"({dx:+.1f}, {dy:+.1f}) [{shift_len:.1f}px]"
            deg_str = f"kat: {angle_deg:5.1f} deg (sin={sin_theta:.2f})"
            print(f"#{c_id:<7d} {n_ev:<9d} ({cx_mean:5.1f}, {cy_mean:5.1f})   {gt_type:<14s} {vec_str:<22s} {deg_str:<20s}")

if __name__ == "__main__":
    run()
