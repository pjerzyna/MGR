#!/usr/bin/env python3
"""Ewaluacja detektora v2 z filtrem epipolarnym (FOE) na sekwencji translacyjnej scene13_01."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path("python").resolve()))

import json
import csv
import argparse
import numpy as np

import evimo2 as ev
from evimo2.geometry import qrot, qconj

SEQ_PATH_DEFAULT = Path("data/npz_left_camera_imo/left_camera/imo/train/scene13_dyn_test_01_000000")
OUT_DIR_DEFAULT = Path("results/v2_epipolar_foe")


def main():
    parser = argparse.ArgumentParser(description="Ewaluacja detektora v2 FOE na sekwencji translacyjnej.")
    parser.add_argument("--seq", type=Path, default=SEQ_PATH_DEFAULT, help="Ścieżka do sekwencji EVIMO2")
    parser.add_argument("--out", type=Path, default=OUT_DIR_DEFAULT, help="Katalog wyjściowy na wyniki CSV/raporty")
    parser.add_argument("--angle-thresh", type=float, default=35.0, help="Próg kątowy w stopniach dla linii epipolarnych")
    parser.add_argument("--min-speed", type=float, default=0.05, help="Minimalna prędkość translacji kamery [m/s]")
    args = parser.parse_args()

    seq_path = args.seq
    out_dir = args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Ładowanie danych z {seq_path}...")
    with open(seq_path / "dataset_intrinsics.json") as f:
        intr = json.load(f)
    info = ev.load_info(seq_path)
    t_ev, xy_ev, p_ev = ev.load_events(seq_path)
    raw_imu = np.load(seq_path / "dataset_imu.npy")
    cal_ts = raw_imu[:, 0]
    cal_w = ev.calibrated_gyro(raw_imu[:, 1:4])

    traj = info["meta"].get("full_trajectory") or []
    cam_ts, cam_q, cam_p = ev.pose_series(traj, "cam")

    series_dict = ev.find_series(info) if hasattr(ev, "find_series") else ec_find_series(info)
    frames_list = series_dict.get("meta/frames") or series_dict.get("frames") or []
    fts = {int(d.get("fid", d.get("frame_id", i))): float(d["ts"]) for i, d in enumerate(frames_list) if "ts" in d}

    _, tracks = ev.load_tracks(info)
    kin = ev.kinematics(tracks, "R*S", 0.025)
    z_mask = np.load(seq_path / "dataset_mask.npz")

    H, W = 480, 640
    tot_tp1, tot_fp1, tot_fn1 = 0, 0, 0
    tot_tp2, tot_fp2, tot_fn2 = 0, 0, 0
    tot_events = 0
    frames_evaluated = 0

    valid_fids = sorted([f for f in fts if f"mask_{f:010d}" in z_mask])
    bench_fids = [f for f in valid_fids if 15 <= f <= 120]
    print(f"Ewaluacja na {len(bench_fids)} klatkach sekwencji {seq_path.name}...")

    for fid in bench_fids:
        tsf = fts[fid]
        if tsf < 0.1 or tsf > t_ev[-1] - 0.1:
            continue

        idx = np.searchsorted(cam_ts, tsf)
        if idx < 2 or idx >= len(cam_ts) - 2:
            continue
        dt = cam_ts[idx+2] - cam_ts[idx-2]
        v_world = (cam_p[idx+2] - cam_p[idx-2]) / dt
        q_wc = cam_q[idx:idx+1]
        q_cw = qconj(q_wc)
        v_cam = qrot(q_cw, v_world[None, :])[0]

        win = ev.prep_window(tsf, 0.040, t_ev, xy_ev, intr)
        xs, ys = ev.warp(win, cal_ts, cal_w, -0.0309)
        te = win["te"]

        excl, fmask = ev.regions_for_frame(z_mask, fid, tsf, kin, slow_v=0.05, fast_v=0.3, w_thr=1.0, dilate_px=6)
        pix_x = np.clip(win["pix"][:, 0].astype(int), 0, 639)
        pix_y = np.clip(win["pix"][:, 1].astype(int), 0, 479)
        posm = fmask[pix_y, pix_x]
        negm = (~excl[pix_y, pix_x])

        # Detector v1 (tau=0.3, blur=3)
        rm = ev.residual_map(xs, ys, te, tsf, 3)
        if rm is None:
            continue
        D, S = rm
        m = (D >= 0.3) & (S >= 1.0 * S.mean())
        g = m.reshape(H // 4, 4, W // 4, 4).mean(axis=(1, 3)) >= 0.3
        from evimo2.detector import label_cc, box_blur
        lab, areas = label_cc(g)

        cx = np.clip(np.floor(xs / 4).astype(int), 0, W // 4 - 1)
        cy = np.clip(np.floor(ys / 4).astype(int), 0, H // 4 - 1)
        ev_lab = lab[cy, cx]

        # Early / Late images for patch cross-correlation
        m_e = te < tsf
        m_l = te >= tsf
        E = np.zeros((H, W), np.float32)
        L = np.zeros((H, W), np.float32)
        np.add.at(E, (np.clip(ys[m_e].astype(int), 0, H - 1), np.clip(xs[m_e].astype(int), 0, W - 1)), 1)
        np.add.at(L, (np.clip(ys[m_l].astype(int), 0, H - 1), np.clip(xs[m_l].astype(int), 0, W - 1)), 1)
        Es = box_blur(E, 3)
        Ls = box_blur(L, 3)

        clusters = []
        for cid in range(1, len(areas)):
            if areas[cid] < 8:
                continue
            cy_i, cx_i = np.where(lab == cid)
            x0 = max(0, int(cx_i.min() * 4 - 4))
            x1 = min(W, int((cx_i.max() + 1) * 4 + 4))
            y0 = max(0, int(cy_i.min() * 4 - 4))
            y1 = min(H, int((cy_i.max() + 1) * 4 + 4))
            clusters.append({
                "cid": cid,
                "bbox": (x0, x1, y0, y1),
                "area": areas[cid]
            })

        # Run epipolar filter v2
        kept_clusters = ev.filter_clusters_epipolar(
            clusters, Es, Ls, v_cam, intr,
            min_speed=args.min_speed, angle_thresh_deg=args.angle_thresh, max_shift=8
        )
        kept_cids = {c["cid"] for c in kept_clusters}

        pred_v1 = (ev_lab > 0)
        pred_v2 = np.isin(ev_lab, list(kept_cids))

        tot_tp1 += int((pred_v1 & posm).sum())
        tot_fp1 += int((pred_v1 & negm).sum())
        tot_fn1 += int(((~pred_v1) & posm).sum())

        tot_tp2 += int((pred_v2 & posm).sum())
        tot_fp2 += int((pred_v2 & negm).sum())
        tot_fn2 += int(((~pred_v2) & posm).sum())

        tot_events += len(xs)
        frames_evaluated += 1

    prec1 = tot_tp1 / (tot_tp1 + tot_fp1 + 1e-6) * 100
    rec1 = tot_tp1 / (tot_tp1 + tot_fn1 + 1e-6) * 100
    f1_1 = 2 * prec1 * rec1 / (prec1 + rec1 + 1e-6)

    prec2 = tot_tp2 / (tot_tp2 + tot_fp2 + 1e-6) * 100
    rec2 = tot_tp2 / (tot_tp2 + tot_fn2 + 1e-6) * 100
    f1_2 = 2 * prec2 * rec2 / (prec2 + rec2 + 1e-6)
    fp_ratio = tot_fp1 / max(1, tot_fp2)

    print("\n" + "=" * 70)
    print(f"PODSUMOWANIE WYNIKÓW {seq_path.name} ({frames_evaluated} klatek, {tot_events:,} zdarzeń)")
    print("=" * 70)
    print(f"{'Metryka':<25s} | {'Detektor v1 (tylko rotcomp)':<25s} | {'Detektor v2 (rotcomp + FOE)':<25s}")
    print("-" * 75)
    print(f"{'True Positives (TP)':<25s} | {tot_tp1:<25d} | {tot_tp2:<25d}")
    print(f"{'False Positives (FP)':<25s} | {tot_fp1:<25d} | {tot_fp2:<25d}")
    print(f"{'False Negatives (FN)':<25s} | {tot_fn1:<25d} | {tot_fn2:<25d}")
    print(f"{'Czułość (Recall)':<25s} | {rec1:<24.2f}% | {rec2:<24.2f}%")
    print(f"{'Precyzja (Precision)':<25s} | {prec1:<24.2f}% | {prec2:<24.2f}%")
    print(f"{'Miara F1':<25s} | {f1_1:<24.2f}% | {f1_2:<24.2f}%")
    print("-" * 75)
    print(f"Wskaźnik redukcji fałszywych alarmów (FP): {fp_ratio:.2f}x (spadek z {tot_fp1:,} do {tot_fp2:,})")
    print(f"Wzrost precyzji: +{prec2 - prec1:.2f} punktów procentowych")
    print("=" * 70)

    # Save to CSV
    csv_file = out_dir / "scene13_01_metrics.csv"
    with open(csv_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["model", "sequence", "frames", "total_events", "tp", "fp", "fn", "recall_pct", "precision_pct", "f1_pct", "fp_reduction_ratio"])
        writer.writerow(["v1_rotcomp", seq_path.name, frames_evaluated, tot_events, tot_tp1, tot_fp1, tot_fn1, f"{rec1:.2f}", f"{prec1:.2f}", f"{f1_1:.2f}", "1.00"])
        writer.writerow(["v2_epipolar_foe", seq_path.name, frames_evaluated, tot_events, tot_tp2, tot_fp2, tot_fn2, f"{rec2:.2f}", f"{prec2:.2f}", f"{f1_2:.2f}", f"{fp_ratio:.2f}"])
    print(f"Zapisano metryki do: {csv_file}")


def ec_find_series(info):
    """Pomocnicze wyszukiwanie serii klatek."""
    for key in ("meta/frames", "frames", "meta"):
        if key in info and isinstance(info[key], (list, dict)):
            val = info[key]
            if isinstance(val, list):
                return {"frames": val}
            if isinstance(val, dict) and "frames" in val:
                return val
    return {}


if __name__ == "__main__":
    main()
