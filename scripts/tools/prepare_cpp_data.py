#!/usr/bin/env python3
"""Generuje dataset_imu.npy, dataset_intrinsics.json i dataset_cam_velocity.npy w katalogach sekwencji dla C++."""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path("python").resolve()))

import numpy as np
import evimo2 as ev
from evimo2.geometry import qrot, qconj


def prepare_seq(seq_dir: Path):
    imu_path = seq_dir / "dataset_imu.npy"
    intr_path = seq_dir / "dataset_intrinsics.json"
    vel_path = seq_dir / "dataset_cam_velocity.npy"

    if imu_path.exists() and intr_path.exists() and vel_path.exists():
        return

    info = ev.load_info(seq_dir)
    if info is None:
        return

    # 1. IMU
    if not imu_path.exists():
        imus = [(p, l) for p, l in ev.find_series(info).items() if ev.is_imu_series(l)]
        imus = [x for x in imus if "left" in x[0]] or imus
        if imus:
            ts, gyro, acc = ev.imu_arrays(imus[0][1])
            imu_mat = np.column_stack([ts, gyro, acc])
            np.save(imu_path, imu_mat)

    # 2. Intrinsics
    if not intr_path.exists():
        intr = ev.find_dict_with(info, ["fx", "fy", "cx", "cy"])
        if intr:
            d = {k: float(intr.get(k, 0.0)) for k in ["fx", "fy", "cx", "cy", "k1", "k2", "k3", "p1", "p2"]}
            intr_path.write_text(json.dumps(d, indent=2))

    # 3. Camera Linear Velocity in camera frame (N, 4): ts, vx_cam, vy_cam, vz_cam
    if not vel_path.exists() and "meta" in info and "full_trajectory" in info["meta"]:
        traj = info["meta"].get("full_trajectory") or []
        res = ev.pose_series(traj, "cam")
        if res is not None:
            cam_ts, cam_q, cam_p = res
            N = len(cam_ts)
            v_cam_arr = np.zeros((N, 4), dtype=np.float64)
            v_cam_arr[:, 0] = cam_ts
            for i in range(N):
                i0 = max(0, i - 2)
                i1 = min(N - 1, i + 2)
                dt = cam_ts[i1] - cam_ts[i0]
                if dt > 1e-6:
                    v_world = (cam_p[i1] - cam_p[i0]) / dt
                    q_wc = cam_q[i:i+1]
                    q_cw = qconj(q_wc)
                    v_cam_arr[i, 1:4] = qrot(q_cw, v_world[None, :])[0]
            np.save(vel_path, v_cam_arr)

    print(f"Przygotowano dane C++ dla: {seq_dir.name}")


def main():
    root = Path("data/npz_left_camera_imo/left_camera/imo")
    if not root.exists():
        root = Path("data/raw")
    if not root.exists():
        sys.exit("Brak katalogu z danymi")

    for p in sorted(root.glob("*/*_000000")):
        if p.is_dir():
            prepare_seq(p)


if __name__ == "__main__":
    main()
