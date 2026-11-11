#!/usr/bin/env python3
"""Generuje dataset_imu.npy i dataset_intrinsics.json w katalogach sekwencji dla szybkiego odczytu przez C++."""
import json
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import evimo2_common as ec

def prepare_seq(seq_dir: Path):
    imu_path = seq_dir / "dataset_imu.npy"
    intr_path = seq_dir / "dataset_intrinsics.json"
    
    if imu_path.exists() and intr_path.exists():
        return
        
    info = ec.load_info(seq_dir)
    if info is None:
        return
        
    # 1. IMU
    imus = [(p, l) for p, l in ec.find_series(info).items() if ec.is_imu_series(l)]
    imus = [x for x in imus if "left" in x[0]] or imus
    if imus:
        ts, gyro, acc = ec.imu_arrays(imus[0][1])
        # Zapisujemy tablicę (N, 7): ts, gx, gy, gz, ax, ay, az
        imu_mat = np.column_stack([ts, gyro, acc])
        np.save(imu_path, imu_mat)
        
    # 2. Intrinsics
    intr = ec.find_dict_with(info, ["fx", "fy", "cx", "cy"])
    if intr:
        d = {k: float(intr.get(k, 0.0)) for k in ["fx", "fy", "cx", "cy", "k1", "k2", "k3", "p1", "p2"]}
        intr_path.write_text(json.dumps(d, indent=2))
        
    print(f"Przygotowano dane C++ dla: {seq_dir.name}")

def main():
    root = Path("data/raw")
    if not root.exists():
        sys.exit("Brak katalogu data/raw")
        
    for p in sorted(root.glob("*/*_000000")):
        if p.is_dir():
            prepare_seq(p)

if __name__ == "__main__":
    main()
