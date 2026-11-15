#!/usr/bin/env python3
"""Orchestrator for comprehensive evaluation across 7 EVIMO2 sequences."""
import os
import sys
import subprocess
import json
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

SEQUENCES = [
    ("scene6_dyn_train_03_000000", "data/npz_left_camera_imo/left_camera/imo/train/scene6_dyn_train_03_000000", "Spoczynek (czyste IMO)", 14.3),
    ("scene14_dyn_test_05_000000", "data/npz_left_camera_imo/left_camera/imo/eval/scene14_dyn_test_05_000000", "Czysta rotacja (wzorzec)", -15.6),
    ("scene14_dyn_test_04_000000", "data/npz_left_camera_imo/left_camera/imo/eval/scene14_dyn_test_04_000000", "Czysta rotacja (mały ruch)", -15.7),
    ("scene13_dyn_test_02_000000", "data/npz_left_camera_imo/left_camera/imo/train/scene13_dyn_test_02_000000", "Szybka rotacja (stress-test FPR)", -21.3),
    ("scene13_dyn_test_01_000000", "data/npz_left_camera_imo/left_camera/imo/train/scene13_dyn_test_01_000000", "Ruch złożony (silna paralaksa)", -30.9),
    ("scene15_dyn_test_02_000000", "data/npz_left_camera_imo/left_camera/imo/eval/scene15_dyn_test_02_000000", "Ruch złożony (duży udział IMO)", -5.6),
    ("scene15_dyn_test_05_000000", "data/npz_left_camera_imo/left_camera/imo/eval/scene15_dyn_test_05_000000", "Ruch złożony (wysokie IMO)", -5.6),
]

EVAL_DIR = Path("results/v1_rotcomp/benchmark")
EVAL_DIR.mkdir(parents=True, exist_ok=True)

PYTHON_EXE = sys.executable

def run_grid_eval():
    for seq_name, seq_path, cat, lag_ms in SEQUENCES:
        out_sub = EVAL_DIR / f"grid_{seq_name}"
        res_file = out_sub / "imo_detect.json"
        if res_file.exists():
            print(f"[{seq_name}] Already computed ({res_file}), skipping run.")
            continue
        print(f"[{seq_name}] Running parameter sweep grid...")
        cmd = [
            PYTHON_EXE, "scripts/tools/imo_detect.py",
            seq_path,
            "--grid",
            "--grid-blur", "2", "3", "5",
            "--grid-tau", "0.2", "0.3", "0.4", "0.5",
            "--max-windows", "40",
            "--out", str(out_sub)
        ]
        ret = subprocess.run(cmd, capture_output=True, text=True)
        if ret.returncode != 0:
            print(f"Error on {seq_name}:\n{ret.stderr}")
        else:
            print(f"[{seq_name}] Finished successfully.")

if __name__ == "__main__":
    run_grid_eval()
