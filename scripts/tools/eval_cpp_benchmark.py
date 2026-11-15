#!/usr/bin/env python3
"""Automated C++ benchmark on the 7 selected EVIMO2 sequences."""
import os
import subprocess
import json
import csv
from pathlib import Path

SEQUENCES = [
    ("scene6_dyn_train_03_000000", "data/npz_left_camera_imo/left_camera/imo/train/scene6_dyn_train_03_000000", "Spoczynek (czyste IMO)", 14.3),
    ("scene14_dyn_test_05_000000", "data/npz_left_camera_imo/left_camera/imo/eval/scene14_dyn_test_05_000000", "Czysta rotacja (wzorzec)", -15.6),
    ("scene14_dyn_test_04_000000", "data/npz_left_camera_imo/left_camera/imo/eval/scene14_dyn_test_04_000000", "Czysta rotacja (mały ruch)", -15.7),
    ("scene13_dyn_test_02_000000", "data/npz_left_camera_imo/left_camera/imo/train/scene13_dyn_test_02_000000", "Szybka rotacja (stress-test FPR)", -21.3),
    ("scene13_dyn_test_01_000000", "data/npz_left_camera_imo/left_camera/imo/train/scene13_dyn_test_01_000000", "Ruch złożony (silna paralaksa)", -30.9),
    ("scene15_dyn_test_02_000000", "data/npz_left_camera_imo/left_camera/imo/eval/scene15_dyn_test_02_000000", "Ruch złożony (duży udział IMO)", -5.6),
    ("scene15_dyn_test_05_000000", "data/npz_left_camera_imo/left_camera/imo/eval/scene15_dyn_test_05_000000", "Ruch złożony (wysokie IMO)", -5.6),
]

OUT_DIR = Path("results/v1_rotcomp/benchmark")
OUT_DIR.mkdir(parents=True, exist_ok=True)

CPP_BIN = "./cpp/bin/imo_detect"

results = []

for seq_name, seq_path, cat, lag_ms in SEQUENCES:
    print(f"=== Running C++ Benchmark for {seq_name} ===")
    json_path = OUT_DIR / f"cpp_bench_{seq_name}.json"
    cmd = [
        CPP_BIN,
        seq_path,
        "--window-ms", "40",
        "--step-ms", "40",
        "--max-windows", "40",
        "--warmup", "2",
        "--blur", "3",
        "--tau", "0.5",
        "--lag-ms", str(lag_ms),
        "--json", str(json_path)
    ]
    ret = subprocess.run(cmd, capture_output=True, text=True)
    if ret.returncode != 0:
        print(f"ERROR running {seq_name}: {ret.stderr}")
        continue
    
    with open(json_path) as jf:
        d = json.load(jf)
    d["category"] = cat
    results.append(d)
    print(f"  Processed {d['processed_windows']} windows, {d['processed_events']} events, "
          f"Throughput: {d['throughput_mev']:.2f} Mev/s, RT-factor: {d['rt_factor']:.2f}x, "
          f"Avg window: {d['avg_window_ms']:.2f} ms")

# Save combined json
with open(OUT_DIR / "benchmark_cpp_7seqs.json", "w") as f:
    json.dump(results, f, indent=2)

# Save combined CSV
csv_cols = [
    "sequence", "category", "processed_windows", "processed_events",
    "total_time_ms", "avg_window_ms", "avg_undist_ms", "avg_warp_ms",
    "avg_detect_ms", "throughput_mev", "rt_factor"
]

with open(OUT_DIR / "benchmark_cpp_7seqs.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=csv_cols)
    writer.writeheader()
    for r in results:
        writer.writerow({k: r[k] for k in csv_cols})

print(f"\nWrote results to {OUT_DIR / 'benchmark_cpp_7seqs.csv'}")
