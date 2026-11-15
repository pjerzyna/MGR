#!/usr/bin/env python3
"""Kompleksowy skrypt do automatycznej ewaluacji i benchmarku na 7 sekwencjach EVIMO2.

Użycie:
    python scripts/run_evaluation.py            # Uruchamia pełny potok end-to-end
    python scripts/run_evaluation.py --bench    # Tylko pomiary wydajności C++
    python scripts/run_evaluation.py --quality  # Tylko pomiary jakości względem Ground Truth
    python scripts/run_evaluation.py --plots    # Tylko generowanie wykresów i raportów z JSON-ów
    python scripts/run_evaluation.py --clean    # Wymusza przeliczenie wszystkiego od zera
"""
import sys
import os
import argparse
import subprocess
import json
import csv
from pathlib import Path
import numpy as np

SEQUENCES = [
    ("scene6_dyn_train_03_000000", "data/npz_left_camera_imo/left_camera/imo/train/scene6_dyn_train_03_000000", "Spoczynek (czyste IMO)", "Kategoria I", 14.3),
    ("scene14_dyn_test_05_000000", "data/npz_left_camera_imo/left_camera/imo/eval/scene14_dyn_test_05_000000", "Czysta rotacja (wzorzec)", "Kategoria II", -15.6),
    ("scene14_dyn_test_04_000000", "data/npz_left_camera_imo/left_camera/imo/eval/scene14_dyn_test_04_000000", "Czysta rotacja (mały ruch)", "Kategoria II", -15.7),
    ("scene13_dyn_test_02_000000", "data/npz_left_camera_imo/left_camera/imo/train/scene13_dyn_test_02_000000", "Szybka rotacja (stress-test FPR)", "Kategoria II", -21.3),
    ("scene13_dyn_test_01_000000", "data/npz_left_camera_imo/left_camera/imo/train/scene13_dyn_test_01_000000", "Ruch złożony (silna paralaksa)", "Kategoria III", -30.9),
    ("scene15_dyn_test_02_000000", "data/npz_left_camera_imo/left_camera/imo/eval/scene15_dyn_test_02_000000", "Ruch złożony (duży udział IMO)", "Kategoria III", -5.6),
    ("scene15_dyn_test_05_000000", "data/npz_left_camera_imo/left_camera/imo/eval/scene15_dyn_test_05_000000", "Ruch złożony (wysokie IMO)", "Kategoria III", -5.6),
]

EVAL_DIR = Path("results/summaries/evaluation")
EVAL_DIR.mkdir(parents=True, exist_ok=True)
PYTHON_EXE = sys.executable
CPP_BIN = Path("./cpp/bin/imo_detect")


def ensure_cpp_binary():
    print("[1/4] Sprawdzanie i kompilacja silnika C++...")
    ret = subprocess.run(["make", "-C", "cpp"], capture_output=True, text=True)
    if ret.returncode != 0:
        print(f"Błąd kompilacji C++:\n{ret.stderr}")
        sys.exit(1)
    if not CPP_BIN.exists():
        print(f"Błąd: Brak pliku binarnego {CPP_BIN}")
        sys.exit(1)
    print("      Silnik C++ gotowy.")


def run_cpp_benchmark(force=False):
    print("\n[2/4] Uruchamianie pomiarów wydajnościowych C++ na 7 sekwencjach...")
    csv_path = EVAL_DIR / "benchmark_cpp_7seqs.csv"
    if csv_path.exists() and not force:
        print(f"      Wyniki już istnieją ({csv_path}). Użyj --clean, aby wymusić.")
        return

    results = []
    for seq_name, seq_path, cat_desc, cat_group, lag_ms in SEQUENCES:
        json_path = EVAL_DIR / f"cpp_bench_{seq_name}.json"
        cmd = [
            str(CPP_BIN), seq_path,
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
            print(f"Błąd uruchomienia C++ dla {seq_name}: {ret.stderr}")
            continue

        with open(json_path) as jf:
            d = json.load(jf)
        d["category"] = cat_desc
        results.append(d)
        print(f"  -> {seq_name[:26]:26s} | {d['processed_windows']:2d} okien | {d['throughput_mev']:5.2f} Mev/s | RT: {d['rt_factor']:4.2f}x | {d['avg_window_ms']:5.2f} ms/okno")

    with open(EVAL_DIR / "benchmark_cpp_7seqs.json", "w") as f:
        json.dump(results, f, indent=2)

    csv_cols = [
        "sequence", "category", "processed_windows", "processed_events",
        "total_time_ms", "avg_window_ms", "avg_undist_ms", "avg_warp_ms",
        "avg_detect_ms", "throughput_mev", "rt_factor"
    ]
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_cols)
        writer.writeheader()
        for r in results:
            writer.writerow({k: r[k] for k in csv_cols})
    print(f"      Zapisano: {csv_path}")


def run_quality_evaluation(force=False):
    print("\n[3/4] Uruchamianie ewaluacji jakości detekcji (siatka parametrów blur x tau) względem GT...")
    for seq_name, seq_path, cat_desc, cat_group, lag_ms in SEQUENCES:
        out_sub = EVAL_DIR / f"grid_{seq_name}"
        res_file = out_sub / "imo_detect.json"
        if res_file.exists() and not force:
            print(f"      [{seq_name[:26]:26s}] Wyniki istnieją, pomijam.")
            continue

        print(f"      [{seq_name[:26]:26s}] Przeliczanie siatki blur{{2,3,5}} x tau{{0.2..0.5}}...")
        cmd = [
            PYTHON_EXE, "python/scripts/imo_detect.py",
            seq_path,
            "--grid",
            "--grid-blur", "2", "3", "5",
            "--grid-tau", "0.2", "0.3", "0.4", "0.5",
            "--max-windows", "40",
            "--out", str(out_sub)
        ]
        ret = subprocess.run(cmd, capture_output=True, text=True)
        if ret.returncode != 0:
            print(f"Błąd dla {seq_name}:\n{ret.stderr}")


def generate_reports():
    print("\n[4/4] Generowanie tabel CSV, wykresów PNG oraz raportu Markdown...")
    subprocess.run([PYTHON_EXE, "scripts/generate_eval_reports.py"], check=True)
    print("      Wszystkie raporty i wykresy zostały wygenerowane w results/summaries/evaluation/")


def main():
    parser = argparse.ArgumentParser(description="Uruchamia zbiorczą ewaluację C++ na całym zestawie 7 sekwencji EVIMO2.")
    parser.add_argument("--bench", action="store_true", help="Tylko benchmark wydajności C++")
    parser.add_argument("--quality", action="store_true", help="Tylko ewaluacja jakości Ground Truth")
    parser.add_argument("--plots", action="store_true", help="Tylko generowanie wykresów z istniejących wyników")
    parser.add_argument("--clean", action="store_true", help="Wymuś przeliczenie wszystkiego od zera")
    args = parser.parse_args()

    if args.plots:
        generate_reports()
        return

    ensure_cpp_binary()

    if args.bench:
        run_cpp_benchmark(force=args.clean)
        generate_reports()
        return

    if args.quality:
        run_quality_evaluation(force=args.clean)
        generate_reports()
        return

    # Default: Run everything end-to-end
    run_cpp_benchmark(force=args.clean)
    run_quality_evaluation(force=args.clean)
    generate_reports()
    print("\nGotowe! Kompletna ewaluacja znajduje się w: results/summaries/evaluation/")


if __name__ == "__main__":
    main()
