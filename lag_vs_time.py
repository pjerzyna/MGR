#!/usr/bin/env python3
"""Opóźnienie IMU względem GT (z gyro_vs_gt.json) w funkcji czasu nagrania (ros_time_offset).

Użycie:
    python3 lag_vs_time.py <katalog_z_sekwencjami> [--out-root DIR] [--csv PLIK] [--max-res 15]

  <katalog_z_sekwencjami>  np. data/npz_left_camera_imo/left_camera/imo (przeszukiwany rekurencyjnie)
  --out-root               gdzie leżą wyniki gyro_vs_gt.json (domyślnie $EVIMO2_OUT albo ./out)

Sekwencje są sortowane po czasie nagrania i dzielone na sesje (przerwa > 1 doba = nowa sesja).
W każdej sesji dopasowywana jest prosta d(t) z sekwencji o residuum < --max-res [%]; nachylenie
w ppm to dryf zegara IMU względem GT (11 ppm = 11 us/s).
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

import evimo2_common as ec

SESSION_GAP_S = 86400.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_root")
    ap.add_argument("--out-root", default=None)
    ap.add_argument("--csv", default=None, help="zapisz tabelę do pliku CSV")
    ap.add_argument("--max-res", type=float, default=15.0, help="próg residuum [%%] dla dopasowania prostej")
    a = ap.parse_args()

    rows, missing = [], []
    for info in sorted(Path(a.data_root).rglob("dataset_info.npz")):
        seq = info.parent.name
        j = ec.find_result_file(seq, "gyro_vs_gt.json", a.out_root)
        if j is None:
            missing.append(seq)
            continue
        with np.load(info, allow_pickle=True) as z:
            ros = float(z["meta"].item()["meta"]["ros_time_offset"])
        d = json.load(open(j))
        rows.append(dict(sequence=seq, ros_time_offset=ros, lag_ms=d["lag_s"] * 1e3,
                         residual_pct=d["residual_rel"] * 100, convention=d.get("convention", "")))
    if not rows:
        roots = ", ".join(str(c) for c in ec.result_dir_candidates("<sekwencja>", a.out_root))
        sys.exit(f"Nie znaleziono żadnego gyro_vs_gt.json (szukano w {roots}).")
    rows.sort(key=lambda r: r["ros_time_offset"])

    session = 0
    for i, r in enumerate(rows):
        if i and r["ros_time_offset"] - rows[i - 1]["ros_time_offset"] > SESSION_GAP_S:
            session += 1
        r["session"] = session
    for s in range(session + 1):
        t0 = min(r["ros_time_offset"] for r in rows if r["session"] == s)
        for r in rows:
            if r["session"] == s:
                r["t_in_session_s"] = r["ros_time_offset"] - t0

    print(f"{'sesja':>5} {'sekwencja':28s} {'t w sesji [s]':>14s} {'d [ms]':>8s} {'res [%]':>8s}")
    for r in rows:
        print(f"{r['session']:5d} {r['sequence']:28s} {r['t_in_session_s']:14.1f} {r['lag_ms']:8.1f} {r['residual_pct']:8.1f}")

    print(f"\nDryf d(t) per sesja (sekwencje z residuum < {a.max_res:.0f}%):")
    for s in range(session + 1):
        pts = [r for r in rows if r["session"] == s and r["residual_pct"] < a.max_res]
        n_all = sum(r["session"] == s for r in rows)
        if len(pts) < 3:
            print(f"  sesja {s}: za mało punktów ({len(pts)} z {n_all})")
            continue
        t = np.array([p["t_in_session_s"] for p in pts])
        d = np.array([p["lag_ms"] for p in pts])
        slope, icpt = np.polyfit(t, d, 1)
        rms = np.sqrt(np.mean((d - (slope * t + icpt)) ** 2))
        span = t.max() - t.min()
        print(f"  sesja {s}: {slope * 1e3:+.1f} ppm, rozrzut wokół prostej {rms:.2f} ms, "
              f"n={len(pts)} z {n_all}, rozpiętość {span:.0f} s")
    if missing:
        print(f"\nBez gyro_vs_gt.json ({len(missing)}): {', '.join(missing)}")

    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)
        print(f"Zapisano {a.csv}")


if __name__ == "__main__":
    main()