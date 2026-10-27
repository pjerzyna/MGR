#!/usr/bin/env python3
"""Punkt A: spójność czasowa zdarzeń / IMU / klatek ground truth w sekwencji EVIMO2.

Użycie:  python3 check_sync.py <katalog_sekwencji> [--masks 5] [--no-xcorr]

Sprawdza: 
(1) szeregi czasowe i ich zakresy, 
(2) czy IMU i zdarzenia dzielą oś czasu,
(3) jakość próbkowania IMU i grawitację, 
(4) szacowany ruch w pikselach na okno,
(5) mapowanie numerów masek/głębi na klatki meta, 
(6) heurystyczną korelację aktywności zdarzeń z |omega| (oszacowanie opóźnienia).

WYWOŁANIE TEJ KOMENDY:
python3 check_sync.py npz_left_camera_imo/left_camera/imo/train/scene15_dyn_test_00_000000/
"""
import argparse
import os
import numpy as np
from evimo2_common import (find_dict_with, find_series, imu_arrays, is_imu_series,
                           load_events, load_info)

OK, BAD, UNK = "[OK] ", "[!!] ", "[?]  "


def overlap_frac(a0, a1, b0, b1):
    return max(0.0, min(a1, b1) - max(a0, b0)) / max(b1 - b0, 1e-12)


def event_rate(t, t0, bin_s):
    nb = int((float(t[-1]) - t0) / bin_s) + 1
    h = np.zeros(nb)
    for i in range(0, len(t), 4_000_000):
        c = np.asarray(t[i:i + 4_000_000], dtype=np.float64)
        h += np.bincount(((c - t0) / bin_s).astype(np.int64), minlength=nb)[:nb]
    return h


def xcorr_lag(rate, bin_s, t0, imu_ts, wmag, max_lag=0.2):
    centers = t0 + (np.arange(len(rate)) + 0.5) * bin_s
    best = (None, -2.0)
    for lag in np.arange(-max_lag, max_lag + bin_s / 2, bin_s):
        w = np.interp(centers + lag, imu_ts, wmag, left=np.nan, right=np.nan)
        m = ~np.isnan(w)
        if m.sum() < 50:
            continue
        r = np.corrcoef(rate[m], w[m])[0, 1]
        if np.isfinite(r) and r > best[1]:
            best = (lag, r)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("seq_dir")
    ap.add_argument("--masks", type=int, default=5, help="ile masek przeanalizować")
    ap.add_argument("--no-xcorr", action="store_true")
    a = ap.parse_args()
    d = a.seq_dir
    verdict = []

    info = load_info(d)
    t, xy, p = load_events(d)
    e0, e1 = float(t[0]), float(t[-1])
    print(f"Zdarzenia: {len(t):,} | {e0:.6f} .. {e1:.6f} s ({e1 - e0:.3f} s)")
    if info is None:
        print(BAD + "brak dataset_info.npz"); return

    intr = find_dict_with(info, ["fx", "fy", "cx", "cy"])
    ros_off = float(intr["ros_time_offset"]) if intr and "ros_time_offset" in intr else None
    if ros_off is not None:
        print(f"ros_time_offset (meta): {ros_off:.6f}")

    # (1) szeregi czasowe
    series = find_series(info)
    print("\n-- Szeregi czasowe w meta --")
    imu_path = None
    frames_path = None
    for path, lst in series.items():
        ts = np.array([x["ts"] for x in lst], dtype=float)
        dt = np.diff(ts)
        kind = "IMU" if is_imu_series(lst) else "inne"
        print(f"  {path} [{kind}] n={len(lst)} ts {ts[0]:.6f}..{ts[-1]:.6f}"
              f" ~{(1 / np.median(dt) if dt.size else float('nan')):.1f} Hz | klucze: {list(lst[0])[:8]}")
        if is_imu_series(lst) and imu_path is None:
            imu_path = path
        elif not is_imu_series(lst) and ("gt_frame" in lst[0] or "id" in lst[0]) and frames_path is None:
            frames_path = path
    if imu_path is None:
        verdict.append((BAD, "nie znaleziono szeregu IMU (lista ze 'ts', 'angular_velocity', 'linear_acceleration')"))
    if frames_path is None:
        verdict.append((UNK, "nie znaleziono listy klatek GT (z 'id'/'gt_frame')"))

    # (2) oś czasu IMU vs zdarzenia
    imu_ts = gyro = acc = None
    if imu_path:
        imu_ts, gyro, acc = imu_arrays(series[imu_path])
        i0, i1 = imu_ts[0], imu_ts[-1]
        print(f"\n-- IMU vs zdarzenia ({imu_path}) --")
        hyps = {"ts IMU bez zmian": 0.0}
        if ros_off:
            hyps["ts IMU - ros_time_offset"] = -ros_off
            hyps["ts IMU + ros_time_offset"] = +ros_off
        best_name, best_ov = None, -1
        for name, sh in hyps.items():
            ov = overlap_frac(i0 + sh, i1 + sh, e0, e1)
            print(f"  hipoteza '{name}': pokrycie zdarzeń przez IMU = {ov * 100:.1f}%")
            if ov > best_ov:
                best_name, best_ov = name, ov
        shift = hyps[best_name]
        if best_ov > 0.95:
            verdict.append((OK, f"IMU pokrywa {best_ov * 100:.0f}% czasu zdarzeń (hipoteza: {best_name})"))
        elif best_ov > 0:
            verdict.append((BAD, f"IMU pokrywa tylko {best_ov * 100:.0f}% czasu zdarzeń (najlepsza hipoteza: {best_name})"))
        else:
            verdict.append((BAD, "oś czasu IMU nie przecina się ze zdarzeniami - wymagane przesunięcie"))
        imu_ts = imu_ts + shift

        # (3) jakość IMU
        dt = np.diff(imu_ts)
        wm = np.linalg.norm(gyro, axis=1)
        am = np.linalg.norm(acc, axis=1)
        print("\n-- Jakość IMU --")
        print(f"  dt: mediana {np.median(dt) * 1e3:.3f} ms, std {dt.std() * 1e3:.3f} ms, max {dt.max() * 1e3:.3f} ms")
        print(f"  |omega| [rad/s]: średnia {wm.mean():.3f}, p95 {np.percentile(wm, 95):.3f}, max {wm.max():.3f}")
        g = acc.mean(axis=0)
        print(f"  |a| średnia {am.mean():.3f} m/s^2 (oczek. ~9.81 lub ~1 g) | wektor średni {np.round(g, 2).tolist()}")
        if (dt <= 0).any():
            verdict.append((BAD, "ts IMU nie jest ściśle rosnący"))
        if dt.max() > 3 * np.median(dt):
            verdict.append((BAD, f"luki w IMU (max dt = {dt.max() * 1e3:.2f} ms)"))
        else:
            verdict.append((OK, "próbkowanie IMU równomierne, bez luk"))
        if intr:
            fx = float(intr["fx"])
            px = wm.max() * fx * 0.010
            print(f"  maks. przesunięcie od rotacji w 10 ms: ~{px:.1f} px (fx={fx:.1f})")
            verdict.append((OK if px > 1 else UNK,
                            f"rotacja daje maks. ~{px:.1f} px/10 ms" + ("" if px > 1 else " - mało ruchu, sekwencja słaba do testu kompensacji")))

        # (6) xcorr
        if not a.no_xcorr and overlap_frac(imu_ts[0], imu_ts[-1], e0, e1) * (e1 - e0) > 0.5:
            bin_s = 0.005
            rate = event_rate(t, e0, bin_s)
            lag, r = xcorr_lag(rate, bin_s, e0, imu_ts, wm)
            print("\n-- Heurystyczna korelacja rate zdarzeń vs |omega| --")
            if lag is None:
                print("  za mało danych")
            else:
                print(f"  najlepsze opóźnienie: {lag * 1e3:+.0f} ms, korelacja {r:.2f}")
                verdict.append((UNK, f"xcorr: lag {lag * 1e3:+.0f} ms, r={r:.2f} (heurystyka; sensowna tylko przy dominującym ruchu kamery)"))

    # (4/5) klatki GT
    gt_ids = None
    for fn in ("dataset_mask.npz", "dataset_depth.npz"):
        pth = os.path.join(d, fn)
        if os.path.isfile(pth):
            with np.load(pth, allow_pickle=True) as z:
                ids = sorted(int(k.split("_")[-1]) for k in z.files if k.split("_")[-1].isdigit())
            print(f"\n{fn}: {len(ids)} klatek, numery {ids[:3]}...{ids[-3:]}" if ids else f"\n{fn}: pusty")
            gt_ids = gt_ids or ids
    if frames_path and gt_ids:
        frames = series[frames_path]
        fts = np.array([f["ts"] for f in frames], dtype=float)
        print(f"\n-- Klatki GT ({frames_path}) --")
        print(f"  {len(frames)} klatek, ts {fts[0]:.4f}..{fts[-1]:.4f}, pola: {list(frames[0])[:10]}")
        idset = set(gt_ids)
        matched = None
        for field, vals in [(k, [f.get(k) for f in frames]) for k in frames[0]
                            if isinstance(frames[0][k], (int, np.integer))]:
            if idset <= set(vals):
                matched = field; break
        if matched is None and idset <= set(range(len(frames))):
            matched = "(indeks na liście)"
        if matched:
            verdict.append((OK, f"numery masek/głębi mapują się na klatki meta przez pole {matched}"))
        else:
            verdict.append((BAD, "nie udało się zmapować numerów masek/głębi na klatki meta"))
        shift_f = 0.0
        inside = np.mean((fts >= e0) & (fts <= e1))
        print(f"  klatki GT w zakresie czasu zdarzeń: {inside * 100:.0f}%")
        verdict.append((OK if inside > 0.95 else BAD, f"{inside * 100:.0f}% klatek GT w zakresie czasu zdarzeń"))

    # maski
    mpath = os.path.join(d, "dataset_mask.npz")
    if os.path.isfile(mpath):
        print("\n-- Maski (próbka) --")
        with np.load(mpath, allow_pickle=True) as z:
            for k in z.files[: a.masks]:
                m = z[k]
                u, c = np.unique(m, return_counts=True)
                desc = ", ".join(f"{v // 1000 if v % 1000 == 0 else v}:{cnt / m.size * 100:.1f}%" for v, cnt in zip(u, c))
                print(f"  {k}: obiekty(id:udział) {desc}")
                if (u[u > 0] % 1000 != 0).any():
                    verdict.append((BAD, f"{k}: wartości maski niebędące wielokrotnością 1000"))

    print("\n===== PODSUMOWANIE =====")
    for tag, msg in verdict:
        print(tag + msg)


if __name__ == "__main__":
    main()
