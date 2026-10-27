#!/usr/bin/env python3
"""Zwięzła inspekcja katalogu sekwencji EVIMO2.

Użycie:  python3 inspect_evimo2.py <katalog_sekwencji> [--depth 3] [--events-only]
"""
import argparse
import os
import numpy as np
from evimo2_common import (describe, find_dict_with, find_series, load_events,
                           load_info)

CHUNK = 4_000_000


def hdr(title):
    print(f"\n{'=' * 6} {title} {'=' * 6}")


def list_files(d):
    hdr("PLIKI")
    for f in sorted(os.listdir(d)):
        p = os.path.join(d, f)
        if os.path.isfile(p):
            print(f"  {f:28s} {os.path.getsize(p) / 2**20:9.2f} MB")


def inspect_events(d, res):
    hdr("ZDARZENIA")
    t, xy, p = load_events(d)
    n = len(t)
    print(f"  N = {n:,} | t: {t.dtype}{t.shape} xy: {xy.dtype}{xy.shape} p: {p.dtype}{p.shape}")
    if n == 0:
        return
    t0, t1 = float(t[0]), float(t[-1])
    dur = t1 - t0
    nb = int(dur * 100) + 2                       # histogram: kubełki 10 ms
    hist = np.zeros(nb, dtype=np.int64)
    bad = 0; max_gap = 0.0; prev = None
    xmin = ymin = 10**9; xmax = ymax = -1
    on = 0; pvals = set()
    for i in range(0, n, CHUNK):
        c = np.asarray(t[i:i + CHUNK], dtype=np.float64)
        if prev is not None:
            bad += int(c[0] < prev); max_gap = max(max_gap, c[0] - prev)
        d = np.diff(c)
        bad += int((d < 0).sum())
        if d.size:
            max_gap = max(max_gap, float(d.max()))
        prev = c[-1]
        hist += np.bincount(((c - t0) * 100).astype(np.int64), minlength=nb)[:nb]
        q = np.asarray(xy[i:i + CHUNK])
        xmin = min(xmin, int(q[:, 0].min())); xmax = max(xmax, int(q[:, 0].max()))
        ymin = min(ymin, int(q[:, 1].min())); ymax = max(ymax, int(q[:, 1].max()))
        pp = np.asarray(p[i:i + CHUNK])
        on += int((pp != 0).sum()); pvals |= set(np.unique(pp).tolist())
    print(f"  czas: {t0:.6f} .. {t1:.6f} s (czas trwania {dur:.3f} s)")
    if t.dtype == np.float32:
        print(f"  UWAGA float32: rozdzielczość czasu przy końcu nagrania ~{np.spacing(np.float32(t1)) * 1e6:.2f} us")
    print(f"  niemonotoniczne: {bad} | najdłuższa przerwa: {max_gap * 1e3:.3f} ms")
    print(f"  średnio {n / max(dur, 1e-9) / 1e6:.3f} Mev/s | szczyt (10 ms): {hist.max() * 100 / 1e6:.3f} Mev/s")
    print(f"  x: {xmin}..{xmax} | y: {ymin}..{ymax}", end="")
    if res:
        ok = xmax < res[0] and ymax < res[1] and xmin >= 0 and ymin >= 0
        print(f"  (rozdzielczość z meta {res[0]}x{res[1]}: {'OK' if ok else 'POZA ZAKRESEM!'})")
    else:
        print()
    print(f"  polaryzacja: wartości {sorted(pvals)}, udział ON(!=0) = {on / n:.3f}")


def inspect_info(d, depth):
    hdr("dataset_info.npz")
    info = load_info(d)
    if info is None:
        print("  brak pliku"); return None, None
    for k, v in info.items():
        describe(v, k, 0, depth)
    series = find_series(info)
    if series:
        print("\n  Znalezione szeregi czasowe (listy ze 'ts'):")
        for path, lst in series.items():
            ts = np.array([x["ts"] for x in lst], dtype=float)
            dt = np.diff(ts)
            rate = 1 / np.median(dt) if dt.size else float("nan")
            print(f"   {path}: n={len(lst)} ts {ts[0]:.4f}..{ts[-1]:.4f} ~{rate:.1f} Hz")
    intr = find_dict_with(info, ["fx", "fy", "cx", "cy"])
    res = None
    if intr:
        print("\n  Intrinsics:", {k: (round(v, 5) if isinstance(v, float) else v) for k, v in intr.items()})
        if "res_x" in intr:
            res = (int(intr["res_x"]), int(intr["res_y"]))
    return info, res


def inspect_gt(d, name, label):
    path = os.path.join(d, name)
    hdr(name)
    if not os.path.isfile(path):
        print("  brak pliku"); return
    with np.load(path, allow_pickle=True) as z:
        keys = z.files
        print(f"  {len(keys)} wpisów ({label}); przykładowe klucze: {keys[:3]}")
        if not keys:
            print("  (pusty)"); return
        ids = sorted(int(k.split("_")[-1]) for k in keys if k.split("_")[-1].isdigit())
        if ids:
            print(f"  numery klatek: {ids[0]} .. {ids[-1]}")
        a = z[keys[0]]
        if a.dtype == object:
            print(f"  zawartość: {short_obj(a)}"); return
        print(f"  {keys[0]}: shape={a.shape} dtype={a.dtype} unikalne={np.unique(a)[:8].tolist()}")


def short_obj(a):
    s = repr(a.item() if a.shape == () else a)
    return s[:100]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("seq_dir")
    ap.add_argument("--depth", type=int, default=3, help="głębokość drzewa meta")
    ap.add_argument("--events-only", action="store_true")
    a = ap.parse_args()

    print(f"=== INSPEKCJA EVIMO2: {os.path.abspath(a.seq_dir)} ===")
    list_files(a.seq_dir)
    res = None
    if not a.events_only:
        _, res = inspect_info(a.seq_dir, a.depth)
    inspect_events(a.seq_dir, res)
    if not a.events_only:
        inspect_gt(a.seq_dir, "dataset_mask.npz", "maski = id obiektu * 1000")
        inspect_gt(a.seq_dir, "dataset_depth.npz", "głębia uint16")
        cl = os.path.join(a.seq_dir, "dataset_classical.npz")
        if os.path.isfile(cl):
            hdr("dataset_classical.npz")
            with np.load(cl, allow_pickle=True) as z:
                print(f"  klucze: {z.files}  (['empty'] = brak klatek klasycznych)")


if __name__ == "__main__":
    main()
