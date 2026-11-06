#!/usr/bin/env python3
"""Jednostki głębi w dataset_depth.npz: porównanie głębi pod maską obiektu z głębią obiektu z jego pozy.

Użycie:
    python3 depth_check.py <katalog_sekwencji> [--cls R*S] [--max-frames 60]

Poza obiektu w układzie kamery (klasa R*S: T_kamera<-obiekt) daje głębokość początku układu obiektu
z_poza = t_z [m, jeśli pozy są w metrach]. Dla każdej klatki z maską i głębią bierzemy medianę głębi
w pikselach maski obiektu i liczymy stosunek mediana_głębi / z_poza. Stosunek ~1000 oznacza głębię
w milimetrach (a pozy w metrach), ~1 - metry, ~10000 - 0,1 mm. Głębia pod maską to powierzchnia obiektu,
a poza dotyczy początku układu obiektu, więc spodziewany rozrzut stosunku to kilkadziesiąt procent
(rozmiar obiektu), nie rzędy wielkości.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

import evimo2_common as ec
import object_motion as om


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("seq_dir")
    ap.add_argument("--cls", default="R*S", choices=om.CLASSES)
    ap.add_argument("--max-frames", type=int, default=60)
    a = ap.parse_args()
    seq = Path(a.seq_dir.rstrip("/"))
    info = ec.load_info(seq)
    if info is None:
        sys.exit("brak dataset_info.npz")
    _, tracks = om.load_tracks(info)
    fts = {int(f["id"]): float(f["ts"]) for f in info["meta"].get("frames", [])}
    mz, dz = np.load(seq / "dataset_mask.npz"), np.load(seq / "dataset_depth.npz")
    ids = sorted(int(k.split("_")[-1]) for k in mz.files if f"depth_{k.split('_')[-1]}" in dz.files and int(k.split("_")[-1]) in fts)
    if not ids:
        sys.exit("brak klatek z maską, głębią i znacznikiem czasu")
    ids = [ids[i] for i in sorted(set(np.linspace(0, len(ids) - 1, min(a.max_frames, len(ids))).astype(int).tolist()))]
    ratios, per_obj = [], {}
    for fid in ids:
        lab = mz[f"mask_{fid:010d}"] // 1000
        depth = dz[f"depth_{fid:010d}"].astype(float)
        for oid in np.unique(lab[lab > 0]):
            k = str(int(oid))
            if k not in tracks:
                continue
            t = tracks[k]
            _, _, tt = (None, None, None)
            # głębokość początku układu obiektu w układzie kamery (klasa R*S: t_s)
            z = float(np.interp(fts[fid], t["ts"], t["ts_"][:, 2])) if a.cls == "R*S" else float("nan")
            d = depth[lab == oid]
            d = d[d > 0]
            if len(d) < 50 or not np.isfinite(z) or z <= 0.05:
                continue
            r = float(np.median(d)) / z
            ratios.append(r)
            per_obj.setdefault(k, []).append((float(np.median(d)), z))
    if not ratios:
        sys.exit("brak obiektów z dodatnią głębią i pozą")
    print(f"Sekwencja: {seq.name} | klatek: {len(ids)} | par (obiekt, klatka): {len(ratios)}")
    print(f"\\n{'id':>4} {'klatek':>7} {'mediana głębi':>14} {'z z pozy [m]':>13} {'stosunek':>9}")
    for k, v in sorted(per_obj.items(), key=lambda x: int(x[0])):
        d, z = np.median([p[0] for p in v]), np.median([p[1] for p in v])
        print(f"{k:>4} {len(v):7d} {d:14.1f} {z:13.3f} {d / z:9.1f}")
    med = float(np.median(ratios))
    p10, p90 = np.percentile(ratios, [10, 90])
    print(f"\\nStosunek mediana_głębi / z_poza: mediana {med:.1f} (p10 {p10:.1f}, p90 {p90:.1f})")
    unit = "milimetry (poza w metrach)" if 700 < med < 1400 else "metry" if 0.7 < med < 1.4 else "0,1 mm" if 7000 < med < 14000 else None
    print("Wniosek: głębia w jednostkach: " + (unit if unit else "NIEJEDNOZNACZNE (stosunek poza oczekiwanymi wartościami 1 / 1000 / 10000)"))
    print("(Zakładam, że translacje w pozach są w metrach - jeśli stosunek wychodzi 1000, to głębia jest w mm.)")


if __name__ == "__main__":
    main()