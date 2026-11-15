#!/usr/bin/env python3
"""Eksportuje referencyjne dane pośrednie jednego okna EVIMO2 do weryfikacji numerycznej C++."""
import json
import sys
from pathlib import Path
import numpy as np

# Dodajemy katalog nadrzędny
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import evimo2_common as ec
import rotcomp as rc
import evimo2.detector as ed

def main():
    seq_dir = Path("data/raw/eval/scene14_dyn_test_05_000000")
    if not seq_dir.exists():
        sys.exit(f"Brak katalogu {seq_dir}")

    info = ec.load_info(seq_dir)
    intr = ec.find_dict_with(info, ["fx", "fy", "cx", "cy"])
    
    imus = [(p, l) for p, l in ec.find_series(info).items() if ec.is_imu_series(l)]
    imus = [x for x in imus if "left" in x[0]] or imus
    cal_ts, gyro, _ = ec.imu_arrays(imus[0][1])
    cal_w = rc.calibrated_gyro(gyro)
    
    t_ev, xy_ev, p_ev = ec.load_events(seq_dir)
    fts = {int(f["id"]): float(f["ts"]) for f in info["meta"].get("frames", [])}
    
    mpath = seq_dir / "dataset_mask.npz"
    z = np.load(mpath)
    frames = sorted(fid for fid in (int(k.split("_")[-1]) for k in z.files) if fid in fts)
    
    wlen = 0.040  # 40 ms
    lo = max(cal_ts[0], float(t_ev[0])) + wlen / 2.0 + 0.06
    hi_t = min(cal_ts[-1], float(t_ev[-1])) - wlen / 2.0 - 0.06
    valid_frames = [fid for fid in frames if lo <= fts[fid] <= hi_t]
    if not valid_frames:
        sys.exit("Brak klatek w przedziale [lo, hi_t]")

    # Wybieramy klatkę w 1/3 sekwencji
    target_fid = valid_frames[len(valid_frames) // 3]
    tsf = fts[target_fid]
    
    print(f"Eksport okna wokół klatki {target_fid}: t_c = {tsf:.6f} s, wlen = {wlen * 1e3:.0f} ms")
    
    win = rc.prep_window(tsf, wlen, t_ev, xy_ev, intr)
    i0, i1 = np.searchsorted(t_ev, [tsf - wlen / 2.0, tsf + wlen / 2.0])
    p_win = p_ev[i0:i1]
    
    d = -0.0157  # standardowy lag
    M = rc.rot_matrices(win, cal_ts, cal_w, d)
    xs, ys = rc.warp(win, cal_ts, cal_w, d)
    
    # Detekcja
    early = win["te"] < tsf
    E = rc.iwe(xs[early], ys[early])
    L = rc.iwe(xs[~early], ys[~early])
    
    ne, nl = E.sum(), L.sum()
    E_norm, L_norm = E / ne, L / nl
    for _ in range(2):
        E_norm = ed.box_blur(E_norm, 3)
        L_norm = ed.box_blur(L_norm, 3)
    S = E_norm + L_norm
    D = np.abs(E_norm - L_norm) / (S + 1e-12)
    
    pred, n_clusters = ed.detect_imo(xs, ys, win["te"], tsf, blur=3, tau=0.5, support=1.0, min_cells=8)
    
    out_dir = Path("cpp/tests/data")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Zapisujemy pliki NPY
    np.save(out_dir / "t_ev.npy", win["te"])
    np.save(out_dir / "pix_ev.npy", win["pix"])
    np.save(out_dir / "p_ev.npy", p_win)
    np.save(out_dir / "bearing.npy", win["bearing"])
    np.save(out_dir / "grid.npy", win["grid"])
    np.save(out_dir / "rot_matrices.npy", M)
    np.save(out_dir / "warped_xs.npy", xs)
    np.save(out_dir / "warped_ys.npy", ys)
    np.save(out_dir / "iwe_E.npy", E)
    np.save(out_dir / "iwe_L.npy", L)
    np.save(out_dir / "D.npy", D)
    np.save(out_dir / "S.npy", S)
    np.save(out_dir / "pred_imo.npy", pred)
    
    # Wycinamy fragment IMU obejmujący okno (+- 0.1 s)
    imu_mask = (cal_ts >= tsf - 0.1) & (cal_ts <= tsf + 0.1)
    np.save(out_dir / "imu_ts.npy", cal_ts[imu_mask])
    np.save(out_dir / "imu_w.npy", cal_w[imu_mask])
    np.save(out_dir / "imu_raw_gyro.npy", gyro[imu_mask])
    
    meta = {
        "sequence": seq_dir.name,
        "sequence_short": "scene14_05",
        "frame_id": int(target_fid),
        "tc": tsf,
        "wlen": wlen,
        "d": d,
        "n_events": int(win["n"]),
        "n_clusters": int(n_clusters),
        "intrinsics": {k: float(intr[k]) for k in ["fx", "fy", "cx", "cy", "k1", "k2", "p1", "p2"]}
    }
    with open(out_dir / "meta.json", "w") as f:
        json.dump(meta, f, indent=2)
        
    print(f"Zapisano dane referencyjne dla {win['n']} zdarzeń do {out_dir}/")

if __name__ == "__main__":
    main()
