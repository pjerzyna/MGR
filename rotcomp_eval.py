#!/usr/bin/env python3
"""Ewaluacja kompensacji rotacji (rotcomp.py): błąd względem orientacji z GT i wrażliwość na opóźnienie.

Użycie:
    python3 rotcomp_eval.py <katalog_sekwencji> [--window-ms 30] [--max-windows 40] [--calib fixed|json]
                            [--lag-const-ms -15.7] [--slow-v 0.05] [--fast-v 0.3] [--out DIR]

Dla okien zdarzeń wokół klatek z maską porównuje kilka wariantów OPÓŹNIENIA IMU (parametr wejściowy
algorytmu; kalibracja R i bias - stałe zestawu z rotcomp.py):
  * oracle d_GT       - opóźnienie z gyro_vs_gt.json dla tej sekwencji (GÓRNA GRANICA, nie realistyczny wariant),
  * stałe d           - jedna wartość dla całego zbioru (mediana, domyślnie -15,7 ms),
  * zero              - brak korekty opóźnienia,
  * oracle +-5/10/25 ms - sztucznie wstrzyknięty błąd opóźnienia.
Metryki (GT służy WYŁĄCZNIE do oceny):
  * błąd kompensacji [px] = odległość między pozycją zdarzenia po kompensacji z IMU a pozycją po
    kompensacji orientacją z poz kamery GT (czysty obrót) - mierzy błąd IMU/kalibracji/opóźnienia, bez
    wpływu translacji i ruchu obiektów; podawane: mediana okien z p50 i p95 błędu w oknie,
  * zysk ostrości tła = wariancja obrazu zdarzeń z kompensacją / bez kompensacji, tylko dla zdarzeń
    spoza masek obiektów ruchomych (v > --slow-v lub omega > 1), po rozszerzeniu masek o 6 px,
  * zysk dla szybkich IMO (v > --fast-v) - kompensacja ego-ruchu nie powinna ostrzyć ruchomych obiektów.
Zapisuje <out>/rotcomp_eval.json.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

import evimo2_common as ec
import imo_labels as il
import object_motion as om
import rotcomp as rc
from gyro_vs_gt import unify_hemisphere
from survey_sequences import pose_series


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("seq_dir")
    ap.add_argument("--window-ms", type=float, default=30.0)
    ap.add_argument("--max-windows", type=int, default=40)
    ap.add_argument("--calib", choices=["fixed", "json"], default="fixed")
    ap.add_argument("--lag-const-ms", type=float, default=rc.LAG_CONST_S * 1e3)
    ap.add_argument("--slow-v", type=float, default=0.05)
    ap.add_argument("--fast-v", type=float, default=0.3)
    ap.add_argument("--w-thr", type=float, default=1.0)
    ap.add_argument("--dilate", type=int, default=6)
    ap.add_argument("--min-events", type=int, default=5000)
    ap.add_argument("--min-fast-events", type=int, default=2000)
    ap.add_argument("--cls", default="R*S", choices=om.CLASSES)
    ap.add_argument("--out", default=None)
    ap.add_argument("--out-root", default=None)
    a = ap.parse_args()

    seq = Path(a.seq_dir.rstrip("/"))
    info = ec.load_info(seq)
    if info is None:
        sys.exit("brak dataset_info.npz")
    intr = ec.find_dict_with(info, ["fx", "fy", "cx", "cy"])
    imus = [(p, l) for p, l in ec.find_series(info).items() if ec.is_imu_series(l)]
    imus = [x for x in imus if "left" in x[0]] or imus
    if not imus:
        sys.exit("brak IMU w meta")
    cal_ts, gyro, _ = ec.imu_arrays(imus[0][1])
    gj_path = ec.find_result_file(seq.name, "gyro_vs_gt.json", a.out_root)
    gj = json.load(open(gj_path)) if gj_path else None
    if a.calib == "json":
        if gj is None:
            sys.exit("--calib json wymaga gyro_vs_gt.json (uruchom gyro_vs_gt.py)")
        cal_w = rc.calibrated_gyro(gyro, np.array(gj["R_imu_to_cam"]), np.array(gj["gyro_bias_imu_rad_s"]))
    else:
        cal_w = rc.calibrated_gyro(gyro)
    d_oracle = float(gj["lag_s"]) if gj else None

    cam = pose_series(info["meta"].get("full_trajectory") or [], "cam")
    if cam is None:
        sys.exit("brak poz kamery w full_trajectory (potrzebne do oceny)")
    ts_gt, q_gt, _ = cam
    q_gt = unify_hemisphere(q_gt)
    _, tracks = om.load_tracks(info)
    kin = il.kinematics(tracks, a.cls, 0.025)
    t_ev, xy_ev, _ = ec.load_events(seq)
    fts = {int(f["id"]): float(f["ts"]) for f in info["meta"].get("frames", [])}

    W = a.window_ms * 1e-3
    lo = max(cal_ts[0], float(t_ev[0]), ts_gt[0]) + W / 2 + 0.06
    hi = min(cal_ts[-1], float(t_ev[-1]), ts_gt[-1]) - W / 2 - 0.06
    mpath = seq / "dataset_mask.npz"
    if not mpath.is_file():
        sys.exit("brak dataset_mask.npz")
    z = np.load(mpath)
    frames = sorted((fid, fts[fid]) for fid in (int(k.split("_")[-1]) for k in z.files) if fid in fts and lo <= fts[fid] <= hi)
    if not frames:
        sys.exit("brak klatek z maską w dopuszczalnym przedziale czasu")
    sel = sorted(set(np.linspace(0, len(frames) - 1, min(a.max_windows, len(frames))).astype(int).tolist()))
    frames = [frames[i] for i in sel]

    modes = []
    if d_oracle is not None:
        modes.append((f"oracle d_GT = {d_oracle * 1e3:+.1f} ms", d_oracle))
    modes.append((f"stałe d = {a.lag_const_ms:+.1f} ms", a.lag_const_ms * 1e-3))
    modes.append(("zero d = 0", 0.0))
    if d_oracle is not None:
        for e in (5, 10, 25):
            modes += [(f"oracle +{e} ms", d_oracle + e * 1e-3), (f"oracle -{e} ms", d_oracle - e * 1e-3)]
    res = {name: dict(p50=[], p95=[], gB=[], gI=[]) for name, _ in modes}
    raw = dict(p50=[], p95=[]); gt_ref = dict(gB=[], gI=[]); n_ev, n_fast = [], 0

    for fid, tsf in frames:
        win = rc.prep_window(tsf, W, t_ev, xy_ev, intr)
        if win["n"] < a.min_events:
            continue
        lab = z[f"mask_{fid:010d}"] // 1000
        slow, fast = [], []
        for oid in np.unique(lab[lab > 0]):
            k = str(int(oid))
            if k not in kin:
                continue
            v = float(np.interp(tsf, kin[k]["ts"], kin[k]["v"])); w = float(np.interp(tsf, kin[k]["ts"], kin[k]["w"]))
            if v > a.slow_v or w > a.w_thr:
                slow.append(int(oid))
            if v > a.fast_v:
                fast.append(int(oid))
        excl = il.dilate(np.isin(lab, slow), a.dilate) if slow else np.zeros(lab.shape, bool)
        fmask = il.dilate(np.isin(lab, fast), a.dilate) if fast else np.zeros(lab.shape, bool)
        x = np.clip(win["pix"][:, 0].astype(int), 0, lab.shape[1] - 1)
        y = np.clip(win["pix"][:, 1].astype(int), 0, lab.shape[0] - 1)
        inB, inI = ~excl[y, x], fmask[y, x]
        useI = int(inI.sum()) >= a.min_fast_events
        n_ev.append(win["n"]); n_fast += int(useI)

        x0, y0 = rc.warp(win, cal_ts, cal_w, None)
        xg, yg = rc.gt_warp(win, ts_gt, q_gt)
        e0 = np.hypot(x0 - xg, y0 - yg)
        raw["p50"].append(np.median(e0)); raw["p95"].append(np.percentile(e0, 95))
        v0B = rc.iwe(x0[inB], y0[inB]).var()
        v0I = rc.iwe(x0[inI], y0[inI]).var() if useI else None
        gt_ref["gB"].append(rc.iwe(xg[inB], yg[inB]).var() / v0B)
        if useI:
            gt_ref["gI"].append(rc.iwe(xg[inI], yg[inI]).var() / v0I)
        for name, d in modes:
            xs, ys = rc.warp(win, cal_ts, cal_w, d)
            e = np.hypot(xs - xg, ys - yg)
            r = res[name]
            r["p50"].append(np.median(e)); r["p95"].append(np.percentile(e, 95))
            r["gB"].append(rc.iwe(xs[inB], ys[inB]).var() / v0B)
            if useI:
                r["gI"].append(rc.iwe(xs[inI], ys[inI]).var() / v0I)

    n_w = len(n_ev)
    if n_w == 0:
        sys.exit("żadne okno nie miało dość zdarzeń")
    med = lambda x: float(np.median(x)) if len(x) else float("nan")
    print(f"Sekwencja: {seq.name} | okien: {n_w} x {a.window_ms:.0f} ms (zdarzeń w oknie: mediana {int(np.median(n_ev))}) | "
          f"z szybkimi IMO (v > {a.fast_v} m/s): {n_fast} | kalibracja: {a.calib}")
    print("Wszystkie wartości: mediana po oknach. Błąd = odległość od kompensacji orientacją z GT [px].\n")
    print(f"{'wariant':28s} {'błąd p50':>9s} {'błąd p95':>9s} {'zysk tło':>9s} {'zysk IMO':>9s}")
    print(f"{'bez kompensacji (vs GT)':28s} {med(raw['p50']):9.2f} {med(raw['p95']):9.2f} {1.0:9.2f} {1.0:9.2f}")
    print(f"{'orientacja z GT (wzorzec)':28s} {0.0:9.2f} {0.0:9.2f} {med(gt_ref['gB']):9.2f} {med(gt_ref['gI']) if gt_ref['gI'] else float('nan'):9.2f}")
    out = dict(sequence=seq.name, windows=n_w, window_ms=a.window_ms, calib=a.calib, windows_with_fast_imo=n_fast,
               raw=dict(p50=med(raw["p50"]), p95=med(raw["p95"])), gt_gain_bg=med(gt_ref["gB"]), modes={})
    for name, d in modes:
        r = res[name]
        print(f"{name:28s} {med(r['p50']):9.2f} {med(r['p95']):9.2f} {med(r['gB']):9.2f} {med(r['gI']):9.2f}")
        out["modes"][name] = dict(lag_ms=d * 1e3, err_p50=med(r["p50"]), err_p95=med(r["p95"]), gain_bg=med(r["gB"]), gain_imo=med(r["gI"]))
    print("\nBłąd mierzy tylko niedokładność rotacji z IMU (kalibracja, bias, opóźnienie, szum); nie obejmuje translacji ani ruchu obiektów."
          "\n'Zysk tło' < 'orientacja z GT' pokazuje, ile ostrości traci się przez błędy IMU; sam zysk z GT jest ograniczony przez translację kamery.")
    if d_oracle is None:
        print("\n[UWAGA] brak gyro_vs_gt.json: pominięto warianty oracle (potrzebne opóźnienie z GT).")
    odir = Path(a.out) if a.out else ec.result_dir(seq.name, a.out_root)
    odir.mkdir(parents=True, exist_ok=True)
    (odir / "rotcomp_eval.json").write_text(json.dumps(out, indent=1))
    print(f"\nZapisano {odir / 'rotcomp_eval.json'}")


if __name__ == "__main__":
    main()
