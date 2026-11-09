#!/usr/bin/env python3
"""Etykiety ruchu obiektów w czasie i maski IMO (niezależnie poruszające się obiekty) z ground truth.

Użycie:
    python3 imo_labels.py <katalog_sekwencji> [--cls R*S] [--v-thr 0.05] [--w-thr 1.0] [--half-ms 25]
                          [--validate] [--act-ms 10] [--dilate 6] [--save-masks] [--out DIR]
    python3 imo_labels.py --all <katalog_z_sekwencjami> [te same opcje] [--csv PLIK]

Dla każdej klatki z maską (dataset_mask.npz, ~60 Hz) liczy prędkość liniową v i kątową ω każdego
obiektu OBECNEGO W MASCE, w układzie świata (pozy obiektów składane z pozą kamery - klasa konwencji
z object_motion.py, domyślnie R*S), różnicą centralną o półszerokości --half-ms. Obiekt jest
'ruchomy' w tej klatce, gdy v > --v-thr [m/s] lub ω > --w-thr [rad/s]. Maska IMO to suma pikseli
ruchomych obiektów. Progi są WSTĘPNE (ok. 3x poziom szumu statycznych obiektów: v p95 <= 1,7 cm/s,
ω p95 do 0,3 rad/s); histogramy prędkości na wydruku służą do ich korekty.
Żadna klatka ani sekwencja nie jest odrzucana; jakość etykiet ocenia opcjonalna walidacja.

--validate: porównuje gęstość zdarzeń (okno +-act-ms wokół klatki) w masce IMO rozszerzonej o --dilate
pikseli i poza nią. Przy nieruchomej kamerze ruchome obiekty powinny generować wyraźnie większą
gęstość zdarzeń niż reszta (stosunek >> 1). Zakłada, że czas GT i zdarzeń jest zgodny (niesprawdzone;
niedokładność rzędu opóźnienia IMU względem GT osłabia wynik).
Zapisuje <out>/imo_labels.json (+ imo_masks.npz z --save-masks).
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

import evimo2_common as ec
import object_motion as om
from survey_sequences import qangle

V_THRESHOLDS = [0.05, 0.1, 0.3, 1.0]      # stałe progi do zestawienia udziału IMO w funkcji prędkości
V_EDGES = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.3, 1.0]
W_EDGES = [0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]


def kinematics(tracks, cls, half_s):
    """v [m/s] i ω [rad/s] każdego obiektu w świecie w chwilach próbek pozy (różnica centralna)."""
    out = {}
    for k, t in tracks.items():
        q, p = om.world_pose(cls, t["qc"], t["tc"], t["qs"], t["ts_"])
        ts = t["ts"]
        n = len(ts)
        s = max(1, int(round(half_s / max(float(np.median(np.diff(ts))), 1e-6))))
        i = np.arange(n)
        a, b = np.clip(i - s, 0, n - 1), np.clip(i + s, 0, n - 1)
        dt = ts[b] - ts[a]
        ok = dt > 0
        v, w = np.zeros(n), np.zeros(n)
        v[ok] = np.linalg.norm(p[b] - p[a], axis=1)[ok] / dt[ok]
        w[ok] = qangle(q[a], q[b])[ok] / dt[ok]
        out[k] = dict(ts=ts, v=v, w=w)
    return out


def dilate(m, r):
    out = m.copy()
    for _ in range(r):
        p = out.copy()
        p[1:] |= out[:-1]
        p[:-1] |= out[1:]
        q = p.copy()
        q[:, 1:] |= p[:, :-1]
        q[:, :-1] |= p[:, 1:]
        out = q
    return out


def hist_line(vals, edges, unit):
    h, _ = np.histogram(vals, bins=edges + [np.inf])
    names = [f"[{a:g},{b:g})" for a, b in zip(edges, edges[1:])] + [f">={edges[-1]:g}"]
    return "  " + "  ".join(f"{nm}: {n}" for nm, n in zip(names, h)) + f"   ({unit}, n = {len(vals)})"


def process(seq, a, cls):
    info = ec.load_info(seq)
    _, tracks = om.load_tracks(info)
    kin = kinematics(tracks, cls, a.half_ms * 1e-3)
    fts = {int(f["id"]): float(f["ts"]) for f in info["meta"].get("frames", [])}
    mpath = Path(seq) / "dataset_mask.npz"
    if not mpath.is_file():
        raise ValueError("brak dataset_mask.npz")

    ev = None
    if a.validate:
        t_ev, xy_ev, _ = ec.load_events(seq)
        ev = (t_ev, xy_ev)
    frames, masks_out = [], {}
    area_by_v = [[] for _ in V_THRESHOLDS]
    per_obj, vs, ws = {}, [], []
    val = dict(n_in=0, n_out=0, a_in=0, a_out=0, frames=0)
    unknown_ids, no_ts = set(), 0
    with np.load(mpath) as z:
        for key in sorted(z.files, key=lambda s: int(s.split("_")[-1])):
            fid = int(key.split("_")[-1])
            tsf = fts.get(fid)
            if tsf is None:
                no_ts += 1
                continue
            m = z[key]
            lab = m // 1000
            ids, cnt = np.unique(lab[lab > 0], return_counts=True)
            moving, moving_v = [], []
            area_v = [0.0] * len(V_THRESHOLDS)
            for oid, c in zip(ids, cnt):
                kk = str(int(oid))
                if kk not in kin:
                    unknown_ids.add(kk)
                    continue
                v = float(np.interp(tsf, kin[kk]["ts"], kin[kk]["v"]))
                w = float(np.interp(tsf, kin[kk]["ts"], kin[kk]["w"]))
                for j, th in enumerate(V_THRESHOLDS):
                    if v > th:
                        area_v[j] += c / m.size
                mv_v, mv_w = v > a.v_thr, w > a.w_thr
                mv = mv_v or mv_w
                vs.append(v); ws.append(w)
                o = per_obj.setdefault(kk, dict(present=0, moving=0, w_only=0, area_sum=0.0, area_mov_sum=0.0, v=[]))
                o["present"] += 1; o["area_sum"] += c / m.size; o["v"].append(v)
                if mv:
                    o["moving"] += 1; o["area_mov_sum"] += c / m.size
                    o["w_only"] += int(not mv_v)
                    moving.append(int(oid))
                if mv_v:
                    moving_v.append(int(oid))
            imo = np.isin(lab, moving) if moving else np.zeros(m.shape, bool)
            for j in range(len(V_THRESHOLDS)):
                area_by_v[j].append(area_v[j])
            imo_v = np.isin(lab, moving_v) if moving_v else np.zeros(m.shape, bool)
            frames.append(dict(id=fid, ts=tsf, moving_ids=[str(i) for i in moving], imo_frac=float(imo.mean()),
                               imo_frac_v_only=float(imo_v.mean())))
            if a.save_masks:
                masks_out[f"imo_{fid:010d}"] = imo.astype(np.uint8)
            if ev is not None and imo.any():
                t_ev, xy_ev = ev
                i0, i1 = np.searchsorted(t_ev, [tsf - a.act_ms * 1e-3, tsf + a.act_ms * 1e-3])
                sel = np.asarray(xy_ev[i0:i1])
                x = np.clip(sel[:, 0].astype(int), 0, m.shape[1] - 1)
                y = np.clip(sel[:, 1].astype(int), 0, m.shape[0] - 1)
                dil = dilate(imo, a.dilate)
                inside = dil[y, x]
                val["n_in"] += int(inside.sum()); val["n_out"] += int((~inside).sum())
                val["a_in"] += int(dil.sum()); val["a_out"] += int(dil.size - dil.sum()); val["frames"] += 1
    if not frames:
        raise ValueError("żadna maska nie ma znacznika czasu w meta['frames']")

    fr = np.array([f["imo_frac"] for f in frames])
    fv = np.array([f["imo_frac_v_only"] for f in frames])
    res = dict(sequence=Path(seq).name, class_used=cls, thresholds=dict(v_thr=a.v_thr, w_thr=a.w_thr, half_ms=a.half_ms),
               n_frames=len(frames), frames_without_ts=no_ts, ids_without_pose=sorted(unknown_ids),
               imo_frac=dict(mean=float(fr.mean()), p95=float(np.percentile(fr, 95)), max=float(fr.max()),
                             frames_with_imo=float((fr > 0).mean())),
               imo_frac_v_only=dict(mean=float(fv.mean()), p95=float(np.percentile(fv, 95)), max=float(fv.max()),
                                    frames_with_imo=float((fv > 0).mean())),
               imo_frac_by_v_thr={str(th): dict(mean=float(np.mean(x)), p95=float(np.percentile(x, 95)),
                                                frames_with_imo=float((np.array(x) > 0).mean()))
                                  for th, x in zip(V_THRESHOLDS, area_by_v)},
               objects={k: dict(frames_present=o["present"], frames_moving=o["moving"],
                                frames_moving_only_by_w=o["w_only"],
                                frac_moving=o["moving"] / o["present"], mean_area_frac=o["area_sum"] / o["present"],
                                mean_area_frac_when_moving=(o["area_mov_sum"] / o["moving"]) if o["moving"] else 0.0,
                                v_p95=float(np.percentile(o["v"], 95))) for k, o in sorted(per_obj.items(), key=lambda x: int(x[0]))},
               validation=None, per_frame=frames)
    if ev is not None and val["frames"] and val["n_out"] and val["a_in"]:
        dens_in = val["n_in"] / val["a_in"]
        dens_out = val["n_out"] / val["a_out"]
        res["validation"] = dict(frames=val["frames"], density_ratio=float(dens_in / dens_out),
                                 frac_events_in_imo=val["n_in"] / (val["n_in"] + val["n_out"]),
                                 frac_area_in_imo=val["a_in"] / (val["a_in"] + val["a_out"]), act_ms=a.act_ms, dilate_px=a.dilate)
    return res, np.array(vs), np.array(ws), masks_out


def print_report(r, vs, ws):
    t = r["thresholds"]
    print(f"Sekwencja: {r['sequence']} | klasa {r['class_used']} | klatek z maską: {r['n_frames']} | "
          f"ruchomy: v > {t['v_thr']} m/s lub ω > {t['w_thr']} rad/s (różnica centralna ±{t['half_ms']:.0f} ms)")
    if r["frames_without_ts"]:
        print(f"[UWAGA] {r['frames_without_ts']} masek bez znacznika czasu w meta['frames'] - pominięte.")
    if r["ids_without_pose"]:
        print(f"[UWAGA] obiekty w maskach bez póz (traktowane jako nie-ruchome): {', '.join(r['ids_without_pose'])}")
    print(f"\n{'id':>4} {'klatek':>7} {'ruchomy':>8} {'w tym tylko ω':>14} {'śr. obszar':>11} {'obszar gdy ruch.':>17} {'v p95':>7}")
    print(f"{'':>4} {'':>7} {'[%]':>8} {'[% klatek]':>14} {'[% obrazu]':>11} {'[% obrazu]':>17} {'[m/s]':>7}")
    for k, o in r["objects"].items():
        print(f"{k:>4} {o['frames_present']:7d} {100 * o['frac_moving']:8.0f} {100 * o['frames_moving_only_by_w'] / o['frames_present']:14.0f} "
              f"{100 * o['mean_area_frac']:11.1f} {100 * o['mean_area_frac_when_moving']:17.1f} {o['v_p95']:7.3f}")
    im = r["imo_frac"]
    print(f"\nIMO: średnio {100 * im['mean']:.2f}% obrazu, p95 {100 * im['p95']:.2f}%, max {100 * im['max']:.2f}%; "
          f"klatek z ruchomym obiektem: {100 * im['frames_with_imo']:.0f}%")
    iv = r["imo_frac_v_only"]
    print(f"IMO tylko z v (bez obiektów oznaczonych wyłącznie przez ω): średnio {100 * iv['mean']:.2f}%, p95 {100 * iv['p95']:.2f}%, "
          f"klatek z ruchomym obiektem: {100 * iv['frames_with_imo']:.0f}%")
    print("Udział IMO (tylko v) w funkcji progu prędkości: " + "; ".join(
        f"v > {th}: śr. {100 * d['mean']:.2f}%, p95 {100 * d['p95']:.2f}%" for th, d in r["imo_frac_by_v_thr"].items()))
    v = r["validation"]
    if v:
        print(f"Walidacja ({v['frames']} klatek z IMO; zdarzenia ±{v['act_ms']} ms, maska + {v['dilate_px']} px): gęstość zdarzeń "
              f"w IMO / poza = {v['density_ratio']:.1f}x; w IMO {100 * v['frac_events_in_imo']:.0f}% zdarzeń na {100 * v['frac_area_in_imo']:.1f}% powierzchni")
    print("\nRozkład prędkości obiektów w chwilach klatek (wszystkie obiekty obecne w maskach):")
    print(hist_line(vs, V_EDGES, "v [m/s]"))
    print(hist_line(ws, W_EDGES, "ω [rad/s]"))


def save(res, masks, out):
    out.mkdir(parents=True, exist_ok=True)
    (out / "imo_labels.json").write_text(json.dumps(res, indent=1))
    if masks:
        np.savez_compressed(out / "imo_masks.npz", **masks)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("seq_dir")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--cls", default="R*S", choices=om.CLASSES)
    ap.add_argument("--v-thr", type=float, default=0.05)
    ap.add_argument("--w-thr", type=float, default=1.0)
    ap.add_argument("--half-ms", type=float, default=25.0)
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--act-ms", type=float, default=10.0)
    ap.add_argument("--dilate", type=int, default=6)
    ap.add_argument("--save-masks", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--out-root", default=None)
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()

    if not a.all:
        seq = Path(a.seq_dir.rstrip("/"))
        res, vs, ws, masks = process(seq, a, a.cls)
        print_report(res, vs, ws)
        out = Path(a.out) if a.out else ec.result_dir(seq.name, a.out_root)
        save(res, masks, out)
        print(f"\nZapisano {out / 'imo_labels.json'}" + (f" i imo_masks.npz ({len(masks)} masek)" if masks else ""))
        return

    seqs = sorted(p.parent for p in Path(a.seq_dir).rglob("dataset_info.npz"))
    if not seqs:
        sys.exit(f"Brak sekwencji w {a.seq_dir}")
    rows, allv, allw = [], [], []
    for s in seqs:
        try:
            res, vs, ws, masks = process(s, a, a.cls)
        except Exception as ex:
            print(f"BŁĄD {s.name}: {ex}", file=sys.stderr)
            continue
        save(res, masks, ec.result_dir(s.name, a.out_root))
        allv.append(vs); allw.append(ws)
        v = res["validation"]
        rows.append(dict(sequence=s.name, frames=res["n_frames"], imo_mean_pct=100 * res["imo_frac"]["mean"],
                         imo_p95_pct=100 * res["imo_frac"]["p95"], frames_with_imo_pct=100 * res["imo_frac"]["frames_with_imo"],
                         imo_v_only_mean_pct=100 * res["imo_frac_v_only"]["mean"], imo_v_only_p95_pct=100 * res["imo_frac_v_only"]["p95"],
                         imo_v01_mean_pct=100 * res["imo_frac_by_v_thr"]["0.1"]["mean"], imo_v01_p95_pct=100 * res["imo_frac_by_v_thr"]["0.1"]["p95"],
                         imo_v03_mean_pct=100 * res["imo_frac_by_v_thr"]["0.3"]["mean"], imo_v03_p95_pct=100 * res["imo_frac_by_v_thr"]["0.3"]["p95"],
                         density_ratio=v["density_ratio"] if v else "", events_in_imo_pct=100 * v["frac_events_in_imo"] if v else ""))
    print(f"{'sekwencja':28s} {'klatek':>7s} {'IMO śr.%':>9s} {'IMO p95%':>9s} {'v>0.05 śr.%':>12s} {'v>0.1 śr.%':>11s} {'v>0.1 p95%':>11s} {'v>0.3 śr.%':>11s} {'v>0.3 p95%':>11s} {'klatki z IMO %':>15s} {'gęstość IMO/tło':>16s}")
    for r in rows:
        dr = f"{r['density_ratio']:.1f}x" if r["density_ratio"] != "" else "-"
        print(f"{r['sequence']:28s} {r['frames']:7d} {r['imo_mean_pct']:9.2f} {r['imo_p95_pct']:9.2f} {r['imo_v_only_mean_pct']:12.2f} "
              f"{r['imo_v01_mean_pct']:11.2f} {r['imo_v01_p95_pct']:11.2f} {r['imo_v03_mean_pct']:11.2f} {r['imo_v03_p95_pct']:11.2f} "
              f"{r['frames_with_imo_pct']:15.0f} {dr:>16s}")
    print("\nRozkład prędkości obiektów (wszystkie sekwencje, wszystkie klatki):")
    print(hist_line(np.concatenate(allv), V_EDGES, "v [m/s]"))
    print(hist_line(np.concatenate(allw), W_EDGES, "ω [rad/s]"))
    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)
        print(f"\nZapisano {a.csv}")


if __name__ == "__main__":
    main()