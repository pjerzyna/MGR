#!/usr/bin/env python3
"""Ruch obiektów w układzie świata: konwencja poz obiektów i lista obiektów ruchomych.

Użycie:
    python3 object_motion.py <katalog_sekwencji> [--out DIR] [--step 5] [--static-deg 2]
                             [--move-m 0.05] [--move-deg 5] [--cls KLASA]
    python3 object_motion.py --all <katalog_z_sekwencjami> [--csv PLIK] [--out-root DIR]

Problem: poza kamery w full_trajectory to T_świat<-kamera (ustalone w gyro_vs_gt.py), ale nie
wiadomo, w jakim układzie zapisano pozy obiektów (świat? kamera? odwrotności?). Dla rotacji kamery
R (= R_wc) i zapisanej rotacji obiektu S obiekt NIERUCHOMY w świecie daje stały dokładnie jeden z
niezmienników:
    S      - S zapisane w układzie świata (T_w<-o)           [ROZRÓŻNIALNE od 'Rt*S' i 'S*R' poniżej]
    R*S    - S zapisane w układzie kamery (T_c<-o)           -> R_wo = R S
    S*Rt   - S zapisane jako T_o<-c (kamera w układzie obiektu) -> R_wo = R S^T
    Rt*S   - obiekt przyklejony do kamery, zapis w świecie (T_w<-o)
    S*R    - obiekt przyklejony do kamery, zapis T_o<-w
Skrypt liczy odchylenie (p95 [deg]) każdego niezmiennika dla każdego obiektu i wybiera klasę,
w której najwięcej obiektów jest stałych. Kamera musi się obracać (>= 10 deg), inaczej klasy
są nierozróżnialne. Następnie składa pozy obiektów do układu świata i liczy ich ruch:
obiekt jest 'ruchomy', gdy p95 odchylenia od pierwszej pozy przekracza --move-m lub --move-deg.
Założenie dla składania translacji: obiekty zapisane tą samą rodziną konwencji co kamera.
Zapisuje <out>/object_motion.json (lista ruchomych obiektów do oceny segmentacji).
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

import evimo2_common as ec
from survey_sequences import qangle, qconj, qmul

CLASSES = ["S", "R*S", "S*R", "Rt*S", "S*Rt"]
MIN_CAM_ROT_DEG = 10.0


def qrot(q, v):
    p = np.concatenate([np.zeros((len(v), 1)), v], axis=1)
    return qmul(qmul(q, p), qconj(q))[:, 1:]


def entry_pose(ent):
    pos = ent.get("pos") if isinstance(ent, dict) else None
    if isinstance(pos, dict) and "q" in pos and "t" in pos:
        return float(ent["ts"]), [pos["q"][k] for k in "wxyz"], [pos["t"][k] for k in "xyz"]
    return None


def load_tracks(info):
    """Kamera i obiekty z full_trajectory, próbki tylko z wpisów, w których jest i kamera, i obiekt."""
    traj = info["meta"].get("full_trajectory") or []
    ids = sorted({k for e in traj for k in e if isinstance(k, str) and k.isdigit()}, key=int)
    rows = {k: [] for k in ids}
    cam = []
    for e in traj:
        c = entry_pose(e.get("cam"))
        if c is None:
            continue
        cam.append(c)
        for k in ids:
            o = entry_pose(e.get(k))
            if o is not None:
                rows[k].append((c[0], c[1], c[2], o[1], o[2]))
    if len(cam) < 3:
        raise ValueError("brak poz kamery w full_trajectory")

    def arr(r, i, n):
        return np.array([x[i] for x in r], dtype=float).reshape(len(r), n)

    def norm(q):
        return q / np.linalg.norm(q, axis=1, keepdims=True)

    cam.sort(key=lambda r: r[0])
    camd = dict(ts=np.array([c[0] for c in cam]), q=norm(np.array([c[1] for c in cam], float)),
                t=np.array([c[2] for c in cam], float))
    tracks = {}
    for k, r in rows.items():
        if len(r) < 3:
            continue
        r.sort(key=lambda x: x[0])
        tracks[k] = dict(ts=arr(r, 0, 1)[:, 0], qc=norm(arr(r, 1, 4)), tc=arr(r, 2, 3), qs=norm(arr(r, 3, 4)), ts_=arr(r, 4, 3))
    return camd, tracks


def invariants(qc, qs):
    return {"S": qs, "R*S": qmul(qc, qs), "S*R": qmul(qs, qc), "Rt*S": qmul(qconj(qc), qs), "S*Rt": qmul(qs, qconj(qc))}


def dev_deg(q):
    return np.degrees(qangle(np.repeat(q[:1], len(q), 0), q))


def world_pose(cls, qc, tc, qs, ts_):
    if cls == "R*S":
        return qmul(qc, qs), qrot(qc, ts_) + tc
    if cls == "S*Rt":
        q = qmul(qc, qconj(qs))
        return q, tc - qrot(q, ts_)
    return qs, ts_                       # 'S' i klasy 'przyklejone do kamery': zapis traktowany jako T_w<-o


def motion_stats(ts, q, p, step):
    st = max(1, min(step, len(ts) - 1))
    dt = ts[st:] - ts[:-st]
    ok = dt > 0
    v = np.linalg.norm(p[st:] - p[:-st], axis=1)[ok] / dt[ok]
    w = qangle(q[:-st], q[st:])[ok] / dt[ok]
    drot, dpos = dev_deg(q), np.linalg.norm(p - p[0], axis=1)
    return dict(rot_p95_deg=float(np.percentile(drot, 95)), pos_p95_m=float(np.percentile(dpos, 95)),
                rot_max_deg=float(drot.max()), pos_max_m=float(dpos.max()),
                v_p95=float(np.percentile(v, 95)), w_p95=float(np.percentile(w, 95)))


def analyze(info, step=5, static_deg=2.0, move_m=0.05, move_deg=5.0, force_cls=None):
    cam, tracks = load_tracks(info)
    cam_rot = float(dev_deg(cam["q"]).max())
    cam_disp = float(np.linalg.norm(cam["t"] - cam["t"][0], axis=1).max())
    dev = {k: {c: float(np.percentile(dev_deg(v), 95)) for c, v in invariants(t["qc"], t["qs"]).items()}
           for k, t in tracks.items()}
    counts = {c: sum(d[c] < static_deg for d in dev.values()) for c in CLASSES}
    order = sorted(CLASSES, key=lambda c: -counts[c])
    best = order[0]
    distinguishable = cam_rot >= MIN_CAM_ROT_DEG and counts[order[0]] - counts[order[1]] >= 1
    used = force_cls or (best if distinguishable else None)
    objs = {}
    if used:
        for k, t in tracks.items():
            q, p = world_pose(used, t["qc"], t["tc"], t["qs"], t["ts_"])
            m = motion_stats(t["ts"], q, p, step)
            m["moving"] = bool(m["rot_p95_deg"] > move_deg or m["pos_p95_m"] > move_m)
            objs[k] = m
    return dict(cam_rot_deg=cam_rot, cam_disp_m=cam_disp, n_obj=len(tracks), dev=dev, counts=counts, best=best,
                distinguishable=bool(distinguishable), class_used=used, class_forced=bool(force_cls), objects=objs)


def print_report(name, r, static_deg):
    print(f"Sekwencja: {name} | kamera: obrót max {r['cam_rot_deg']:.1f} deg, przemieszczenie max {r['cam_disp_m']:.2f} m | obiektów: {r['n_obj']}")
    print(f"\nOdchylenie niezmienników rotacji, p95 [deg]; obiekt nieruchomy w świecie ma ~0 dla właściwej klasy (próg {static_deg} deg):")
    print(f"{'id':>4} " + " ".join(f"{c:>8}" for c in CLASSES))
    for k, d in r["dev"].items():
        print(f"{k:>4} " + " ".join(f"{d[c]:8.2f}" + ("*" if d[c] < static_deg else " ") for c in CLASSES))
    print(f"{'stałe':>4} " + " ".join(f"{r['counts'][c]:8d} " for c in CLASSES))
    if r["cam_rot_deg"] < MIN_CAM_ROT_DEG:
        print(f"\n[UWAGA] kamera obraca się tylko o {r['cam_rot_deg']:.1f} deg (< {MIN_CAM_ROT_DEG:.0f}): klasy są nierozróżnialne.")
    elif not r["distinguishable"]:
        print("\n[UWAGA] dwie klasy mają tyle samo stałych obiektów: konwencja niejednoznaczna dla tej sekwencji.")
    else:
        print(f"\nNajlepsza klasa: {r['best']} ({r['counts'][r['best']]} z {r['n_obj']} obiektów stałych)")
    if not r["class_used"]:
        print("Brak klasy do złożenia poz (podaj --cls albo użyj --all, żeby klasa wyszła z głosowania).")
        return
    print(f"\nRuch w układzie świata (klasa {r['class_used']}{', wymuszona' if r['class_forced'] else ''}); "
          "translacja składana przy założeniu tej samej rodziny konwencji co kamera:")
    print(f"{'id':>4} {'rot p95':>8} {'poz p95':>8} {'v p95':>8} {'ω p95':>8}  ruchomy")
    print(f"{'':>4} {'[deg]':>8} {'[cm]':>8} {'[m/s]':>8} {'[rad/s]':>8}")
    for k, m in r["objects"].items():
        print(f"{k:>4} {m['rot_p95_deg']:8.1f} {100 * m['pos_p95_m']:8.1f} {m['v_p95']:8.3f} {m['w_p95']:8.3f}  {'TAK' if m['moving'] else 'nie'}")
    mv = [k for k, m in r["objects"].items() if m["moving"]]
    print(f"\nRuchome: {len(mv)} z {len(r['objects'])}: {', '.join(mv) if mv else '-'}")


def save_json(path, name, r, args):
    out = dict(sequence=name, class_used=r["class_used"], class_best=r["best"], distinguishable=r["distinguishable"],
               cam_rot_deg=r["cam_rot_deg"], cam_disp_m=r["cam_disp_m"], class_counts=r["counts"],
               thresholds=dict(static_deg=args.static_deg, move_m=args.move_m, move_deg=args.move_deg),
               moving_ids=[k for k, m in r["objects"].items() if m["moving"]], objects=r["objects"])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2))


def run_all(a):
    seqs = sorted(Path(a.seq_dir).rglob("dataset_info.npz"))
    if not seqs:
        sys.exit(f"Brak sekwencji w {a.seq_dir}")
    infos, res = {}, {}
    for p in seqs:
        try:
            infos[p.parent.name] = ec.load_info(p.parent)
            res[p.parent.name] = analyze(infos[p.parent.name], a.step, a.static_deg, a.move_m, a.move_deg, a.cls)
        except Exception as ex:
            print(f"BŁĄD {p.parent.name}: {ex}", file=sys.stderr)
    votes = {}
    for r in res.values():
        if r["distinguishable"]:
            votes[r["best"]] = votes.get(r["best"], 0) + 1
    n_dist = sum(votes.values())
    print(f"Głosowanie klasy konwencji (sekwencje z rozróżnialnymi klasami: {n_dist} z {len(res)}):")
    for c, n in sorted(votes.items(), key=lambda x: -x[1]):
        print(f"  {c}: {n}")
    if not votes:
        sys.exit("Żadna sekwencja nie ma rozróżnialnych klas - nie da się ustalić konwencji (podaj --cls).")
    unanimous = len(votes) == 1
    ds_cls = next(iter(votes)) if unanimous else None
    if not unanimous:
        print("  [UWAGA] głosy niejednomyślne: konwencja nie jest jednolita w zbiorze, więc sekwencji o nierozróżnialnych "
              "klasach NIE klasyfikuję (podaj --cls, jeśli znasz konwencję).")
    rows = []
    for name, r in res.items():
        if r["class_used"] is None and ds_cls:              # nierozróżnialne -> klasa zbiorcza (tylko przy jednomyślności)
            r = analyze(infos[name], a.step, a.static_deg, a.move_m, a.move_deg, ds_cls)
            res[name] = r
        if r["class_used"]:
            save_json(ec.result_dir(name, a.out_root) / "object_motion.json", name, r, a)
        mv = [k for k, m in r["objects"].items() if m["moving"]]
        rows.append(dict(sequence=name, cam_rot_deg=r["cam_rot_deg"], n_obj=r["n_obj"], class_used=r["class_used"] or "?",
                         class_source=("wymuszona/głosowanie" if r["class_forced"] else "własna") if r["class_used"] else "nieustalona",
                         n_moving=len(mv) if r["class_used"] else "-", moving_ids=" ".join(mv) if r["class_used"] else "-"))
    print(f"\n{'sekwencja':28s} {'rot kam':>8s} {'obiekty':>8s} {'klasa':>6s} {'źródło klasy':>20s} {'ruchome':>8s}  id ruchomych")
    for r in sorted(rows, key=lambda x: x["sequence"]):
        print(f"{r['sequence']:28s} {r['cam_rot_deg']:8.1f} {r['n_obj']:8d} {r['class_used']:>6s} {r['class_source']:>20s} {str(r['n_moving']):>8s}  {r['moving_ids']}")
    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)
        print(f"\nZapisano {a.csv}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("seq_dir")
    ap.add_argument("--all", action="store_true", help="przejrzyj wszystkie sekwencje w katalogu")
    ap.add_argument("--out", default=None)
    ap.add_argument("--out-root", default=None)
    ap.add_argument("--csv", default=None)
    ap.add_argument("--step", type=int, default=5)
    ap.add_argument("--static-deg", type=float, default=2.0)
    ap.add_argument("--move-m", type=float, default=0.05)
    ap.add_argument("--move-deg", type=float, default=5.0)
    ap.add_argument("--cls", choices=CLASSES, default=None, help="wymuś klasę konwencji")
    a = ap.parse_args()
    if a.all:
        return run_all(a)
    seq = Path(a.seq_dir.rstrip("/"))
    info = ec.load_info(seq)
    if info is None:
        sys.exit("brak dataset_info.npz")
    r = analyze(info, a.step, a.static_deg, a.move_m, a.move_deg, a.cls)
    print_report(seq.name, r, a.static_deg)
    if r["class_used"]:
        out = Path(a.out) if a.out else ec.result_dir(seq.name, a.out_root)
        save_json(out / "object_motion.json", seq.name, r, a)
        print(f"Zapisano {out / 'object_motion.json'}")


if __name__ == "__main__":
    main()
