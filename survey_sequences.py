#!/usr/bin/env python3
"""Przegląd wszystkich sekwencji EVIMO2 w katalogu: tabela do wyboru sekwencji do pracy.

Użycie:
    python3 survey_sequences.py <katalog_z_danymi> [--csv survey.csv] [--sort rot_px_10ms]
                                [--step 5] [--move-thr-m 0.05] [--move-thr-deg 5]

Szuka rekurencyjnie katalogów z dataset_info.npz. Dla każdej sekwencji liczy:
  * zdarzenia: czas trwania, przepustowość w oknach 10 ms (średnia/p95/max) - te same
    okna co `evimo2_viewer --stats-only`, więc wartości są porównywalne,
  * IMU lewej kamery: |omega| (średnia/p95/max) i maks. przesunięcie od rotacji w 10 ms [px],
  * ground truth z full_trajectory (200 Hz): ruch kamery i obiektów.

UWAGA - konwencja poz nie jest jeszcze potwierdzona (T_świat<-kamera czy odwrotnie):
  * kąty obrotu i prędkość kątową kamery liczymy z kąta względnego obrotu - nie zależy
    od konwencji,
  * prędkość liniową kamery podajemy w dwóch hipotezach: `cam_v_max` (t = pozycja kamery)
    i `cam_v_max_alt` (pozycja = -R^T t); jeśli się mocno różnią, konwencja ma znaczenie,
  * `n_obj_moving` zakłada, że pozy obiektów są w stałym układzie świata (obiekt statyczny =
    stała poza); to hipoteza, nie ustalony fakt.
Prędkości liczone są na odstępie --step próbek (domyślnie 5 = 25 ms) dla redukcji szumu.
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

import evimo2_common as ec

CHUNK = 4_000_000

PRINT_COLS = [  # (klucz, nagłówek, format)
    ("name", "sekwencja", "{}"), ("dur_s", "czas[s]", "{:.2f}"),
    ("rate_mean_mev_s", "śr.Mev/s", "{:.2f}"), ("rate_max_mev_s", "max Mev/s", "{:.2f}"),
    ("imu_w_max", "ω_imu max", "{:.3f}"), ("rot_px_10ms", "rot px/10ms", "{:.2f}"),
    ("cam_w_max_gt", "ω_gt max", "{:.3f}"), ("cam_v_max", "v_cam max", "{:.3f}"),
    ("cam_rot_deg", "rot_cam[°]", "{:.1f}"), ("n_obj", "obiekty", "{}"),
    ("n_obj_moving", "ruchome*", "{}"),
]


# ---------- zdarzenia ----------
def event_stats(seq):
    t, _, p = ec.load_events(seq)
    n = len(t)
    t0, t1 = float(t[0]), float(t[-1])
    dur = t1 - t0
    nb = int(dur * 100) + 2
    hist = np.zeros(nb, dtype=np.int64)
    on = 0
    for i in range(0, n, CHUNK):
        c = np.asarray(t[i:i + CHUNK], dtype=np.float64)
        hist += np.bincount(((c - t0) * 100).astype(np.int64), minlength=nb)[:nb]
        on += int((np.asarray(p[i:i + CHUNK]) != 0).sum())
    full = hist[:max(1, int(dur * 100))]           # tylko pełne okna 10 ms
    rate = full * 100 / 1e6                         # zdarzenia/10ms -> Mev/s
    return dict(dur_s=dur, n_events=n, rate_mean_mev_s=n / max(dur, 1e-9) / 1e6,
                rate_p95_mev_s=float(np.percentile(rate, 95)), rate_max_mev_s=float(rate.max()),
                on_ratio=on / n)


# ---------- IMU ----------
def imu_stats(info, fx):
    cand = [(p, l) for p, l in ec.find_series(info).items() if ec.is_imu_series(l)]
    left = [c for c in cand if "left" in c[0]] or cand
    if not left:
        raise ValueError("brak IMU w meta")
    _, gyro, _ = ec.imu_arrays(left[0][1])
    w = np.linalg.norm(gyro, axis=1)
    return dict(imu_w_mean=float(w.mean()), imu_w_p95=float(np.percentile(w, 95)),
                imu_w_max=float(w.max()), rot_px_10ms=float(w.max() * fx * 0.01))


# ---------- pozy ----------
def qconj(q):
    return q * np.array([1.0, -1.0, -1.0, -1.0])


def qmul(a, b):
    aw, ax, ay, az = a.T
    bw, bx, by, bz = b.T
    return np.stack([aw * bw - ax * bx - ay * by - az * bz,
                     aw * bx + ax * bw + ay * bz - az * by,
                     aw * by - ax * bz + ay * bw + az * bx,
                     aw * bz + ax * by - ay * bx + az * bw], axis=1)


def qangle(a, b):
    """Kąt względnego obrotu [rad]; ta sama wartość dla obu konwencji poz."""
    r = qmul(qconj(a), b)
    return 2 * np.arctan2(np.linalg.norm(r[:, 1:], axis=1), np.abs(r[:, 0]))


def rot_inv(q, t):
    """R^T t policzone kwaternionami."""
    p = np.concatenate([np.zeros((len(t), 1)), t], axis=1)
    return qmul(qmul(qconj(q), p), q)[:, 1:]


def pose_series(traj, key):
    rows = []
    for e in traj:
        ent = e.get(key)
        pos = ent.get("pos") if isinstance(ent, dict) else None
        if isinstance(pos, dict) and "q" in pos and "t" in pos:
            rows.append((float(ent["ts"]), [pos["q"][k] for k in "wxyz"], [pos["t"][k] for k in "xyz"]))
    if len(rows) < 2:
        return None
    rows.sort(key=lambda r: r[0])
    ts = np.array([r[0] for r in rows])
    q = np.array([r[1] for r in rows], dtype=float)
    t = np.array([r[2] for r in rows], dtype=float)
    return ts, q / np.linalg.norm(q, axis=1, keepdims=True), t


def motion(ts, q, t, step):
    step = max(1, min(step, len(ts) - 1))
    dt = ts[step:] - ts[:-step]
    ok = dt > 0
    w = qangle(q[:-step], q[step:])[ok] / dt[ok]
    v = np.linalg.norm(t[step:] - t[:-step], axis=1)[ok] / dt[ok]
    c = -rot_inv(q, t)                               # hipoteza: poza = T_kamera<-świat
    v_alt = np.linalg.norm(c[step:] - c[:-step], axis=1)[ok] / dt[ok]
    return dict(w_max=float(w.max()), v_max=float(v.max()), v_alt_max=float(v_alt.max()),
                rot_deg=float(np.degrees(qangle(np.repeat(q[:1], len(q), 0), q)).max()),
                disp_m=float(np.linalg.norm(t - t[0], axis=1).max()))


def gt_stats(info, seq, step, thr_m, thr_deg):
    meta = info["meta"]
    traj = meta.get("full_trajectory") or []
    out = dict(gt_frames=len(meta.get("frames", [])))
    for fn, key in (("dataset_mask.npz", "mask_frames"), ("dataset_depth.npz", "depth_frames")):
        pth = Path(seq) / fn
        if pth.is_file():
            with np.load(pth) as z:
                out[key] = len(z.files)
    cam = pose_series(traj, "cam")
    if cam is None:
        raise ValueError("brak poz kamery w full_trajectory")
    m = motion(*cam, step)
    out.update(cam_w_max_gt=m["w_max"], cam_v_max=m["v_max"], cam_v_max_alt=m["v_alt_max"],
               cam_rot_deg=m["rot_deg"], cam_disp_m=m["disp_m"])
    obj_ids = sorted({k for e in traj for k in e if isinstance(k, str) and k.isdigit()}, key=int)
    moving, disp_max, rot_max, n_obj = 0, 0.0, 0.0, 0
    for k in obj_ids:
        s = pose_series(traj, k)
        if s is None:
            continue
        n_obj += 1
        mo = motion(*s, step)
        disp_max, rot_max = max(disp_max, mo["disp_m"]), max(rot_max, mo["rot_deg"])
        moving += int(mo["disp_m"] > thr_m or mo["rot_deg"] > thr_deg)
    out.update(n_obj=n_obj, n_obj_moving=moving, obj_disp_max_m=disp_max, obj_rot_max_deg=rot_max)
    return out


def survey_sequence(seq, step=5, thr_m=0.05, thr_deg=5.0):
    seq = Path(seq)
    info = ec.load_info(seq)
    if info is None:
        raise ValueError("brak dataset_info.npz")
    intr = ec.find_dict_with(info, ["fx", "fy", "cx", "cy"])
    if intr is None:
        raise ValueError("brak intrinsics w meta")
    row = dict(name=seq.name, split=seq.parent.name)
    row.update(event_stats(seq))
    row.update(imu_stats(info, float(intr["fx"])))
    row.update(gt_stats(info, seq, step, thr_m, thr_deg))
    return row


# ---------- CLI ----------
def print_table(rows):
    cells = [[h for _, h, _ in PRINT_COLS]]
    for r in rows:
        cells.append([(f.format(r[k]) if k in r else "-") for k, _, f in PRINT_COLS])
    widths = [max(len(c[i]) for c in cells) for i in range(len(PRINT_COLS))]
    for i, c in enumerate(cells):
        print("  ".join(x.ljust(w) if j == 0 else x.rjust(w) for j, (x, w) in enumerate(zip(c, widths))))
        if i == 0:
            print("  ".join("-" * w for w in widths))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--csv", default="survey.csv")
    ap.add_argument("--sort", default="rot_px_10ms", help="klucz kolumny, malejąco")
    ap.add_argument("--step", type=int, default=5)
    ap.add_argument("--move-thr-m", type=float, default=0.05)
    ap.add_argument("--move-thr-deg", type=float, default=5.0)
    a = ap.parse_args()

    seqs = sorted(p.parent for p in Path(a.root).rglob("dataset_info.npz"))
    if not seqs:
        sys.exit(f"Brak sekwencji (dataset_info.npz) w {a.root}")
    rows, errors = [], []
    for i, s in enumerate(seqs, 1):
        print(f"[{i}/{len(seqs)}] {s.name}", file=sys.stderr)
        try:
            rows.append(survey_sequence(s, a.step, a.move_thr_m, a.move_thr_deg))
        except Exception as ex:                       # jedna zła sekwencja nie przerywa przeglądu
            errors.append((s.name, str(ex)))
    rows.sort(key=lambda r: r.get(a.sort, float("-inf")), reverse=True)

    print()
    print_table(rows)
    print("\n* ruchome = obiekty o przemieszczeniu > %.2f m lub obrocie > %.0f° od pierwszej pozy,"
          " przy założeniu, że pozy są w stałym układzie świata (niepotwierdzone)."
          % (a.move_thr_m, a.move_thr_deg))
    for name, err in errors:
        print(f"BŁĄD {name}: {err}")

    if rows:
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with open(a.csv, "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=keys)
            wr.writeheader()
            wr.writerows(rows)
        print(f"\nZapisano {a.csv} ({len(rows)} sekwencji, {len(keys)} kolumn)")


if __name__ == "__main__":
    main()