#!/usr/bin/env python3
"""Opóźnienie IMU względem ZDARZEŃ metodą maksymalizacji kontrastu (prototyp w Pythonie).

Użycie:
    python3 event_lag_check.py <katalog_sekwencji> [--json PLIK] [--out DIR] [--out-root DIR] [--windows 10]
                               [--window-ms 30] [--lag-range -60 30] [--lag-step 2] [--no-plot]
    python3 event_lag_check.py --summary [--out-root DIR] [--glob WZORZEC]

Wymaga wyniku gyro_vs_gt.py (szukany w $EVIMO2_OUT/<sceneXX_NN>/gyro_vs_gt.json, potem w
$EVIMO2_OUT/<pełna_nazwa>/; albo --json): rotacji R IMU->kamera, biasu i opóźnienia d_GT.
Wyniki zapisywane są obok gyro_vs_gt.json (albo w --out).

Metoda: wybieramy okna (domyślnie 10 x 30 ms) w chwilach największego |omega|. Dla każdego
kandydata d z siatki: omega_kamery(t) = R (omega_imu(t + d) - b), zdarzenia okna są obracane do
chwili środkowej okna (czysty obrót, z uwzględnieniem zniekształceń radtan), a z nich budowany
jest obraz (głosowanie dwuliniowe). Miarą ostrości jest wariancja obrazu. Optimum d_ev
(z interpolacją paraboliczną) jest opóźnieniem IMU względem zdarzeń - ta sama definicja
co d_GT w gyro_vs_gt.py: omega_kamery(t) = R (omega_imu(t + d) - b).

Interpretacja (potrzebne >= 2-3 sekwencje o różnym d_GT, patrz --summary):
  d_ev ~ d_GT        -> zdarzenia i GT dzielą oś czasu, dryfuje IMU,
  d_ev ~ const       -> to GT jest przesunięte względem zdarzeń (i IMU).
Ograniczenia: model czystego obrotu (bez translacji), obiekty ruchome zaburzają kontrast.
"""
import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np

import evimo2_common as ec

W, H = 640, 480
I3 = np.eye(3)


# ---------- geometria ----------
def undistort(px, py, intr, iters=8):
    fx, fy, cx, cy = (float(intr[k]) for k in ("fx", "fy", "cx", "cy"))
    k1, k2, k3, p1, p2 = (float(intr.get(k, 0.0)) for k in ("k1", "k2", "k3", "p1", "p2"))
    xd, yd = (px - cx) / fx, (py - cy) / fy
    x, y = xd.copy(), yd.copy()
    for _ in range(iters):
        r2 = x * x + y * y
        rad = 1 + k1 * r2 + k2 * r2 ** 2 + k3 * r2 ** 3
        x, y = (xd - (2 * p1 * x * y + p2 * (r2 + 2 * x * x))) / rad, (yd - (p1 * (r2 + 2 * y * y) + 2 * p2 * x * y)) / rad
    return x, y


def rodrigues_batch(rv):
    th = np.linalg.norm(rv, axis=1)
    small = th < 1e-6
    ths = np.where(small, 1.0, th)
    a = np.where(small, 1 - th ** 2 / 6, np.sin(ths) / ths)
    b = np.where(small, 0.5 - th ** 2 / 24, (1 - np.cos(ths)) / ths ** 2)
    K = np.zeros((len(rv), 3, 3))
    K[:, 0, 1], K[:, 0, 2], K[:, 1, 0] = -rv[:, 2], rv[:, 1], rv[:, 2]
    K[:, 1, 2], K[:, 2, 0], K[:, 2, 1] = -rv[:, 0], -rv[:, 1], rv[:, 0]
    return I3 + a[:, None, None] * K + b[:, None, None] * (K @ K)


# ---------- okno zdarzeń ----------
def prep_window(tc, wlen, t, xy, intr):
    i0, i1 = np.searchsorted(t, [tc - wlen / 2, tc + wlen / 2])
    te = np.asarray(t[i0:i1], dtype=np.float64)
    pix = np.asarray(xy[i0:i1], dtype=np.float64)
    bx, by = undistort(pix[:, 0], pix[:, 1], intr)
    nst = max(2, int(round(wlen / 0.001)))
    return dict(tc=tc, te=te, bearing=np.column_stack([bx, by, np.ones(len(te))]),
                grid=np.linspace(tc - wlen / 2, tc + wlen / 2, nst + 1), kref=nst // 2, n=len(te), intr=intr)


def warp(win, cal_ts, cal_w, d):
    """Pozycje pikseli zdarzeń po obrocie do chwili środkowej okna; d=None -> bez kompensacji."""
    fx, fy, cx, cy = (float(win["intr"][k]) for k in ("fx", "fy", "cx", "cy"))
    b = win["bearing"]
    if d is None:
        return fx * b[:, 0] + cx, fy * b[:, 1] + cy
    g, n = win["grid"], len(win["grid"]) - 1
    om = np.column_stack([np.interp(0.5 * (g[:-1] + g[1:]) + d, cal_ts, cal_w[:, k]) for k in range(3)])
    steps = rodrigues_batch(om * np.diff(g)[:, None])
    Rg = np.empty((n + 1, 3, 3))
    Rg[0] = I3
    for k in range(n):
        Rg[k + 1] = Rg[k] @ steps[k]
    ki = np.clip(np.searchsorted(g, win["te"], side="right") - 1, 0, n - 1)
    E = rodrigues_batch(om[ki] * (win["te"] - g[ki])[:, None])
    Ri = np.einsum("nij,njk->nik", Rg[ki], E)                       # R_wc(t_i)
    M = np.einsum("ji,njk->nik", Rg[win["kref"]], Ri)               # R_wc(t_ref)^T R_wc(t_i)
    b2 = np.einsum("nij,nj->ni", M, b)
    return fx * b2[:, 0] / b2[:, 2] + cx, fy * b2[:, 1] / b2[:, 2] + cy


def iwe(xs, ys):
    x0, y0 = np.floor(xs), np.floor(ys)
    ax, ay = xs - x0, ys - y0
    x0, y0 = x0.astype(np.int64), y0.astype(np.int64)
    ok = (x0 >= 0) & (x0 < W - 1) & (y0 >= 0) & (y0 < H - 1)
    x0, y0, ax, ay = x0[ok], y0[ok], ax[ok], ay[ok]
    idx = y0 * W + x0
    n = W * H
    img = (np.bincount(idx, (1 - ax) * (1 - ay), n) + np.bincount(idx + 1, ax * (1 - ay), n)
           + np.bincount(idx + W, (1 - ax) * ay, n) + np.bincount(idx + W + 1, ax * ay, n))
    return img.reshape(H, W)


def parabola_peak(grid, y):
    i = int(np.argmax(y))
    if 0 < i < len(y) - 1:
        y0, y1, y2 = y[i - 1:i + 2]
        den = y0 - 2 * y1 + y2
        if den < 0:
            return float(grid[i] + 0.5 * (y0 - y2) / den * (grid[1] - grid[0]))
    return float(grid[i])


def pick_windows(cal_ts, cal_w, t0, t1, k, wlen, min_omega, margin=0.1, min_sep=0.25):
    wm = np.linalg.norm(cal_w, axis=1)
    nb = max(1, int(round(wlen / np.median(np.diff(cal_ts)))))
    sm = np.convolve(wm, np.ones(nb) / nb, mode="same")
    valid = (cal_ts > max(cal_ts[0], t0) + margin) & (cal_ts < min(cal_ts[-1], t1) - margin)
    chosen = []
    for i in np.argsort(-np.where(valid, sm, -1)):
        if not valid[i] or sm[i] < min_omega:
            break
        if all(abs(cal_ts[i] - c) >= min_sep for c in chosen):
            chosen.append(float(cal_ts[i]))
            if len(chosen) == k:
                break
    return sorted(chosen)


# ---------- raport / wykres ----------
def make_plot(path, name, wins, lags, curves, d_gt, d_ev, cal_ts, cal_w, best, per_win):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1], hspace=0.28, wspace=0.12)
    win = wins[best]
    cases = [("bez kompensacji", None), (f"d = d_GT = {d_gt * 1e3:+.1f} ms", d_gt), (f"d = d_ev = {d_ev * 1e3:+.1f} ms", d_ev)]
    for j, (lab, d) in enumerate(cases):
        ax = fig.add_subplot(gs[0, j])
        img = iwe(*warp(win, cal_ts, cal_w, d))
        ax.imshow(img, cmap="gray_r", vmin=0, vmax=max(np.percentile(img, 99.7), 1e-9))
        ax.set_title(f"{lab}\nwariancja = {img.var():.3f}", fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
    ax = fig.add_subplot(gs[1, 0:2])
    for c in curves:
        ax.plot(lags * 1e3, c / c.max(), color="0.7", lw=0.8)
    ax.plot(lags * 1e3, np.mean([c / c.max() for c in curves], axis=0), "k", lw=2, label="średnia po oknach")
    ax.axvline(d_gt * 1e3, color="tab:blue", ls="--", label=f"d_GT = {d_gt * 1e3:+.1f} ms")
    ax.axvline(d_ev * 1e3, color="tab:red", ls="-", label=f"d_ev = {d_ev * 1e3:+.1f} ms")
    ax.set_xlabel("opóźnienie d [ms]"); ax.set_ylabel("kontrast / max (po oknach)")
    ax.set_title(f"{name}: ostrość obrazu zdarzeń vs opóźnienie IMU"); ax.legend(); ax.grid(alpha=0.3)
    ax = fig.add_subplot(gs[1, 2])
    ax.plot([p["omega"] for p in per_win], [p["d_opt_ms"] for p in per_win], "o")
    ax.axhline(d_gt * 1e3, color="tab:blue", ls="--"); ax.axhline(d_ev * 1e3, color="tab:red")
    ax.set_xlabel("|ω| w oknie [rad/s]"); ax.set_ylabel("optimum okna [ms]")
    ax.set_title("Optimum per okno"); ax.grid(alpha=0.3)
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def summary(pattern):
    rows = []
    for f in sorted(glob.glob(pattern)):
        j = json.load(open(f))
        rows.append((j["sequence"], j["d_gt_ms"], j["d_ev_ms"], j["d_ev_std_ms"]))
    if not rows:
        sys.exit(f"Brak plików: {pattern}\nevent_lag_check.json powstaje po uruchomieniu skryptu na pojedynczej "
                 "sekwencji (bez --summary); sprawdź też EVIMO2_OUT / --out-root.")
    print(f"{'sekwencja':28s} {'d_GT [ms]':>10s} {'d_ev [ms]':>10s} {'różnica':>8s} {'std okien':>10s}")
    for s, g, e, sd in rows:
        print(f"{s:28s} {g:10.1f} {e:10.1f} {e - g:8.1f} {sd:10.1f}")
    if len(rows) >= 3:
        g = np.array([r[1] for r in rows]); e = np.array([r[2] for r in rows])
        if np.ptp(g) < 5:
            print("\nZakres d_GT < 5 ms: za mały, by odróżnić 'd_ev ~ d_GT' od 'd_ev ~ const'.")
            return
        slope, icpt = np.polyfit(g, e, 1)
        print(f"\nd_ev = {slope:.2f} * d_GT {icpt:+.1f} ms (zakres d_GT: {np.ptp(g):.0f} ms, n={len(rows)})")
        if abs(slope - 1) < 0.3:
            print("-> nachylenie bliskie 1: zdarzenia i GT dzielą oś czasu; opóźnienie dotyczy IMU.")
        elif abs(slope) < 0.3:
            print("-> nachylenie bliskie 0: d_ev jest stałe, więc GT jest przesunięte względem zdarzeń/IMU.")
        else:
            print("-> nachylenie pośrednie: wynik niejednoznaczny (mało sekwencji, szum albo oba efekty).")


# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("seq_dir", nargs="?")
    ap.add_argument("--json", default=None, help="wynik gyro_vs_gt.py (domyślnie out/<seq>/gyro_vs_gt.json)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--windows", type=int, default=10)
    ap.add_argument("--window-ms", type=float, default=30.0)
    ap.add_argument("--min-omega", type=float, default=0.5, help="minimalne |ω| okna [rad/s]")
    ap.add_argument("--lag-range", type=float, nargs=2, default=[-60.0, 30.0], metavar=("MIN", "MAX"))
    ap.add_argument("--lag-step", type=float, default=2.0, help="krok siatki [ms]")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--summary", action="store_true")
    ap.add_argument("--out-root", default=None, help="katalog główny wyników (domyślnie $EVIMO2_OUT albo ./out)")
    ap.add_argument("--glob", default=None, help="wzorzec plików dla --summary (domyślnie <out-root>/*/event_lag_check.json)")
    a = ap.parse_args()
    if a.summary:
        root = Path(a.out_root) if a.out_root else ec.default_out_root()
        return summary(a.glob or str(root / "*" / "event_lag_check.json"))
    if not a.seq_dir:
        ap.error("podaj katalog sekwencji albo użyj --summary")

    seq = Path(a.seq_dir.rstrip("/"))
    jpath = Path(a.json) if a.json else ec.find_result_file(seq.name, "gyro_vs_gt.json", a.out_root)
    if jpath is None or not jpath.is_file():
        where = ", ".join(str(c) for c in ec.result_dir_candidates(seq.name, a.out_root))
        sys.exit(f"Nie znaleziono gyro_vs_gt.json (szukano w: {where}).\n"
                 "Najpierw uruchom gyro_vs_gt.py dla tej sekwencji, ustaw EVIMO2_OUT/--out-root albo podaj --json.")
    out = Path(a.out) if a.out else jpath.parent
    gj = json.load(open(jpath))
    R, bias, d_gt = np.array(gj["R_imu_to_cam"]), np.array(gj["gyro_bias_imu_rad_s"]), float(gj["lag_s"])

    info = ec.load_info(seq)
    intr = ec.find_dict_with(info, ["fx", "fy", "cx", "cy"])
    imus = [(p, l) for p, l in ec.find_series(info).items() if ec.is_imu_series(l)]
    imus = [x for x in imus if "left" in x[0]] or imus
    cal_ts, gyro, _ = ec.imu_arrays(imus[0][1])
    cal_w = (gyro - bias) @ R.T                          # omega kamery (w czasie IMU)
    t, xy, _ = ec.load_events(seq)

    centers = pick_windows(cal_ts, cal_w, float(t[0]), float(t[-1]), a.windows, a.window_ms * 1e-3, a.min_omega)
    if not centers:
        sys.exit(f"Brak okien z |ω| >= {a.min_omega} rad/s - sekwencja ma za mało ruchu do tego testu.")
    lags = np.arange(a.lag_range[0], a.lag_range[1] + 1e-9, a.lag_step) * 1e-3
    print(f"Sekwencja: {seq.name} | d_GT = {d_gt * 1e3:+.1f} ms | {len(centers)} okien x {a.window_ms:.0f} ms, "
          f"{len(lags)} opóźnień ({lags[0] * 1e3:.0f}..{lags[-1] * 1e3:.0f} ms)")
    print(f"{'okno':>4} {'t [s]':>7} {'|ω|':>6} {'zdarzeń':>8} {'d_opt [ms]':>11} {'zysk kontrastu':>15}")

    wins, curves, per_win = [], [], []
    for j, tc in enumerate(centers):
        win = prep_window(tc, a.window_ms * 1e-3, t, xy, intr)
        v0 = iwe(*warp(win, cal_ts, cal_w, None)).var()
        c = np.array([iwe(*warp(win, cal_ts, cal_w, d)).var() for d in lags])
        wm = float(np.linalg.norm(cal_w[np.searchsorted(cal_ts, tc)]))
        p = dict(t_center=tc, omega=wm, n_events=win["n"], d_opt_ms=parabola_peak(lags, c) * 1e3, gain=float(c.max() / v0))
        print(f"{j:4d} {tc:7.3f} {wm:6.2f} {win['n']:8d} {p['d_opt_ms']:11.1f} {p['gain']:15.2f}")
        wins.append(win); curves.append(c); per_win.append(p)

    mean_curve = np.mean([c / c.max() for c in curves], axis=0)
    d_ev = parabola_peak(lags, mean_curve)
    opts = np.array([p["d_opt_ms"] for p in per_win])
    edge = d_ev <= lags[0] + 1e-9 or d_ev >= lags[-1] - 1e-9
    print(f"\nd_ev (optimum średniej krzywej) = {d_ev * 1e3:+.1f} ms | okna: mediana {np.median(opts):+.1f}, "
          f"std {opts.std():.1f} ms | d_GT = {d_gt * 1e3:+.1f} ms | różnica d_ev - d_GT = {(d_ev - d_gt) * 1e3:+.1f} ms")
    if edge:
        print("[UWAGA] optimum na brzegu przeszukiwanego zakresu - rozszerz --lag-range.")
    if opts.std() > 8:
        print("[UWAGA] duży rozrzut optimów między oknami: wynik niepewny (ruchome obiekty, mało ruchu, szum).")
    print("Rozstrzygnięcie wymaga porównania >= 2-3 sekwencji o różnym d_GT: python3 event_lag_check.py --summary")

    out.mkdir(parents=True, exist_ok=True)
    (out / "event_lag_check.json").write_text(json.dumps(dict(
        sequence=seq.name, d_gt_ms=d_gt * 1e3, d_ev_ms=d_ev * 1e3, d_ev_median_ms=float(np.median(opts)),
        d_ev_std_ms=float(opts.std()), n_windows=len(per_win), window_ms=a.window_ms,
        lag_grid_ms=(lags * 1e3).tolist(), mean_curve=mean_curve.tolist(), windows=per_win), indent=2))
    print(f"Zapisano {out / 'event_lag_check.json'}")
    if not a.no_plot:
        try:
            best = int(np.argmax([p["gain"] for p in per_win]))
            make_plot(out / "event_lag_check.png", seq.name, wins, lags, curves, d_gt, d_ev, cal_ts, cal_w, best, per_win)
            print(f"Zapisano {out / 'event_lag_check.png'}")
        except ImportError:
            print("Brak matplotlib - wykres pominięty")


if __name__ == "__main__":
    main()