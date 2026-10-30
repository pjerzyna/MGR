#!/usr/bin/env python3
"""Żyroskop vs ground truth: opóźnienie IMU, rotacja IMU->kamera, bias, konwencja poz.

Użycie:
    python3 gyro_vs_gt.py <katalog_sekwencji> [--out DIR | --out-root DIR] [--step 2] [--max-lag 0.1] [--no-plot]

Metoda:
  1. Z kwaternionów kamery w full_trajectory liczymy prędkość kątową w układzie kamery
     w DWÓCH hipotezach konwencji poz:
        A: poza = T_świat<-kamera (R_wc),  B: poza = T_kamera<-świat (R_cw).
  2. Model:  omega_gt(t) = R (omega_imu(t + d) - b),  R in SO(3), d - opóźnienie, b - bias.
     Dla każdego d dopasowujemy R i b (Kabsch/SVD na danych scentrowanych), wybieramy d
     o najmniejszym niewyjaśnionym residuum (+ interpolacja paraboliczna).
     d > 0: ts IMU jest o d większe niż ts zdarzeń/GT dla tego samego ruchu
            (aby skorygować IMU: ts_imu - d).
  3. Niepewność R: bootstrap po blokach czasowych (0.25 s) + wartości osobliwe macierzy H.
  4. Konwencję rozstrzyga residuum (zła konwencja daje ω w innym układzie i znak,
     którego stała rotacja nie wchłonie) oraz niezależny test grawitacji: a_świat(t) =
     R_wc(t) R a_imu(t) powinno mieć stały kierunek.
  5. Tabela ruchu per obiekt (prędkości w układzie świata zakładając, że obiekty mają tę
     samą konwencję co kamera - ZAŁOŻENIE).
Zapisuje: <out>/gyro_vs_gt.png oraz <out>/gyro_vs_gt.json (R, opóźnienie, bias).
  Domyślnie <out> = $EVIMO2_OUT/<sceneXX_NN> (albo ./out/<sceneXX_NN>, jeśli zmienna nie jest ustawiona).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

import evimo2_common as ec
from survey_sequences import pose_series, qangle, qconj, qmul, rot_inv   # pomocnicze funkcje poz

MARGIN = 0.02      # [s] odcinek brzegowy pomijany w dopasowaniu
BLOCK = 0.25       # [s] długość bloku w bootstrapie


# ---------- geometria ----------
def q2R(q):
    w, x, y, z = q.T
    R = np.empty((len(q), 3, 3))
    R[:, 0, 0] = 1 - 2 * (y * y + z * z); R[:, 0, 1] = 2 * (x * y - z * w); R[:, 0, 2] = 2 * (x * z + y * w)
    R[:, 1, 0] = 2 * (x * y + z * w); R[:, 1, 1] = 1 - 2 * (x * x + z * z); R[:, 1, 2] = 2 * (y * z - x * w)
    R[:, 2, 0] = 2 * (x * z - y * w); R[:, 2, 1] = 2 * (y * z + x * w); R[:, 2, 2] = 1 - 2 * (x * x + y * y)
    return R


def rotvec(a, b):
    """Wektor obrotu log(a^* ⊗ b) dla tablic kwaternionów (N,4)."""
    r = qmul(qconj(a), b)
    r = r * np.where(r[:, :1] < 0, -1.0, 1.0)
    v = r[:, 1:]
    n = np.linalg.norm(v, axis=1)
    ang = 2 * np.arctan2(n, r[:, 0])
    f = np.where(n > 1e-12, ang / np.maximum(n, 1e-12), 2.0)
    return v * f[:, None]


def R_to_rotvec(R):
    ang = np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))
    v = 0.5 * np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    s = np.sin(ang)
    return v * (ang / s) if s > 1e-9 else v


def unify_hemisphere(q):
    q = q.copy()
    dots = np.sum(q[1:] * q[:-1], axis=1)
    q[1:] *= np.cumprod(np.where(dots < 0, -1.0, 1.0))[:, None]
    return q


def box(x, width):
    if width <= 1:
        return x
    k = np.ones(width) / width
    return np.column_stack([np.convolve(x[:, i], k, mode="same") for i in range(x.shape[1])])


# ---------- estymacja ----------
def gt_omega(ts, q, step, conv):
    """Prędkość kątowa kamery w jej układzie [rad/s] w chwilach pośrednich."""
    dt = ts[step:] - ts[:-step]
    w = rotvec(q[:-step], q[step:]) if conv == "A" else rotvec(qconj(q[:-step]), qconj(q[step:]))
    ok = dt > 0
    return ((ts[:-step] + ts[step:]) / 2)[ok], w[ok] / dt[ok, None]


def kabsch(X, Y):
    """y = R x + c, R in SO(3) (min. kwadratów). Zwraca R, c, wartości osobliwe H, residua."""
    xm, ym = X.mean(0), Y.mean(0)
    Xc, Yc = X - xm, Y - ym
    U, S, Vt = np.linalg.svd(Xc.T @ Yc)
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    return R, ym - R @ xm, S, Yc - Xc @ R.T


def fit_at(tm, Y, ts_i, X_i, d):
    tq = tm + d
    ok = (tq >= ts_i[0] + MARGIN) & (tq <= ts_i[-1] - MARGIN)
    X = np.column_stack([np.interp(tq[ok], ts_i, X_i[:, k]) for k in range(3)])
    R, c, S, res = kabsch(X, Y[ok])
    Yc = Y[ok] - Y[ok].mean(0)
    return dict(rel=float(np.sqrt((res ** 2).sum() / (Yc ** 2).sum())), R=R, c=c, S=S, ok=ok, X=X)


def best_lag(tm, Y, ts_i, X_i, max_lag):
    grid = np.arange(-max_lag, max_lag + 1e-9, 0.001)
    rels = np.array([fit_at(tm, Y, ts_i, X_i, d)["rel"] for d in grid])
    i = int(np.argmin(rels))
    d = grid[i]
    if 0 < i < len(grid) - 1:
        y0, y1, y2 = rels[i - 1:i + 2]
        den = y0 - 2 * y1 + y2
        if den > 0:
            d += 0.5 * (y0 - y2) / den * 0.001
    return d, grid, rels, bool(i in (0, len(grid) - 1))


def bootstrap(tm, X, Y, R_full, n=200, seed=0):
    ids = np.floor((tm - tm[0]) / BLOCK).astype(int)
    blocks = [np.where(ids == b)[0] for b in np.unique(ids)]
    if len(blocks) < 4:
        return None
    rng = np.random.default_rng(seed)
    angs, vecs = [], []
    for _ in range(n):
        idx = np.concatenate([blocks[i] for i in rng.integers(0, len(blocks), len(blocks))])
        v = R_to_rotvec(kabsch(X[idx], Y[idx])[0] @ R_full.T)
        vecs.append(np.degrees(v)); angs.append(np.degrees(np.linalg.norm(v)))
    return dict(axis_std_deg=np.std(vecs, axis=0).tolist(), angle_p95_deg=float(np.percentile(angs, 95)),
                n_blocks=len(blocks))


def unconstrained_det(X, Y):
    M = np.linalg.lstsq(X - X.mean(0), Y - Y.mean(0), rcond=None)[0]
    return float(np.linalg.det(M))


def gravity_check(conv, ts_gt, q_unif, ts_i, acc, R, d):
    tg = ts_i - d
    ok = (tg >= ts_gt[0]) & (tg <= ts_gt[-1])
    qi = np.column_stack([np.interp(tg[ok], ts_gt, q_unif[:, k]) for k in range(4)])
    qi /= np.linalg.norm(qi, axis=1, keepdims=True)
    Rs = q2R(qi)
    Rwc = Rs if conv == "A" else np.transpose(Rs, (0, 2, 1))
    a_w = np.einsum("nij,nj->ni", Rwc, acc[ok] @ R.T)
    u = a_w / np.linalg.norm(a_w, axis=1, keepdims=True)
    m = u.mean(0)
    m /= np.linalg.norm(m)
    return tg[ok], np.degrees(np.arccos(np.clip(u @ m, -1, 1)))


def analyse(conv, ts, q, ts_i, gyro_s, acc_s, step, max_lag):
    tm, Y = gt_omega(ts, q, step, conv)
    d, grid, rels, at_edge = best_lag(tm, Y, ts_i, gyro_s, max_lag)
    fit = fit_at(tm, Y, ts_i, gyro_s, d)
    fit0 = fit_at(tm, Y, ts_i, gyro_s, 0.0)
    ok = fit["ok"]
    R, c = fit["R"], fit["c"]
    S = fit["S"]
    return dict(conv=conv, tm=tm, Y=Y, d=d, grid=grid, rels=rels, at_edge=at_edge, fit=fit,
                rel0=fit0["rel"], R=R, c=c, bias=(-R.T @ c), S=S,
                boot=bootstrap(tm[ok], fit["X"], Y[ok], R),
                det_free=unconstrained_det(fit["X"], Y[ok]),
                grav=gravity_check(conv, ts, unify_hemisphere(q), ts_i, acc_s, R, d))


# ---------- raport ----------
def report(a):
    np.set_printoptions(precision=4, suppress=True)
    name = {"A": "A: poza = T_świat<-kamera", "B": "B: poza = T_kamera<-świat"}[a["conv"]]
    print(f"\n== Konwencja {name} ==")
    print(f"  opóźnienie d = {a['d'] * 1e3:+.1f} ms" + ("  [UWAGA: na brzegu zakresu przeszukiwania]" if a["at_edge"] else ""))
    print(f"  niewyjaśnione residuum: {a['fit']['rel'] * 100:.1f}% (przy d = 0: {a['rel0'] * 100:.1f}%)")
    S = a["S"] / a["S"][0]
    print(f"  wartości osobliwe H (znorm.): {S[0]:.2f} {S[1]:.2f} {S[2]:.2f}", end="")
    print("  [UWAGA: ruch prawie wokół jednej osi, R słabo wyznaczona]" if S[1] < 0.15 else "")
    print(f"  det swobodnego dopasowania liniowego: {a['det_free']:+.2f} (ujemny = układ IMU wygląda na lustrzany względem kamery)")
    rv = np.degrees(R_to_rotvec(a["R"]))
    print("  R_imu->cam =")
    for row in a["R"]:
        print("    [" + "  ".join(f"{x:8.4f}" for x in row) + " ]")
    print(f"  wektor obrotu R: {rv} deg (|.| = {np.linalg.norm(rv):.1f} deg)")
    print(f"  bias żyroskopu (układ IMU): {a['bias']} rad/s")
    if a["boot"]:
        b = a["boot"]
        print(f"  niepewność R (bootstrap, {b['n_blocks']} bloków): std osi kamery = "
              f"{np.round(b['axis_std_deg'], 2)} deg, p95 kąta = {b['angle_p95_deg']:.2f} deg")
        if max(b["axis_std_deg"]) > 2:
            print("  [UWAGA: niepewność jednej z osi > 2 deg]")
    else:
        print("  bootstrap: za mało bloków czasowych (< 4), pominięty")
    g = a["grav"][1]
    print(f"  test grawitacji: odchylenie kierunku a_świat od średniej: mediana {np.median(g):.1f} deg, p95 {np.percentile(g, 95):.1f} deg")


def verdict(A, B):
    ra, rb = A["fit"]["rel"], B["fit"]["rel"]
    win, lose = (A, B) if ra <= rb else (B, A)
    ratio = max(ra, rb) / max(min(ra, rb), 1e-12)
    gw, gl = np.median(win["grav"][1]), np.median(lose["grav"][1])
    print("\n== Wniosek ==")
    strength = "wyraźna" if ratio >= 1.3 else "słaba" if ratio >= 1.1 else "BRAK (konwencje nierozróżnialne)"
    print(f"  lepsze dopasowanie ω: konwencja {win['conv']} (residuum {win['fit']['rel'] * 100:.1f}% vs "
          f"{lose['fit']['rel'] * 100:.1f}%, stosunek {ratio:.2f}) - przewaga: {strength}")
    print(f"  test grawitacji (mediana odchylenia): {win['conv']} {gw:.1f} deg vs {lose['conv']} {gl:.1f} deg"
          f" -> {'zgodny z ω' if gw <= gl else 'SPRZECZNY z ω'}")
    if win["fit"]["rel"] > 0.35:
        print("  [UWAGA] nawet najlepsze dopasowanie wyjaśnia mało wariancji: wynik niepewny "
              "(za mały ruch, szum, albo inny problem z danymi).")
    return win


def object_table(traj, conv, step):
    ids = sorted({k for e in traj for k in e if isinstance(k, str) and k.isdigit()}, key=int)
    rows = []
    for key in ["cam"] + ids:
        s = pose_series(traj, key)
        if s is None:
            continue
        ts, q, t = s
        st = max(1, min(step, len(ts) - 1))
        p = t if conv == "A" else -rot_inv(q, t)
        dt = ts[st:] - ts[:-st]
        ok = dt > 0
        v = np.linalg.norm(p[st:] - p[:-st], axis=1)[ok] / dt[ok]
        w = qangle(q[:-st], q[st:])[ok] / dt[ok]
        rows.append((key, np.median(v), np.percentile(v, 95), v.max(), np.percentile(w, 95), float((v > 0.05).mean())))
    print("\n== Ruch w układzie świata (obiekty: ZAŁOŻENIE tej samej konwencji co kamera) ==")
    print(f"  {'id':>4}  {'v med':>7}  {'v p95':>7}  {'v max':>7}  {'ω p95':>7}  {'czas v>5cm/s':>12}")
    print(f"  {'':>4}  {'[m/s]':>7}  {'[m/s]':>7}  {'[m/s]':>7}  {'[rad/s]':>7}")
    for key, med, p95, mx, wp, fr in sorted(rows, key=lambda r: (r[0] != "cam", -r[2])):
        print(f"  {('cam' if key == 'cam' else key):>4}  {med:7.3f}  {p95:7.3f}  {mx:7.3f}  {wp:7.3f}  {fr * 100:11.0f}%")
    return rows


# ---------- wykres ----------
def make_plot(path, seq_name, A, B, win, ts_i, gyro, step):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(15, 9.5))
    gs = fig.add_gridspec(3, 2, hspace=0.38, wspace=0.18)
    f = win["fit"]
    pred = f["X"] @ win["R"].T + win["c"]
    tt = win["tm"][f["ok"]]
    for k, lab in enumerate("xyz"):
        ax = fig.add_subplot(gs[k, 0])
        ax.plot(tt, win["Y"][f["ok"], k], "k", lw=1.4, label="GT (poza)")
        ax.plot(tt, pred[:, k], "r", lw=1.0, label=r"IMU: $R(\omega_{imu}-b)$")
        ax.set_ylabel(f"ω_{lab} kamery [rad/s]")
        ax.grid(alpha=0.3)
        if k == 0:
            ax.set_title(f"{seq_name} - konwencja {win['conv']}, d = {win['d'] * 1e3:+.1f} ms, "
                         f"residuum {f['rel'] * 100:.1f}%")
            ax.legend(loc="upper right", ncol=2)
        if k == 2:
            ax.set_xlabel("czas [s]")

    ax = fig.add_subplot(gs[0, 1])
    for a, col in ((A, "tab:blue"), (B, "tab:orange")):
        ax.plot(a["grid"] * 1e3, a["rels"] * 100, color=col, label=f"konwencja {a['conv']}")
        ax.plot(a["d"] * 1e3, a["fit"]["rel"] * 100, "o", color=col)
    ax.set_xlabel("opóźnienie d [ms]"); ax.set_ylabel("niewyjaśnione residuum [%]")
    ax.set_title("Residuum vs opóźnienie"); ax.legend(); ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1, 1])
    tg = win["tm"]
    ax.plot(tg, np.linalg.norm(win["Y"], axis=1), "k", lw=1.4, label="|ω| GT")
    ax.plot(ts_i - win["d"], np.linalg.norm(gyro, axis=1), "r", lw=0.8, label="|ω| IMU (surowe, ts - d)")
    ax.set_xlim(tg[0], tg[-1]); ax.set_xlabel("czas [s]"); ax.set_ylabel("|ω| [rad/s]")
    ax.set_title("Moduł ω (nie zależy od konwencji ani od R)"); ax.legend(); ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[2, 1])
    for a, col in ((A, "tab:blue"), (B, "tab:orange")):
        ax.plot(a["grav"][0], a["grav"][1], color=col, lw=0.8, label=f"konwencja {a['conv']}")
    ax.set_xlabel("czas [s]"); ax.set_ylabel("odchylenie [deg]")
    ax.set_title("Test grawitacji: kąt a_świat od średniego kierunku"); ax.legend(); ax.grid(alpha=0.3)
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("seq_dir")
    ap.add_argument("--out", default=None, help="katalog wyjściowy tej sekwencji")
    ap.add_argument("--out-root", default=None, help="katalog główny wyników (domyślnie $EVIMO2_OUT albo ./out)")
    ap.add_argument("--step", type=int, default=2, help="odstęp próbek GT do różnicowania (2 = 10 ms)")
    ap.add_argument("--max-lag", type=float, default=0.1, help="zakres przeszukiwania opóźnienia [s]")
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()

    seq = Path(a.seq_dir.rstrip("/"))
    out = Path(a.out) if a.out else ec.result_dir(seq.name, a.out_root)
    info = ec.load_info(seq)
    if info is None:
        sys.exit("brak dataset_info.npz")
    traj = info["meta"].get("full_trajectory") or []
    cam = pose_series(traj, "cam")
    if cam is None:
        sys.exit("brak poz kamery w full_trajectory")
    ts, q, _ = cam
    imus = [(p, l) for p, l in ec.find_series(info).items() if ec.is_imu_series(l)]
    imus = [x for x in imus if "left" in x[0]] or imus
    if not imus:
        sys.exit("brak IMU w meta")
    ts_i, gyro, acc = ec.imu_arrays(imus[0][1])

    st = max(1, min(a.step, len(ts) - 1))
    width = max(1, int(round(st * np.median(np.diff(ts)) / np.median(np.diff(ts_i)))))
    gyro_s = box(gyro, width)
    acc_s = box(acc, 21)

    print(f"Sekwencja: {seq.name}")
    print(f"GT: {len(ts)} póz kamery ({1 / np.median(np.diff(ts)):.0f} Hz), IMU: {len(ts_i)} próbek "
          f"({1 / np.median(np.diff(ts_i)):.0f} Hz); różnicowanie co {st} próbek GT "
          f"({st * np.median(np.diff(ts)) * 1e3:.0f} ms), wygładzanie IMU: {width} próbek")

    A = analyse("A", ts, q, ts_i, gyro_s, acc_s, st, a.max_lag)
    B = analyse("B", ts, q, ts_i, gyro_s, acc_s, st, a.max_lag)
    wmag = np.linalg.norm(A["Y"], axis=1)
    print(f"|ω| z GT: p95 {np.percentile(wmag, 95):.3f}, max {wmag.max():.3f} rad/s")
    report(A)
    report(B)
    win = verdict(A, B)
    object_table(traj, win["conv"], st)

    out.mkdir(parents=True, exist_ok=True)
    lose = B if win is A else A
    res = dict(sequence=seq.name, convention=win["conv"], lag_s=float(win["d"]),
               R_imu_to_cam=win["R"].tolist(), gyro_bias_imu_rad_s=win["bias"].tolist(),
               residual_rel=float(win["fit"]["rel"]), residual_rel_other_convention=float(lose["fit"]["rel"]),
               bootstrap=win["boot"], step_gt_samples=st)
    (out / "gyro_vs_gt.json").write_text(json.dumps(res, indent=2))
    print(f"\nZapisano {out / 'gyro_vs_gt.json'}")
    if not a.no_plot:
        try:
            make_plot(out / "gyro_vs_gt.png", seq.name, A, B, win, ts_i, gyro, st)
            print(f"Zapisano {out / 'gyro_vs_gt.png'}")
        except ImportError:
            print("Brak matplotlib - wykres pominięty (pip install matplotlib)")


if __name__ == "__main__":
    main()