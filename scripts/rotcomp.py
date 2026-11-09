#!/usr/bin/env python3
"""Kompensacja rotacji zdarzeń na podstawie żyroskopu (prototyp v0; referencja do portu na C++).

Model: czysty obrót kamery. Dla zdarzenia o czasie t_i i piksela (x, y):
  1. korekcja zniekształceń (radtan) -> promień (u, v, 1) w układzie kamery w chwili t_i,
  2. orientacja kamery R(t) z całkowania prędkości kątowej  omega_cam(t) = R_imu->cam (omega_imu(t + d) - b)
     (siatka 1 ms, w obrębie kroku obrót dokładny z własnym znacznikiem czasu zdarzenia),
  3. obrót promienia do chwili odniesienia: v_ref = R(t_ref)^T R(t_i) v_i, rzut z powrotem na obraz.
Parametry układu czujników (R_IMU2CAM, GYRO_BIAS) są stałymi kalibracyjnymi zestawu - NIE zależą od GT
w trakcie przetwarzania. Opóźnienie d jest parametrem wejściowym (zmienia się między nagraniami).
Funkcje gt_warp / nlerp służą wyłącznie do ewaluacji (orientacja z ground truth).
"""
import numpy as np

# Stałe kalibracyjne zestawu (z gyro_vs_gt.py na 4 sekwencjach z wyraźnym ruchem; odchylenia R < 1 deg,
# bias +-0,001 rad/s między sekwencjami). Zapisane tu raz, jako kalibracja czujników.
R_IMU2CAM = np.diag([-1.0, -1.0, 1.0])
GYRO_BIAS = np.array([-0.0295, 0.0154, 0.0180])      # rad/s, w układzie IMU
LAG_CONST_S = -0.0157                                 # mediana opóźnienia IMU względem GT po 19 sekwencjach


def calibrated_gyro(gyro, R=R_IMU2CAM, bias=GYRO_BIAS):
    """Prędkość kątowa kamery [rad/s] z surowego żyroskopu (N, 3): R (omega - b)."""
    return (np.asarray(gyro) - bias) @ np.asarray(R).T


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
    return dict(tc=tc, te=te, pix=pix, bearing=np.column_stack([bx, by, np.ones(len(te))]),
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

# ---------- orientacja z ground truth (tylko ewaluacja) ----------
def nlerp(ts, q, t):
    """Interpolacja kwaternionów (nlerp) w chwilach t; q musi być ciągłe półsferowo."""
    t = np.atleast_1d(np.asarray(t, dtype=float))
    i = np.clip(np.searchsorted(ts, t, side="right") - 1, 0, len(ts) - 2)
    a = np.clip((t - ts[i]) / (ts[i + 1] - ts[i]), 0.0, 1.0)[:, None]
    out = (1 - a) * q[i] + a * q[i + 1]
    return out / np.linalg.norm(out, axis=1, keepdims=True)


def gt_warp(win, ts_gt, q_gt):
    """Jak warp(), ale orientacja z poz kamery GT (konwencja A: q = R_wc). Tylko do porównań."""
    from gyro_vs_gt import q2R
    fx, fy, cx, cy = (float(win["intr"][k]) for k in ("fx", "fy", "cx", "cy"))
    Ri = q2R(nlerp(ts_gt, q_gt, win["te"]))
    Rr = q2R(nlerp(ts_gt, q_gt, [win["grid"][win["kref"]]]))[0]
    M = np.einsum("ji,njk->nik", Rr, Ri)                           # R_wc(t_ref)^T R_wc(t_i)
    b2 = np.einsum("nij,nj->ni", M, win["bearing"])
    return fx * b2[:, 0] / b2[:, 2] + cx, fy * b2[:, 1] / b2[:, 2] + cy
