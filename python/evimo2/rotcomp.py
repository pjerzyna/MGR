"""Kompensacja rotacji zdarzeń na podstawie prędkości kątowej żyroskopu (SO(3), RadTan, IWE)."""
from typing import Any, Dict, Optional, Tuple

import numpy as np

from .calibration import GYRO_BIAS, LAG_CONST_S, R_IMU2CAM, calibrated_gyro
from .geometry import I3, nlerp, q2R, rodrigues_batch, undistort

W_DEFAULT, H_DEFAULT = 640, 480


def prep_window(tc: float, wlen: float, t: np.ndarray, xy: np.ndarray, intr: Dict[str, Any]) -> Dict[str, Any]:
    """Przygotowuje okno zdarzeń wokół chwili tc: odwraca dystorsję, tworzy promienie i siatkę całkowania."""
    i0, i1 = np.searchsorted(t, [tc - wlen / 2.0, tc + wlen / 2.0])
    te = np.asarray(t[i0:i1], dtype=np.float64)
    pix = np.asarray(xy[i0:i1], dtype=np.float64)
    bx, by = undistort(pix[:, 0], pix[:, 1], intr)
    nst = max(2, int(round(wlen / 0.001)))
    return dict(
        tc=tc,
        te=te,
        pix=pix,
        bearing=np.column_stack([bx, by, np.ones(len(te))]),
        grid=np.linspace(tc - wlen / 2.0, tc + wlen / 2.0, nst + 1),
        kref=nst // 2,
        n=len(te),
        intr=intr,
    )


def rot_matrices(win: Dict[str, Any], cal_ts: np.ndarray, cal_w: np.ndarray, d: float) -> np.ndarray:
    """Macierze obrotu M_i = R(t_ref)^T R(t_i) (N, 3, 3) z całkowania prędkości kątowej (opóźnienie IMU d)."""
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
    return np.einsum("ji,njk->nik", Rg[win["kref"]], Ri)            # R_wc(t_ref)^T R_wc(t_i)


def warp(win: Dict[str, Any], cal_ts: np.ndarray, cal_w: np.ndarray, d: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray]:
    """Pozycje pikseli zdarzeń po obrocie do chwili środkowej okna; d=None -> bez kompensacji."""
    fx, fy, cx, cy = (float(win["intr"][k]) for k in ("fx", "fy", "cx", "cy"))
    b = win["bearing"]
    if d is None:
        return fx * b[:, 0] + cx, fy * b[:, 1] + cy

    b2 = np.einsum("nij,nj->ni", rot_matrices(win, cal_ts, cal_w, d), b)
    return fx * b2[:, 0] / b2[:, 2] + cx, fy * b2[:, 1] / b2[:, 2] + cy


def iwe(xs: np.ndarray, ys: np.ndarray, W: int = W_DEFAULT, H: int = H_DEFAULT) -> np.ndarray:
    """Tworzy obraz Image of Warped Events (IWE) metodą biliniowej akumulacji (splattingu)."""
    x0, y0 = np.floor(xs), np.floor(ys)
    ax, ay = xs - x0, ys - y0
    x0, y0 = x0.astype(np.int64), y0.astype(np.int64)
    ok = (x0 >= 0) & (x0 < W - 1) & (y0 >= 0) & (y0 < H - 1)
    x0, y0, ax, ay = x0[ok], y0[ok], ax[ok], ay[ok]
    idx = y0 * W + x0
    n = W * H
    img = (
        np.bincount(idx, (1.0 - ax) * (1.0 - ay), n)
        + np.bincount(idx + 1, ax * (1.0 - ay), n)
        + np.bincount(idx + W, (1.0 - ax) * ay, n)
        + np.bincount(idx + W + 1, ax * ay, n)
    )
    return img.reshape(H, W)


def gt_rot_matrices(win: Dict[str, Any], ts_gt: np.ndarray, q_gt: np.ndarray) -> np.ndarray:
    """M_i = R_wc(t_ref)^T R_wc(t_i) z orientacji poz kamery GT (konwencja A: q = R_wc). Tylko do ewaluacji."""
    Ri = q2R(nlerp(ts_gt, q_gt, win["te"]))
    Rr = q2R(nlerp(ts_gt, q_gt, [win["grid"][win["kref"]]]))[0]
    return np.einsum("ji,njk->nik", Rr, Ri)


def gt_warp(win: Dict[str, Any], ts_gt: np.ndarray, q_gt: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Warping orientacją z poz kamery GT (wzorzec oracle do ewaluacji)."""
    fx, fy, cx, cy = (float(win["intr"][k]) for k in ("fx", "fy", "cx", "cy"))
    b2 = np.einsum("nij,nj->ni", gt_rot_matrices(win, ts_gt, q_gt), win["bearing"])
    return fx * b2[:, 0] / b2[:, 2] + cx, fy * b2[:, 1] / b2[:, 2] + cy


def warp_translation(
    win: Dict[str, Any],
    M: np.ndarray,
    ts_gt: np.ndarray,
    q_gt: np.ndarray,
    p_gt: np.ndarray,
    depth_z: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Kompensacja rotacji i translacji (ORACLE: ruch kamery z poz GT, głębia z map GT).

    Punkt statyczny widziany w chwili t_i jako P_i = Z_i * (u, v, 1) (Z - głębia wzdłuż osi optycznej) ma
    w układzie kamery w chwili odniesienia t_ref współrzędne
        P_ref = M_i P_i + R_wc(t_ref)^T (p_c(t_i) - p_c(t_ref)),     M_i = R_wc(t_ref)^T R_wc(t_i),
    gdzie p_c to pozycja kamery w świecie (konwencja A: t z pozy kamery). M może pochodzić z IMU
    (rot_matrices) albo z GT (gt_rot_matrices). depth_z: głębia [m] każdego zdarzenia (N,).
    Obiekty poruszające się niezależnie NIE są dopasowane do tego modelu (to właśnie ich sygnał).
    """
    fx, fy, cx, cy = (float(win["intr"][k]) for k in ("fx", "fy", "cx", "cy"))
    t_ref = win["grid"][win["kref"]]
    Rr = q2R(nlerp(ts_gt, q_gt, [t_ref]))[0]
    p_i = np.column_stack([np.interp(win["te"], ts_gt, p_gt[:, k]) for k in range(3)])
    p_r = np.array([np.interp(t_ref, ts_gt, p_gt[:, k]) for k in range(3)])
    T = (p_i - p_r) @ Rr                                              # R_ref^T (p_i - p_ref)
    P = np.einsum("nij,nj->ni", M, depth_z[:, None] * win["bearing"]) + T
    ok = P[:, 2] > 1e-3
    z = np.where(ok, P[:, 2], 1.0)
    return np.where(ok, fx * P[:, 0] / z + cx, -1e6), np.where(ok, fy * P[:, 1] / z + cy, -1e6)
