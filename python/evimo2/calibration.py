"""Kalibracja czujników IMU względem kamery, estymacja biasu, ekstrynsek i opóźnienia."""
from typing import Any, Dict, Optional, Tuple

import numpy as np

from .geometry import R_to_rotvec, q2R, qconj, rotvec

# Stałe kalibracyjne zestawu sensorów Prophesee (z sekwencji EVIMO2 o wyraźnym ruchu)
R_IMU2CAM = np.diag([-1.0, -1.0, 1.0])
GYRO_BIAS = np.array([-0.0295, 0.0154, 0.0180])      # rad/s, w układzie IMU
LAG_CONST_S = -0.0157                                 # mediana opóźnienia IMU względem GT po 19 sekwencjach


def calibrated_gyro(gyro: np.ndarray, R: np.ndarray = R_IMU2CAM, bias: np.ndarray = GYRO_BIAS) -> np.ndarray:
    """Prędkość kątowa kamery [rad/s] z surowego żyroskopu (N, 3): R (omega - b)."""
    return (np.asarray(gyro) - bias) @ np.asarray(R).T


def gt_omega(ts: np.ndarray, q: np.ndarray, step: int, conv: str = "A") -> Tuple[np.ndarray, np.ndarray]:
    """Prędkość kątowa kamery w jej układzie [rad/s] w chwilach pośrednich dla konwencji A lub B."""
    dt = ts[step:] - ts[:-step]
    w = rotvec(q[:-step], q[step:]) if conv == "A" else rotvec(qconj(q[:-step]), qconj(q[step:]))
    ok = dt > 0
    return ((ts[:-step] + ts[step:]) / 2.0)[ok], w[ok] / dt[ok, None]


def kabsch(X: np.ndarray, Y: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Dopasowanie macierzy rotacji SO(3) i translacji: Y = R X + c (min. kwadratów)."""
    xm, ym = X.mean(0), Y.mean(0)
    Xc, Yc = X - xm, Y - ym
    U, S, Vt = np.linalg.svd(Xc.T @ Yc)
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    return R, ym - R @ xm, S, Yc - Xc @ R.T


def fit_at(tm: np.ndarray, Y: np.ndarray, ts_i: np.ndarray, X_i: np.ndarray, d: float, margin: float = 0.02) -> Dict[str, Any]:
    """Dopasowanie kabsch dla zadanego opóźnienia d: omega_gt(t) = R (omega_imu(t + d) - b)."""
    tq = tm + d
    ok = (tq >= ts_i[0] + margin) & (tq <= ts_i[-1] - margin)
    X = np.column_stack([np.interp(tq[ok], ts_i, X_i[:, k]) for k in range(3)])
    R, c, S, res = kabsch(X, Y[ok])
    Yc = Y[ok] - Y[ok].mean(0)
    rel_res = float(np.sqrt((res ** 2).sum() / (Yc ** 2).sum()))
    return dict(rel=rel_res, R=R, c=c, S=S, ok=ok, X=X)


def best_lag(tm: np.ndarray, Y: np.ndarray, ts_i: np.ndarray, X_i: np.ndarray, max_lag: float = 0.1) -> Tuple[float, np.ndarray, np.ndarray, bool]:
    """Wyszukuje optymalne opóźnienie d minimalizujące residuum dopasowania Kabscha."""
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


def bootstrap(tm: np.ndarray, X: np.ndarray, Y: np.ndarray, R_full: np.ndarray, n: int = 200, seed: int = 0, block_s: float = 0.25) -> Optional[Dict[str, Any]]:
    """Oszacowanie niepewności macierzy rotacji R metodą bootstrapu blokowego."""
    ids = np.floor((tm - tm[0]) / block_s).astype(int)
    blocks = [np.where(ids == b)[0] for b in np.unique(ids)]
    if len(blocks) < 4:
        return None
    rng = np.random.default_rng(seed)
    angs, vecs = [], []
    for _ in range(n):
        idx = np.concatenate([blocks[i] for i in rng.integers(0, len(blocks), len(blocks))])
        v = R_to_rotvec(kabsch(X[idx], Y[idx])[0] @ R_full.T)
        vecs.append(np.degrees(v))
        angs.append(np.degrees(np.linalg.norm(v)))
    return dict(axis_std_deg=np.std(vecs, axis=0).tolist(), angle_p95_deg=float(np.percentile(angs, 95)), n_blocks=len(blocks))


def unconstrained_det(X: np.ndarray, Y: np.ndarray) -> float:
    """Wyznacznik macierzy regresji liniowej bez ograniczenia ortogonalności."""
    M = np.linalg.lstsq(X - X.mean(0), Y - Y.mean(0), rcond=None)[0]
    return float(np.linalg.det(M))


def gravity_check(conv: str, ts_gt: np.ndarray, q_unif: np.ndarray, ts_i: np.ndarray, acc: np.ndarray, R: np.ndarray, d: float) -> Tuple[float, float, float]:
    """Test stałości kierunku wektora grawitacji w układzie świata."""
    tg = ts_i - d
    ok = (tg >= ts_gt[0]) & (tg <= ts_gt[-1])
    qi = np.column_stack([np.interp(tg[ok], ts_gt, q_unif[:, k]) for k in range(4)])
    qi /= np.linalg.norm(qi, axis=1, keepdims=True)
    Rs = q2R(qi)
    Rwc = Rs if conv == "A" else np.transpose(Rs, (0, 2, 1))
    a_w = np.einsum("nij,nj->ni", Rwc, acc[ok] @ R.T)
    u = a_w / np.linalg.norm(a_w, axis=1, keepdims=True)
    m = u.mean(0)
    norm_m = float(np.linalg.norm(m))
    cos_dev = np.clip(np.sum(u * (m / norm_m), axis=1), -1.0, 1.0)
    dev_p95 = float(np.degrees(np.arccos(cos_dev)).max())
    a_norm_mean = float(np.linalg.norm(a_w, axis=1).mean())
    return norm_m, dev_p95, a_norm_mean
