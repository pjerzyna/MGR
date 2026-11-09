"""Geometria widzenia, algebra kwaternionów, SO(3) i zniekształcenia radtan."""
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

I3 = np.eye(3)


# ---------- Kwaterniony i SO(3) ----------
def qconj(q: np.ndarray) -> np.ndarray:
    """Koniugacja kwaternionu q = [w, x, y, z] -> [w, -x, -y, -z]."""
    return q * np.array([1.0, -1.0, -1.0, -1.0])


def qmul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Iloczyn kwaternionów Hamiltona a ⊗ b dla tablic (N, 4)."""
    aw, ax, ay, az = a.T
    bw, bx, by, bz = b.T
    return np.stack([
        aw * bw - ax * bx - ay * by - az * bz,
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
    ], axis=1)


def qangle(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Kąt względnego obrotu [rad] między kwaternionami a i b."""
    r = qmul(qconj(a), b)
    return 2.0 * np.arctan2(np.linalg.norm(r[:, 1:], axis=1), np.abs(r[:, 0]))


def qrot(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Obrót wektorów v (N, 3) przez kwaterniony q (N, 4): q ⊗ (0, v) ⊗ q*."""
    p = np.concatenate([np.zeros((len(v), 1)), v], axis=1)
    return qmul(qmul(q, p), qconj(q))[:, 1:]


def rot_inv(q: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Obrót odwrotny R^T t policzony kwaternionami: q* ⊗ (0, t) ⊗ q."""
    p = np.concatenate([np.zeros((len(t), 1)), t], axis=1)
    return qmul(qmul(qconj(q), p), q)[:, 1:]


def q2R(q: np.ndarray) -> np.ndarray:
    """Konwersja tablicy kwaternionów (N, 4) do macierzy rotacji SO(3) (N, 3, 3)."""
    w, x, y, z = q.T
    R = np.empty((len(q), 3, 3))
    R[:, 0, 0] = 1 - 2 * (y * y + z * z); R[:, 0, 1] = 2 * (x * y - z * w); R[:, 0, 2] = 2 * (x * z + y * w)
    R[:, 1, 0] = 2 * (x * y + z * w);     R[:, 1, 1] = 1 - 2 * (x * x + z * z); R[:, 1, 2] = 2 * (y * z - x * w)
    R[:, 2, 0] = 2 * (x * z - y * w);     R[:, 2, 1] = 2 * (y * z + x * w);     R[:, 2, 2] = 1 - 2 * (x * x + y * y)
    return R


def rotvec(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Wektor obrotu log(a^* ⊗ b) dla tablic kwaternionów (N, 4)."""
    r = qmul(qconj(a), b)
    r = r * np.where(r[:, :1] < 0, -1.0, 1.0)
    v = r[:, 1:]
    n = np.linalg.norm(v, axis=1)
    ang = 2.0 * np.arctan2(n, r[:, 0])
    f = np.where(n > 1e-12, ang / np.maximum(n, 1e-12), 2.0)
    return v * f[:, None]


def R_to_rotvec(R: np.ndarray) -> np.ndarray:
    """Konwersja pojedynczej macierzy rotacji SO(3) (3, 3) na wektor obrotu axis-angle."""
    ang = np.arccos(np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0))
    v = 0.5 * np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    s = np.sin(ang)
    return v * (ang / s) if s > 1e-9 else v


def unify_hemisphere(q: np.ndarray) -> np.ndarray:
    """Uzgadnia półsferę kwaternionów w szeregu czasowym (ciągłość q_i dot q_{i+1} >= 0)."""
    q_out = q.copy()
    dots = np.sum(q_out[1:] * q_out[:-1], axis=1)
    q_out[1:] *= np.cumprod(np.where(dots < 0, -1.0, 1.0))[:, None]
    return q_out


def nlerp(ts: np.ndarray, q: np.ndarray, t: Union[float, np.ndarray, List[float]]) -> np.ndarray:
    """Interpolacja kwaternionów (nlerp) w chwilach t; q musi być ciągłe półsferowo."""
    t_arr = np.atleast_1d(np.asarray(t, dtype=float))
    i = np.clip(np.searchsorted(ts, t_arr, side="right") - 1, 0, len(ts) - 2)
    dt = ts[i + 1] - ts[i]
    a = np.clip((t_arr - ts[i]) / np.maximum(dt, 1e-12), 0.0, 1.0)[:, None]
    out = (1.0 - a) * q[i] + a * q[i + 1]
    return out / np.linalg.norm(out, axis=1, keepdims=True)


def rodrigues_batch(rv: np.ndarray) -> np.ndarray:
    """Wzór Rodrigues dla paczki wektorów obrotu (N, 3) -> macierze SO(3) (N, 3, 3)."""
    th = np.linalg.norm(rv, axis=1)
    small = th < 1e-6
    ths = np.where(small, 1.0, th)
    a = np.where(small, 1.0 - th ** 2 / 6.0, np.sin(ths) / ths)
    b = np.where(small, 0.5 - th ** 2 / 24.0, (1.0 - np.cos(ths)) / ths ** 2)
    K = np.zeros((len(rv), 3, 3))
    K[:, 0, 1], K[:, 0, 2], K[:, 1, 0] = -rv[:, 2], rv[:, 1], rv[:, 2]
    K[:, 1, 2], K[:, 2, 0], K[:, 2, 1] = -rv[:, 0], -rv[:, 1], rv[:, 0]
    return I3 + a[:, None, None] * K + b[:, None, None] * (K @ K)


# ---------- Model kamery i dystorsja RadTan ----------
def undistort(px: np.ndarray, py: np.ndarray, intr: Dict[str, Any], iters: int = 8) -> Tuple[np.ndarray, np.ndarray]:
    """Inwersja modelu radtan OpenCV (Brown-Conrady) metodą punktu stałego.

    Zwraca znormalizowane współrzędne promienia (x, y) na płaszczyźnie rzutowej z=1.
    """
    fx, fy, cx, cy = (float(intr[k]) for k in ("fx", "fy", "cx", "cy"))
    k1, k2, k3, p1, p2 = (float(intr.get(k, 0.0)) for k in ("k1", "k2", "k3", "p1", "p2"))
    xd, yd = (px - cx) / fx, (py - cy) / fy
    x, y = xd.copy(), yd.copy()
    for _ in range(iters):
        r2 = x * x + y * y
        rad = 1.0 + k1 * r2 + k2 * r2 ** 2 + k3 * r2 ** 3
        x = (xd - (2.0 * p1 * x * y + p2 * (r2 + 2.0 * x * x))) / rad
        y = (yd - (p1 * (r2 + 2.0 * y * y) + 2.0 * p2 * x * y)) / rad
    return x, y


# ---------- Parsowanie poz z full_trajectory ----------
def entry_pose(ent: Any) -> Optional[Tuple[float, List[float], List[float]]]:
    """Zwraca (ts, [w,x,y,z], [x,y,z]) z wpisu trajektorii."""
    pos = ent.get("pos") if isinstance(ent, dict) else None
    if isinstance(pos, dict) and "q" in pos and "t" in pos:
        return float(ent["ts"]), [pos["q"][k] for k in "wxyz"], [pos["t"][k] for k in "xyz"]
    return None


def pose_series(traj: List[Dict[str, Any]], key: str) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Wyciąga szereg czasowy pozy (ts, q, t) dla danego klucza ('cam', '9', '22'...)."""
    rows = []
    for e in traj:
        ent = e.get(key)
        pos = entry_pose(ent)
        if pos is not None:
            rows.append(pos)
    if len(rows) < 2:
        return None
    rows.sort(key=lambda r: r[0])
    ts = np.array([r[0] for r in rows], dtype=np.float64)
    q = np.array([r[1] for r in rows], dtype=np.float64)
    t = np.array([r[2] for r in rows], dtype=np.float64)
    q_norm = q / np.linalg.norm(q, axis=1, keepdims=True)
    return ts, q_norm, t


def motion(ts: np.ndarray, q: np.ndarray, t: np.ndarray, step: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Wyznacza w [rad/s] oraz prędkości kamery w dwóch hipotezach konwencji."""
    st = max(1, min(step, len(ts) - 1))
    dt = ts[st:] - ts[:-st]
    ok = dt > 0
    w = qangle(q[:-st], q[st:])[ok] / dt[ok]
    disp = np.linalg.norm(t[st:] - t[:-st], axis=1)[ok] / dt[ok]
    t_alt = rot_inv(q, t)
    disp_alt = np.linalg.norm(t_alt[st:] - t_alt[:-st], axis=1)[ok] / dt[ok]
    return w, disp, disp_alt
