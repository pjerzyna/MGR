"""Kinematyka obiektów w świecie, składanie poz (T_w<-o = T_w<-c * T_c<-o) i etykietowanie ruchu."""
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .geometry import entry_pose, qangle, qmul, qrot

CLASSES = ["S", "R*S", "S*R", "Rt*S", "S*Rt"]
MIN_CAM_ROT_DEG = 10.0
V_THRESHOLDS = [0.05, 0.1, 0.3, 1.0]
V_EDGES = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.3, 1.0]
W_EDGES = [0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]


def load_tracks(info: Dict[str, Any]) -> Tuple[Dict[str, np.ndarray], Dict[str, Dict[str, np.ndarray]]]:
    """Kamera i obiekty z full_trajectory, tylko ze wspólnych próbek czasowych."""
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
    camd = dict(
        ts=np.array([c[0] for c in cam], dtype=np.float64),
        q=norm(np.array([c[1] for c in cam], float)),
        t=np.array([c[2] for c in cam], float),
    )
    tracks = {}
    for k, r in rows.items():
        if len(r) < 3:
            continue
        r.sort(key=lambda x: x[0])
        tracks[k] = dict(
            ts=arr(r, 0, 1)[:, 0],
            qc=norm(arr(r, 1, 4)),
            tc=arr(r, 2, 3),
            qs=norm(arr(r, 3, 4)),
            ts_=arr(r, 4, 3),
        )
    return camd, tracks


def invariants(qc: np.ndarray, qs: np.ndarray) -> Dict[str, np.ndarray]:
    """Niezmienniki rotacji dla obiektu nieruchomego w świecie."""
    from .geometry import qconj
    return {
        "S": qs,
        "R*S": qmul(qc, qs),
        "S*R": qmul(qs, qc),
        "Rt*S": qmul(qconj(qc), qs),
        "S*Rt": qmul(qs, qconj(qc)),
    }


def dev_deg(q: np.ndarray) -> np.ndarray:
    """Odchylenie kątowe względem pierwszej próbki szeregu [deg]."""
    return np.degrees(qangle(np.repeat(q[:1], len(q), 0), q))


def world_pose(cls: str, qc: np.ndarray, tc: np.ndarray, qs: np.ndarray, ts_: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Złożenie pozy obiektu do układu świata w zależności od klasy konwencji."""
    from .geometry import qconj
    if cls == "R*S":
        return qmul(qc, qs), qrot(qc, ts_) + tc
    if cls == "S*Rt":
        q = qmul(qc, qconj(qs))
        return q, tc - qrot(q, ts_)
    return qs, ts_                       # 'S' i klasy przyklejone do kamery


def kinematics(tracks: Dict[str, Dict[str, np.ndarray]], cls: str = "R*S", half_s: float = 0.025) -> Dict[str, Dict[str, np.ndarray]]:
    """Prędkość liniowa v [m/s] i kątowa ω [rad/s] każdego obiektu w świecie (różnica centralna)."""
    out = {}
    for k, t in tracks.items():
        q, p = world_pose(cls, t["qc"], t["tc"], t["qs"], t["ts_"])
        ts = t["ts"]
        n = len(ts)
        med_dt = max(float(np.median(np.diff(ts))), 1e-6) if n > 1 else 0.005
        s = max(1, int(round(half_s / med_dt)))
        i = np.arange(n)
        a, b = np.clip(i - s, 0, n - 1), np.clip(i + s, 0, n - 1)
        dt = ts[b] - ts[a]
        ok = dt > 0
        v, w = np.zeros(n), np.zeros(n)
        v[ok] = np.linalg.norm(p[b] - p[a], axis=1)[ok] / dt[ok]
        w[ok] = qangle(q[a], q[b])[ok] / dt[ok]
        out[k] = dict(ts=ts, v=v, w=w)
    return out


def motion_stats(ts: np.ndarray, q: np.ndarray, p: np.ndarray, step: int) -> Dict[str, Any]:
    st = max(1, min(step, len(ts) - 1))
    dt = ts[st:] - ts[:-st]
    ok = dt > 0
    v = np.linalg.norm(p[st:] - p[:-st], axis=1)[ok] / dt[ok]
    w = qangle(q[:-st], q[st:])[ok] / dt[ok]
    drot, dpos = dev_deg(q), np.linalg.norm(p - p[0], axis=1)
    return dict(
        rot_p95_deg=float(np.percentile(drot, 95)),
        pos_p95_m=float(np.percentile(dpos, 95)),
        rot_max_deg=float(drot.max()),
        pos_max_m=float(dpos.max()),
        v_p95=float(np.percentile(v, 95)),
        w_p95=float(np.percentile(w, 95)),
    )


def dilate(m: np.ndarray, r: int) -> np.ndarray:
    """Szybka dylatacja logiczna maski binarnej o r pikseli (4-sąsiedztwo)."""
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


def hist_line(vals: np.ndarray, edges: List[float], unit: str) -> str:
    """Formatowanie jednoliniowego histogramu tekstowego."""
    h, _ = np.histogram(vals, bins=edges + [np.inf])
    names = [f"[{a:g},{b:g})" for a, b in zip(edges, edges[1:])] + [f">={edges[-1]:g}"]
    return "  " + "  ".join(f"{nm}: {n}" for nm, n in zip(names, h)) + f"   ({unit}, n = {len(vals)})"
