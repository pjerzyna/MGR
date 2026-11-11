"""Detektor niezależnie poruszających się obiektów (IMO v0): mapa niespójności i klasteryzacja."""
from collections import deque
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np

from .kinematics import dilate
from .rotcomp import iwe

CELL = 4


def box_blur(img: np.ndarray, r: int) -> np.ndarray:
    """Szybki filtr pudełkowy 2D z tablicą sum skumulowanych."""
    k = 2 * r + 1
    p = np.pad(img, r, mode="edge")
    c = np.zeros((p.shape[0] + 1, p.shape[1] + 1), dtype=np.float64)
    c[1:, 1:] = p.cumsum(0).cumsum(1)
    return (c[k:, k:] - c[:-k, k:] - c[k:, :-k] + c[:-k, :-k]) / (k * k)


def residual_map(xs: np.ndarray, ys: np.ndarray, te: np.ndarray, tc: float, r: int, eps: float = 1e-12) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Generuje mapę niespójności D = |E - L| / (E + L) pomiędzy połówkami okna zdarzeń."""
    early = te < tc
    E = iwe(xs[early], ys[early])
    L = iwe(xs[~early], ys[~early])
    ne, nl = E.sum(), L.sum()
    if ne < 1 or nl < 1:
        return None
    E, L = E / ne, L / nl
    for _ in range(2):
        E, L = box_blur(E, r), box_blur(L, r)
    S = E + L
    return np.abs(E - L) / (S + eps), S


def label_cc(g: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Składowe spójne (8-sąsiedztwo) maski logicznej: (etykiety, pola; indeks 0 = tło)."""
    H, W = g.shape
    lab = np.zeros((H, W), np.int32)
    areas = [0]
    n = 0
    for y0, x0 in zip(*np.nonzero(g)):
        if lab[y0, x0]:
            continue
        n += 1
        lab[y0, x0] = n
        q, cnt = deque([(y0, x0)]), 0
        while q:
            y, x = q.popleft()
            cnt += 1
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    yy, xx = y + dy, x + dx
                    if 0 <= yy < H and 0 <= xx < W and g[yy, xx] and not lab[yy, xx]:
                        lab[yy, xx] = n
                        q.append((yy, xx))
        areas.append(cnt)
    return lab, np.array(areas)


def fill_holes(g: np.ndarray) -> np.ndarray:
    """Wypełnia dziury maski logicznej (komórki tła niedostępne z brzegu siatki, 4-sąsiedztwo)."""
    H, W = g.shape
    reach = np.zeros_like(g, bool)
    q = deque()
    for y in range(H):
        for x in (0, W - 1):
            if not g[y, x] and not reach[y, x]:
                reach[y, x] = True; q.append((y, x))
    for x in range(W):
        for y in (0, H - 1):
            if not g[y, x] and not reach[y, x]:
                reach[y, x] = True; q.append((y, x))
    while q:
        y, x = q.popleft()
        for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            yy, xx = y + dy, x + dx
            if 0 <= yy < H and 0 <= xx < W and not g[yy, xx] and not reach[yy, xx]:
                reach[yy, xx] = True; q.append((yy, xx))
    return g | ~reach


def detect_imo(
    xs: np.ndarray,
    ys: np.ndarray,
    te: np.ndarray,
    tc: float,
    blur: int = 3,
    tau: float = 0.5,
    support: float = 1.0,
    min_cells: int = 8,
    fill: bool = False,
    close: int = 0,
) -> Tuple[np.ndarray, int]:
    """Detekcja IMO: zwraca maskę logiczną zdarzeń (N,) zaklasyfikowanych jako IMO oraz liczbę klastrów."""
    pred = np.zeros(len(xs), bool)
    rm = residual_map(xs, ys, te, tc, blur)
    if rm is None:
        return pred, 0
    D, S = rm
    m = (D >= tau) & (S >= support * S.mean())
    H, W = m.shape
    g = m.reshape(H // CELL, CELL, W // CELL, CELL).mean(axis=(1, 3)) >= 0.3
    if close:
        g = ~dilate(~dilate(g, close), close)
    lab, areas = label_cc(g)
    keep = np.nonzero(areas >= min_cells)[0]
    keep = keep[keep > 0]
    if len(keep) == 0:
        return pred, 0
    kept = np.isin(lab, keep)
    if fill:
        kept = fill_holes(kept)
    cx, cy = np.floor(xs / CELL).astype(int), np.floor(ys / CELL).astype(int)
    ok = (cx >= 0) & (cx < g.shape[1]) & (cy >= 0) & (cy < g.shape[0])
    pred[ok] = kept[cy[ok], cx[ok]]
    return pred, len(keep)


def outline(m: np.ndarray) -> np.ndarray:
    """Brzeg maski logicznej (piksele maski sąsiadujące z tłem)."""
    return m & ~(~dilate(~m, 1))


def regions_for_frame(
    z: Any,
    fid: int,
    tsf: float,
    kin: Dict[str, Any],
    slow_v: float = 0.05,
    fast_v: float = 0.3,
    w_thr: float = 1.0,
    dilate_px: int = 6,
) -> Tuple[np.ndarray, np.ndarray]:
    """Maski GT dla klatki: (wszystkie obiekty ruchome v > slow_v lub w > w_thr) oraz (szybkie v > fast_v)."""
    lab = z[f"mask_{fid:010d}"] // 1000
    slow, fast = [], []
    for oid in np.unique(lab[lab > 0]):
        k = str(int(oid))
        if k not in kin:
            continue
        v = float(np.interp(tsf, kin[k]["ts"], kin[k]["v"]))
        w = float(np.interp(tsf, kin[k]["ts"], kin[k]["w"]))
        if v > slow_v or w > w_thr:
            slow.append(int(oid))
        if v > fast_v:
            fast.append(int(oid))
    excl = dilate(np.isin(lab, slow), dilate_px) if slow else np.zeros(lab.shape, bool)
    fmask = dilate(np.isin(lab, fast), dilate_px) if fast else np.zeros(lab.shape, bool)
    return excl, fmask


def render_events(
    x: np.ndarray,
    y: np.ndarray,
    posm: np.ndarray,
    negm: np.ndarray,
    pred: np.ndarray,
    fm: Optional[np.ndarray] = None,
    ex: Optional[np.ndarray] = None,
    H: int = 480,
    W: int = 640,
) -> np.ndarray:
    """Wizualizacja klasyfikacji zdarzeń: TP zielony, FP czerwony, FN niebieski, TN szary, pominięte pomarańczowy."""
    img = np.ones((H, W, 3), np.float32)
    amb, tn = ~posm & ~negm, negm & ~pred
    fn, fp, tp = posm & ~pred, negm & pred, posm & pred
    for m, c in (
        (amb, (1.0, 0.85, 0.5)),
        (tn, (0.82, 0.82, 0.82)),
        (fn, (0.1, 0.3, 0.95)),
        (fp, (0.9, 0.1, 0.1)),
        (tp, (0.1, 0.65, 0.2)),
    ):
        img[y[m], x[m]] = c
    if ex is not None:
        img[outline(ex)] = (0.75, 0.1, 0.75)              # obrys ruchomych v > 0.05, fiolet
    if fm is not None:
        img[outline(fm)] = (0.0, 0.0, 0.0)                # obrys szybkich IMO v > 0.3, czerń
    return img


def depth_for_events(
    zd: Any,
    zdfiles: Set[str],
    fts: Dict[int, float],
    fid: int,
    te: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
) -> Tuple[Optional[np.ndarray], float]:
    """Głębia [m] w pikselu każdego zdarzenia z najbliższej w czasie klatki z mapą głębi (mm w pliku);
    brakujące piksele (0) zastępowane medianą głębi okna. Zwraca (głębia albo None, udział zdarzeń z poprawną głębią)."""
    cand = [f for f in (fid - 1, fid, fid + 1) if f in fts and f"depth_{f:010d}" in zdfiles]
    if not cand:
        return None, 0.0
    tsc = np.array([fts[f] for f in cand])
    maps = np.stack([zd[f"depth_{f:010d}"] for f in cand]).astype(np.float32) / 1000.0
    fi = np.abs(te[:, None] - tsc[None, :]).argmin(1)
    z = maps[fi, y, x].astype(np.float64)
    valid = z > 0.05
    if not valid.any():
        return None, 0.0
    return np.where(valid, z, float(np.median(z[valid]))), float(valid.mean())

