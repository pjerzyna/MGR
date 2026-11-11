#!/usr/bin/env python3
"""Detekcja IMO v0: kompensacja rotacji -> mapa zdarzeń niewyjaśnionych ruchem własnym -> klasteryzacja.

Użycie:
    python3 imo_detect.py <katalog_sekwencji> [--window-ms 40] [--max-windows 40] [--blur 3] [--tau 0.5]
                          [--support 1.0] [--min-cells 8] [--lag-const-ms X] [--hi-px 2.0] [--out DIR]

Algorytm (dla okna zdarzeń o środku t_c i długości W):
  1. kompensacja rotacji z IMU (rotcomp.py; stałe kalibracyjne, opóźnienie jako parametr) do chwili t_c,
  2. podział okna na połowy (wcześniejsza / późniejsza) i dwa obrazy zdarzeń E, L po kompensacji,
     znormalizowane liczbą zdarzeń i wygładzone (filtr pudełkowy 2x o promieniu --blur),
  3. mapa niespójności  D = |E - L| / (E + L): tło statyczne po kompensacji pokrywa się w obu połowach
     (D ~ 0), obiekt poruszający się niezależnie przesuwa się (D ~ 1),
  4. maska = (D >= --tau) i (E + L >= --support * średnia): odrzuca szum o małej gęstości,
  5. klasteryzacja: składowe spójne (8-sąsiedztwo) na siatce komórek 4x4 px (opcjonalnie po domknięciu --close),
     odrzucenie składowych mniejszych niż --min-cells komórek; zdarzenie jest IMO, gdy leży w zachowanej
     składowej (--fill: także w jej wnętrzu, tj. w dziurach otoczonych obrysem).
Brak parametrów uczonych ani GT w detekcji. Model zakłada czysty obrót: translacja kamery daje w tle
niespójność proporcjonalną do paralaksy (znane ograniczenie v0).

Ocena na poziomie zdarzeń (GT: maski, pozy -> prędkości z imo_labels.py; tylko do oceny):
  etykiety zdarzeń (--gt-mode temporal, domyślnie): każde zdarzenie dostaje etykietę z maski NAJBLIŻSZEJ w czasie
  klatki (<= 12 ms); --gt-mode center używa jednej maski w środku okna (zaniża FPR/czułość dla szybkich obiektów).
  pozytywne = zdarzenia w masce SZYBKICH IMO (v > 0,3 m/s, maska + 6 px), negatywne = poza maską
  wszystkich obiektów ruchomych (v > 0,05 m/s lub omega > 1 rad/s); zdarzenia w maskach obiektów powolnych są pomijane.
  czułość = TP/P, FPR = FP/N, precyzja = TP/(TP+FP), IoU zdarzeń = TP/(P + FP) (suma po oknach z >= --min-pos
  pozytywnymi; FPR po wszystkich oknach). Warianty: brak kompensacji, oracle d_GT, stałe d (mediana opóźnień
  z INNYCH sekwencji w $EVIMO2_OUT), zero. Wyniki osobno dla okien z wyraźną ego-rotacją (p95 > --hi-px).
--translation: dodaje warianty ORACLE (GT w przetwarzaniu - tylko górna granica, nie algorytm): "rotacja IMU + translacja"
(translacja kamery z poz GT i głębia z map GT per zdarzenie, M z IMU), "rotacja GT" i "rotacja + translacja GT"
(pełne 6-DoF z GT). Porównanie z wariantem tylko-rotacja pokazuje, ile fałszywych alarmów wynika z translacji.
--grid: zamiast jednego zestawu parametrów liczy siatkę --grid-blur x --grid-tau (domyślnie {3,4,5,6,8} x {0,2; 0,3; 0,4}; warianty: brak kompensacji, stałe) -
pokazuje kompromis między czułością a odpornością na błędy kompensacji (drobne rozmycie = większa czułość,
ale większe wymagania co do dokładności kompensacji).
Zapisuje <out>/imo_detect.json.
"""
import argparse
import json
import sys
from collections import deque
from pathlib import Path

import numpy as np

import evimo2_common as ec
import imo_labels as il
import object_motion as om
import rotcomp as rc
from event_lag_check import gt_speed_series, spearman, window_features
from gyro_vs_gt import unify_hemisphere
from survey_sequences import pose_series

CELL = 4


# ---------- detektor ----------
def box_blur(img, r):
    k = 2 * r + 1
    p = np.pad(img, r, mode="edge")
    c = np.zeros((p.shape[0] + 1, p.shape[1] + 1))
    c[1:, 1:] = p.cumsum(0).cumsum(1)
    return (c[k:, k:] - c[:-k, k:] - c[k:, :-k] + c[:-k, :-k]) / (k * k)


def residual_map(xs, ys, te, tc, r, eps=1e-12):
    early = te < tc
    E = rc.iwe(xs[early], ys[early])
    L = rc.iwe(xs[~early], ys[~early])
    ne, nl = E.sum(), L.sum()
    if ne < 1 or nl < 1:
        return None
    E, L = E / ne, L / nl
    for _ in range(2):
        E, L = box_blur(E, r), box_blur(L, r)
    S = E + L
    return np.abs(E - L) / (S + eps), S


def label_cc(g):
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


def fill_holes(g):
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


def detect(xs, ys, te, tc, a):
    """Etykieta IMO dla każdego zdarzenia i liczba zachowanych składowych."""
    pred = np.zeros(len(xs), bool)
    rm = residual_map(xs, ys, te, tc, a.blur)
    if rm is None:
        return pred, 0
    D, S = rm
    m = (D >= a.tau) & (S >= a.support * S.mean())
    H, W = m.shape
    g = m.reshape(H // CELL, CELL, W // CELL, CELL).mean(axis=(1, 3)) >= 0.3
    close = getattr(a, "close", 0)
    if close:                                          # domknięcie: dylatacja, potem erozja (kwadrat, close komórek)
        g = ~il.dilate(~il.dilate(g, close), close)
    lab, areas = label_cc(g)
    keep = np.nonzero(areas >= a.min_cells)[0]
    keep = keep[keep > 0]
    if len(keep) == 0:
        return pred, 0
    kept = np.isin(lab, keep)
    if getattr(a, "fill", False):                      # wnętrza obiektów otoczonych wykrytym obrysem
        kept = fill_holes(kept)
    cx, cy = np.floor(xs / CELL).astype(int), np.floor(ys / CELL).astype(int)
    ok = (cx >= 0) & (cx < g.shape[1]) & (cy >= 0) & (cy < g.shape[0])
    pred = np.zeros(len(xs), bool)
    pred[ok] = kept[cy[ok], cx[ok]]
    return pred, len(keep)


# ---------- wizualizacja ----------
def outline(m):
    """Brzeg maski logicznej (piksele maski sąsiadujące z tłem)."""
    return m & ~(~il.dilate(~m, 1))


def render_events(x, y, posm, negm, pred, fm=None, ex=None):
    """Obraz zdarzeń we współrzędnych czujnika: TP zielony, FP czerwony, FN niebieski, TN szary, pominięte pomarańczowy."""
    img = np.ones((480, 640, 3), np.float32)
    amb, tn = ~posm & ~negm, negm & ~pred
    fn, fp, tp = posm & ~pred, negm & pred, posm & pred
    for m, c in ((amb, (1.0, 0.85, 0.5)), (tn, (0.82, 0.82, 0.82)), (fn, (0.1, 0.3, 0.95)), (fp, (0.9, 0.1, 0.1)), (tp, (0.1, 0.65, 0.2))):
        img[y[m], x[m]] = c
    if ex is not None:
        img[outline(ex)] = (0.75, 0.1, 0.75)              # obrys obiektów ruchomych (v > 0,05 m/s), fiolet
    if fm is not None:
        img[outline(fm)] = (0.0, 0.0, 0.0)                # obrys szybkich IMO (v > 0,3 m/s), czerń
    return img


def save_viz(store, a, odir, names):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("Brak matplotlib - pominięto wizualizacje.")
        return
    vdir = odir / "imo_detect_viz"
    vdir.mkdir(parents=True, exist_ok=True)
    pick = sorted(set(np.linspace(0, len(store) - 1, min(a.viz, len(store))).astype(int).tolist()))
    for i in pick:
        w = store[i]
        fig, axes = plt.subplots(1, len(names), figsize=(6.4 * len(names), 4.9))
        for ax, nm in zip(np.atleast_1d(axes), names):
            pred = w["pred"][nm]
            tp, fp = int((pred & w["pos"]).sum()), int((pred & w["neg"]).sum())
            rec = tp / max(int(w["pos"].sum()), 1)
            fpr = fp / max(int(w["neg"].sum()), 1)
            ax.imshow(render_events(w["x"], w["y"], w["pos"], w["neg"], pred, w.get("fm"), w.get("ex")), interpolation="nearest")
            ax.set_title(f"{nm[:34]}\nczułość {rec:.2f}, FPR {fpr:.3f}", fontsize=10)
            ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle(f"{w['seq']} | klatka {w['fid']} | t = {w['ts']:.3f} s | zielony TP, czerwony FP, niebieski FN, "
                     "szary TN, pomarańczowy pominięte (powolne obiekty); obrys: czarny = szybkie IMO, fiolet = wszystkie ruchome (maski GT)", fontsize=8)
        fig.savefig(vdir / f"win_{ec.seq_label(w['seq'])}_{w['fid']:05d}.png", dpi=110, bbox_inches="tight")
        plt.close(fig)
    print(f"Zapisano {len(pick)} wizualizacji do {vdir}/")


def depth_for_events(zd, zdfiles, fts, fid, te, x, y):
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


# ---------- ocena ----------
def regions_for_frame(z, fid, tsf, kin, a):
    """Maski GT: (obiekty ruchome v > slow_v lub omega > w_thr, rozszerzone) i (szybkie v > fast_v, rozszerzone)."""
    lab = z[f"mask_{fid:010d}"] // 1000
    slow, fast = [], []
    for oid in np.unique(lab[lab > 0]):
        k = str(int(oid))
        if k not in kin:
            continue
        v = float(np.interp(tsf, kin[k]["ts"], kin[k]["v"])); w = float(np.interp(tsf, kin[k]["ts"], kin[k]["w"]))
        if v > a.slow_v or w > a.w_thr:
            slow.append(int(oid))
        if v > a.fast_v:
            fast.append(int(oid))
    excl = il.dilate(np.isin(lab, slow), a.dilate) if slow else np.zeros(lab.shape, bool)
    fmask = il.dilate(np.isin(lab, fast), a.dilate) if fast else np.zeros(lab.shape, bool)
    return excl, fmask



class Acc:
    def __init__(self):
        self.tp = self.fp = self.p = self.n = 0
        self.det_fp = self.det_n = 0
        self.comps = []
        self.w_recall = []

    def add(self, tp, fp, p, n, comps, min_pos):
        self.n += n; self.fp += fp
        self.comps.append(comps)
        if p >= min_pos:
            self.tp += tp; self.p += p; self.det_fp += fp; self.det_n += n; self.w_recall.append(tp / p)

    def metrics(self):
        r = lambda x: float(x)
        rec = self.tp / self.p if self.p else float("nan")
        prec = self.tp / (self.tp + self.det_fp) if (self.tp + self.det_fp) else float("nan")
        iou = self.tp / (self.p + self.det_fp) if self.p else float("nan")
        fpr = self.fp / self.n if self.n else float("nan")
        return dict(recall=rec, precision=prec, iou=iou, fpr=fpr, comps_median=float(np.median(self.comps)) if self.comps else float("nan"),
                    windows_det=len(self.w_recall), base_rate=(self.p / self.det_n) if self.det_n else float("nan"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("seq_dir")
    ap.add_argument("--window-ms", type=float, default=40.0)
    ap.add_argument("--max-windows", type=int, default=40)
    ap.add_argument("--blur", type=int, default=3)
    ap.add_argument("--tau", type=float, default=0.5)
    ap.add_argument("--support", type=float, default=1.0)
    ap.add_argument("--min-cells", type=int, default=8)
    ap.add_argument("--min-pos", type=int, default=300, help="min. liczba pozytywnych zdarzeń w oknie do oceny czułości")
    ap.add_argument("--lag-const-ms", type=float, default=None, help="stałe opóźnienie (domyślnie mediana z innych sekwencji)")
    ap.add_argument("--hi-px", type=float, default=2.0)
    ap.add_argument("--slow-v", type=float, default=0.05)
    ap.add_argument("--fast-v", type=float, default=0.3)
    ap.add_argument("--w-thr", type=float, default=1.0)
    ap.add_argument("--dilate", type=int, default=6)
    ap.add_argument("--gt-mode", choices=["temporal", "center"], default="temporal",
                    help="etykiety zdarzeń: z maski NAJBLIŻSZEJ w czasie klatki (temporal) albo z jednej maski w środku okna (center)")
    ap.add_argument("--grid", action="store_true", help="siatka blur x tau dla wariantów: brak kompensacji, oracle, stałe")
    ap.add_argument("--grid-blur", type=int, nargs="+", default=[3, 4, 5, 6, 8], help="promienie rozmycia w --grid")
    ap.add_argument("--grid-tau", type=float, nargs="+", default=[0.2, 0.3, 0.4], help="progi tau w --grid")
    ap.add_argument("--viz", type=int, default=0,
                    help="zapisz tyle wizualizacji okien (zdarzenia: TP/FP/FN/TN, bez kompensacji vs stałe d) do <out>/imo_detect_viz/")
    ap.add_argument("--translation", action="store_true",
                    help="dodaj warianty ORACLE z translacją: ruch kamery z poz GT i głębia z map GT (dataset_depth.npz, mm)")
    ap.add_argument("--fill", action="store_true", help="wypełniaj dziury w wykrytych składowych (wnętrza obiektów otoczone obrysem)")
    ap.add_argument("--close", type=int, default=0, help="domknięcie morfologiczne maski przed klasteryzacją [komórki 4x4 px]")
    ap.add_argument("--min-events", type=int, default=5000)
    ap.add_argument("--cls", default="R*S", choices=om.CLASSES)
    ap.add_argument("--out", default=None)
    ap.add_argument("--out-root", default=None)
    a = ap.parse_args()

    seq = Path(a.seq_dir.rstrip("/"))
    info = ec.load_info(seq)
    if info is None:
        sys.exit("brak dataset_info.npz")
    intr = ec.find_dict_with(info, ["fx", "fy", "cx", "cy"])
    imus = [(p, l) for p, l in ec.find_series(info).items() if ec.is_imu_series(l)]
    imus = [x for x in imus if "left" in x[0]] or imus
    cal_ts, gyro, _ = ec.imu_arrays(imus[0][1])
    cal_w = rc.calibrated_gyro(gyro)
    gj_path = ec.find_result_file(seq.name, "gyro_vs_gt.json", a.out_root)
    d_oracle = float(json.load(open(gj_path))["lag_s"]) if gj_path else None
    if a.lag_const_ms is not None:
        d_const, n_const = a.lag_const_ms * 1e-3, -1
    else:
        d_const, n_const = ec.lag_const_loo(a.out_root, seq.name, rc.LAG_CONST_S)

    cam = pose_series(info["meta"].get("full_trajectory") or [], "cam")
    if cam is None:
        sys.exit("brak poz kamery w full_trajectory (potrzebne do oceny)")
    ts_gt, q_gt, p_gt = cam
    q_gt = unify_hemisphere(q_gt)
    _, tracks = om.load_tracks(info)
    kin = il.kinematics(tracks, a.cls, 0.025)
    t_ev, xy_ev, _ = ec.load_events(seq)
    fts = {int(f["id"]): float(f["ts"]) for f in info["meta"].get("frames", [])}
    W = a.window_ms * 1e-3
    lo = max(cal_ts[0], float(t_ev[0]), ts_gt[0]) + W / 2 + 0.06
    hi_t = min(cal_ts[-1], float(t_ev[-1]), ts_gt[-1]) - W / 2 - 0.06
    mpath = seq / "dataset_mask.npz"
    if not mpath.is_file():
        sys.exit("brak dataset_mask.npz")
    z = np.load(mpath)
    zfiles = set(z.files)
    zd, zdfiles = None, set()
    if a.translation:
        dpath = seq / "dataset_depth.npz"
        if not dpath.is_file():
            sys.exit("--translation wymaga dataset_depth.npz")
        zd = np.load(dpath)
        zdfiles = set(zd.files)
    frames = sorted((fid, fts[fid]) for fid in (int(k.split("_")[-1]) for k in z.files) if fid in fts and lo <= fts[fid] <= hi_t)
    if not frames:
        sys.exit("brak klatek z maską w dopuszczalnym przedziale czasu")
    sel = sorted(set(np.linspace(0, len(frames) - 1, min(a.max_windows, len(frames))).astype(int).tolist()))
    frames = [frames[i] for i in sel]

    modes = [("brak kompensacji", None)]
    if d_oracle is not None:
        modes.append((f"oracle d_GT = {d_oracle * 1e3:+.1f} ms", d_oracle))
    label_const = f"stałe d = {d_const * 1e3:+.1f} ms" + (f" (mediana z {n_const} innych)" if n_const > 0 else "")
    modes += [(label_const, d_const), ("zero d = 0", 0.0)]
    if a.grid:
        psets = [(b, t) for b in a.grid_blur for t in a.grid_tau]
        modes = [m for m in modes if m[0] == "brak kompensacji" or m[0].startswith("stałe")]
    else:
        psets = [(a.blur, a.tau)]
    if a.translation and not a.grid:
        modes += [("rotacja IMU + translacja (oracle)", dict(kind="rt_imu", d=d_const)),
                  ("rotacja GT (oracle)", dict(kind="r_gt")),
                  ("rotacja+translacja GT (oracle)", dict(kind="rt_gt"))]
    acc = {(name, ps): dict(all=Acc(), hi=Acc()) for name, _ in modes for ps in psets}
    n_win = n_hi = n_det = 0
    speed = gt_speed_series(info)
    wrec = []                                    # cechy ego-ruchu i FPR okna (do korelacji)
    vstore = []                                  # okna do wizualizacji (--viz)
    dfracs = []                                  # udział zdarzeń z poprawną głębią GT (--translation)
    viz_names = [m[0] for m in modes if m[0] == "brak kompensacji" or m[0].startswith("stałe")]

    for fid, tsf in frames:
        win = rc.prep_window(tsf, W, t_ev, xy_ev, intr)
        if win["n"] < a.min_events:
            continue
        x = np.clip(win["pix"][:, 0].astype(int), 0, 639)
        y = np.clip(win["pix"][:, 1].astype(int), 0, 479)
        if a.gt_mode == "center":
            excl, fmask = regions_for_frame(z, fid, tsf, kin, a)
            cen_ex, cen_fm = excl, fmask
            negm, posm = ~excl[y, x], fmask[y, x]
        else:                                    # każde zdarzenie dostaje etykietę z najbliższej w czasie klatki z maską
            cand = [f for f in (fid - 1, fid, fid + 1) if f in fts and f"mask_{f:010d}" in zfiles]
            ts_c = np.array([fts[f] for f in cand])
            regs = [regions_for_frame(z, f, fts[f], kin, a) for f in cand]
            ex, fm = np.stack([r[0] for r in regs]), np.stack([r[1] for r in regs])
            cen_ex, cen_fm = regs[cand.index(fid)] if fid in cand else regs[0]
            dd = np.abs(win["te"][:, None] - ts_c[None, :])
            fi, valid = dd.argmin(1), dd.min(1) <= 0.012
            negm, posm = (~ex[fi, y, x]) & valid, fm[fi, y, x] & valid
        n_pos, n_neg = int(posm.sum()), int(negm.sum())
        depth_z, dfrac = depth_for_events(zd, zdfiles, fts, fid, win["te"], x, y) if a.translation else (None, None)
        if a.translation:
            dfracs.append(dfrac)
        x0, y0 = rc.warp(win, cal_ts, cal_w, None)
        xg, yg = rc.gt_warp(win, ts_gt, q_gt)
        is_hi = np.percentile(np.hypot(x0 - xg, y0 - yg), 95) >= a.hi_px
        n_win += 1; n_hi += int(is_hi); n_det += int(n_pos >= a.min_pos)
        rec = dict(raw_p95=float(np.percentile(np.hypot(x0 - xg, y0 - yg), 95)), n_neg=n_neg,
                   v_cam=window_features(tsf, W, speed, cal_ts, cal_w, float(intr["fx"]), 1.0).get("v_mps", float("nan")), fpr={})
        wrec.append(rec)
        vrec = dict(seq=seq.name, fid=fid, ts=tsf, x=x.astype(np.int16), y=y.astype(np.int16), pos=posm, neg=negm, pred={}, ex=cen_ex, fm=cen_fm) \
            if (a.viz and len(psets) == 1 and n_pos >= a.min_pos) else None
        for name, d in modes:
            if isinstance(d, dict):                                   # warianty oracle z GT
                if d["kind"] != "r_gt" and depth_z is None:
                    continue
                if d["kind"] == "r_gt":
                    xs, ys = xg, yg
                elif d["kind"] == "rt_imu":
                    xs, ys = rc.warp_translation(win, rc.rot_matrices(win, cal_ts, cal_w, d["d"]), ts_gt, q_gt, p_gt, depth_z)
                else:
                    xs, ys = rc.warp_translation(win, rc.gt_rot_matrices(win, ts_gt, q_gt), ts_gt, q_gt, p_gt, depth_z)
            else:
                xs, ys = (x0, y0) if d is None else rc.warp(win, cal_ts, cal_w, d)
            for ps in psets:
                pa = argparse.Namespace(blur=ps[0], tau=ps[1], support=a.support, min_cells=a.min_cells, fill=a.fill, close=a.close)
                pred, comps = detect(xs, ys, win["te"], tsf, pa)
                tp, fp = int((pred & posm).sum()), int((pred & negm).sum())
                if len(psets) == 1 and n_neg >= 2000:
                    rec["fpr"][name] = fp / n_neg
                if vrec is not None and name in viz_names:
                    vrec["pred"][name] = pred
                acc[(name, ps)]["all"].add(tp, fp, n_pos, n_neg, comps, a.min_pos)
                if is_hi:
                    acc[(name, ps)]["hi"].add(tp, fp, n_pos, n_neg, comps, a.min_pos)
        if vrec is not None:
            vstore.append(vrec)

    if n_win == 0:
        sys.exit("żadne okno nie miało dość zdarzeń")
    print(f"Sekwencja: {seq.name} | okien: {n_win} x {a.window_ms:.0f} ms | z szybkimi IMO (>= {a.min_pos} zdarzeń): {n_det} | "
          f"z ego-rotacją >= {a.hi_px:g} px: {n_hi}")
    out = dict(sequence=seq.name, windows=n_win, windows_with_fast_imo=n_det, windows_hi_rotation=n_hi,
               params=dict(window_ms=a.window_ms, support=a.support, min_cells=a.min_cells), modes={})
    strata = (("all", "WSZYSTKIE OKNA"), ("hi", f"OKNA Z WYRAŹNĄ EGO-ROTACJĄ (p95 przemieszczenia >= {a.hi_px:g} px)"))
    if not a.grid:
        ps = psets[0]
        print(f"Parametry: rozmycie {ps[0]} (2x), próg tau = {ps[1]}, wsparcie {a.support}, min. {a.min_cells} komórek {CELL}x{CELL} px. "
              "Czułość/precyzja/IoU: okna z szybkimi IMO; FPR: wszystkie okna.\n")
        out["params"].update(blur=ps[0], tau=ps[1])
        for stratum, title in strata:
            print(title)
            print(f"{'wariant':40s} {'czułość':>8s} {'precyzja':>9s} {'IoU zdarz.':>11s} {'FPR tło':>8s} {'skł./okno':>10s} {'okien det.':>11s} {'P/N':>6s}")
            for name, _ in modes:
                m = acc[(name, ps)][stratum].metrics()
                out["modes"].setdefault(name, {})[stratum] = m
                print(f"{name:40s} {m['recall']:8.2f} {m['precision']:9.2f} {m['iou']:11.2f} {m['fpr']:8.3f} {m['comps_median']:10.1f} {m['windows_det']:11d} {m['base_rate']:6.3f}")
            print()
        print("P/N = stosunek liczby zdarzeń pozytywnych (szybkie IMO) do negatywnych w oknach detekcji: przy P/N rzędu 0,05 nawet FPR 5% daje precyzję poniżej 0,5.")
        for name, _ in modes:
            pts = [(r["fpr"][name], r["raw_p95"], r["v_cam"]) for r in wrec if name in r["fpr"] and np.isfinite(r["v_cam"])]
            if len(pts) >= 8:
                f_, rr, vv = ([p[i] for p in pts] for i in range(3))
                out["modes"][name]["fpr_corr"] = dict(n=len(pts), vs_raw_rotation=spearman(f_, rr), vs_cam_speed=spearman(f_, vv))
                print(f"FPR okna vs ego-ruch [{name[:24]}]: Spearman z ego-rotacją (p95 px) = {spearman(f_, rr):+.2f}, "
                      f"z prędkością kamery |v| = {spearman(f_, vv):+.2f} (n = {len(pts)}, próg istotności ok. ±{2 / np.sqrt(len(pts)):.2f})")
    else:
        print("Siatka parametrów (blur = promień rozmycia, tau = próg niespójności); czułość/IoU: okna z szybkimi IMO, FPR: wszystkie okna.\n")
        for stratum, title in strata:
            print(title)
            head = f"{'blur':>4s} {'tau':>4s} | " + " | ".join(f"{n[:26]:^26s}" for n, _ in modes)
            print(head)
            print(f"{'':>4s} {'':>4s} | " + " | ".join(f"{'czułość  FPR   IoU':^26s}" for _ in modes))
            for ps in psets:
                cells = []
                for name, _ in modes:
                    m = acc[(name, ps)][stratum].metrics()
                    out["modes"].setdefault(name, {}).setdefault(f"blur{ps[0]}_tau{ps[1]}", {})[stratum] = m
                    cells.append(f"{m['recall']:7.2f} {m['fpr']:6.3f} {m['iou']:5.2f}".center(26))
                print(f"{ps[0]:4d} {ps[1]:4.1f} | " + " | ".join(cells))
            print()
    print("IoU zdarzeń = TP/(P + FP) w oknach z szybkimi IMO. Wartości nan: brak takich okien w danej grupie.")
    if a.translation and dfracs:
        out["depth_valid_fraction_median"] = float(np.median(dfracs))
        print(f"Głębia GT dostępna dla {100 * np.median(dfracs):.0f}% zdarzeń (mediana okien; min {100 * np.min(dfracs):.0f}%); "
              "pozostałe dostają medianę głębi okna. Warianty oracle używają GT w przetwarzaniu - to górna granica, nie algorytm.")
    odir = Path(a.out) if a.out else ec.result_dir(seq.name, a.out_root)
    odir.mkdir(parents=True, exist_ok=True)
    (odir / "imo_detect.json").write_text(json.dumps(out, indent=1))
    if a.viz and vstore:
        save_viz(vstore, a, odir, viz_names)
    print(f"Zapisano {odir / 'imo_detect.json'}")


if __name__ == "__main__":
    main()