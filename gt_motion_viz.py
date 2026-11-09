#!/usr/bin/env python3
"""Wizualizacja Ground Truth ruchu obiektów i etykiet zdarzeń (czyste TP z GT).

Użycie:
    python3 gt_motion_viz.py <katalog_sekwencji> [--frames FID ...] [--stride 10]
                            [--window-ms 40] [--slow-v 0.05] [--fast-v 0.3]
                            [--w-thr 1.0] [--dilate 6] [--out DIR]

Tworzy dwupanelowe wizualizacje czystego Ground Truth:
  * Panel lewy: Zdarzenia (polaryzacja ON/OFF) + obrysy obiektów z GT wraz z prędkościami:
      - Szary: obiekty statyczne (v <= --slow-v i w <= --w-thr)
      - Fiolet/Pomarańcz: obiekty wolne (--slow-v < v <= --fast-v)
      - Zielony: szybkie IMO (v > --fast-v) - niezależnie poruszające się obiekty
  * Panel prawy: Czyste etykiety zdarzeń według Ground Truth:
      - Zielony: Prawdziwe Zdarzenia IMO (True Positives z GT, v > --fast-v)
      - Pomarańczowy: Zdarzenia obiektów wolnych (strefa niejednoznaczna)
      - Morski/Szaroniebieski: Zdarzenia obiektów statycznych
      - Jasnoszary: Zdarzenia tła
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import evimo2_common as ec
import imo_labels as il
import object_motion as om


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("seq_dir", help="Ścieżka do sekwencji EVIMO2")
    ap.add_argument("--frames", type=int, nargs="+", default=None, help="Konkretne ID klatek do wizualizacji")
    ap.add_argument("--stride", type=int, default=15, help="Krok wyboru klatek, jeśli nie podano --frames")
    ap.add_argument("--max-viz", type=int, default=8, help="Maksymalna liczba zapisanych klatek")
    ap.add_argument("--window-ms", type=float, default=40.0, help="Szerokość okna zdarzeń [ms]")
    ap.add_argument("--slow-v", type=float, default=0.05, help="Próg prędkości ruchu [m/s]")
    ap.add_argument("--fast-v", type=float, default=0.3, help="Próg szybkiego IMO [m/s]")
    ap.add_argument("--w-thr", type=float, default=1.0, help="Próg prędkości kątowej [rad/s]")
    ap.add_argument("--dilate", type=int, default=6, help="Dylatacja masek obiektów [px]")
    ap.add_argument("--cls", default="R*S", choices=om.CLASSES, help="Klasa konwencji poz obiektów")
    ap.add_argument("--out", default=None, help="Katalog wyjściowy (domyślnie <out>/gt_motion_viz)")
    ap.add_argument("--out-root", default=None)
    return ap.parse_args()


def get_outline(mask):
    """Zwraca kontur maski binarnej (1 px)."""
    return mask & ~(~il.dilate(~mask, 1))


def main():
    args = parse_args()
    seq = Path(args.seq_dir.rstrip("/"))
    info = ec.load_info(seq)
    if info is None:
        sys.exit(f"Błąd: brak dataset_info.npz w {seq}")

    intr = ec.find_dict_with(info, ["fx", "fy", "cx", "cy", "res_x", "res_y"])
    H = int(intr.get("res_y", 480))
    W = int(intr.get("res_x", 640))

    _, tracks = om.load_tracks(info)
    kin = il.kinematics(tracks, args.cls, 0.025)
    t_ev, xy_ev, p_ev = ec.load_events(seq)
    fts = {int(f["id"]): float(f["ts"]) for f in info["meta"].get("frames", [])}

    mpath = seq / "dataset_mask.npz"
    if not mpath.is_file():
        sys.exit(f"Błąd: brak dataset_mask.npz w {seq}")
    z = np.load(mpath)

    w_s = args.window_ms * 1e-3
    t_min = float(t_ev[0]) + w_s / 2
    t_max = float(t_ev[-1]) - w_s / 2

    avail_frames = sorted(
        (int(k.split("_")[-1]), fts[int(k.split("_")[-1])])
        for k in z.files
        if int(k.split("_")[-1]) in fts and t_min <= fts[int(k.split("_")[-1])] <= t_max
    )
    if not avail_frames:
        sys.exit("Brak klatek z maską w dopuszczalnym zakresie czasowym zdarzeń")

    if args.frames:
        target_fids = set(args.frames)
        frames_to_render = [f for f in avail_frames if f[0] in target_fids]
        if not frames_to_render:
            sys.exit(f"Żadna z podanych klatek {args.frames} nie jest dostępna w poprawnym oknie czasu")
    else:
        indices = np.linspace(0, len(avail_frames) - 1, min(args.max_viz, len(avail_frames))).astype(int)
        frames_to_render = [avail_frames[i] for i in sorted(set(indices))]

    out_dir = Path(args.out) if args.out else ec.result_dir(seq.name, args.out_root) / "gt_motion_viz"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Sekwencja: {seq.name}")
    print(f"Klatek do wyrenderowania: {len(frames_to_render)} (okno {args.window_ms:.0f} ms)")
    print(f"Katalog zapisu: {out_dir}/\n")

    for fid, tsf in frames_to_render:
        lab = z[f"mask_{fid:010d}"] // 1000
        obj_ids = [int(o) for o in np.unique(lab[lab > 0])]

        obj_stats = {}
        for oid in obj_ids:
            k = str(oid)
            if k in kin:
                v = float(np.interp(tsf, kin[k]["ts"], kin[k]["v"]))
                w = float(np.interp(tsf, kin[k]["ts"], kin[k]["w"]))
            else:
                v, w = 0.0, 0.0

            if v > args.fast_v:
                status = "IMO"
            elif v > args.slow_v or w > args.w_thr:
                status = "SLOW"
            else:
                status = "STATIC"

            obj_stats[oid] = dict(v=v, w=w, status=status)

        # Maski dla kategorii
        fast_ids = [oid for oid, s in obj_stats.items() if s["status"] == "IMO"]
        slow_ids = [oid for oid, s in obj_stats.items() if s["status"] == "SLOW"]
        static_ids = [oid for oid, s in obj_stats.items() if s["status"] == "STATIC"]

        fast_mask = np.isin(lab, fast_ids) if fast_ids else np.zeros((H, W), bool)
        slow_mask = np.isin(lab, slow_ids) if slow_ids else np.zeros((H, W), bool)
        static_mask = np.isin(lab, static_ids) if static_ids else np.zeros((H, W), bool)

        # Dylatacja dla stref zdarzeniowych (zgodnie z metodologią pracy)
        fast_dil = il.dilate(fast_mask, args.dilate) if fast_ids else np.zeros((H, W), bool)
        slow_dil = il.dilate(slow_mask, args.dilate) if slow_ids else np.zeros((H, W), bool)
        static_dil = il.dilate(static_mask, args.dilate) if static_ids else np.zeros((H, W), bool)

        # Pobieranie zdarzeń w oknie
        i0, i1 = np.searchsorted(t_ev, [tsf - w_s / 2, tsf + w_s / 2])
        xe = np.clip(xy_ev[i0:i1, 0].astype(int), 0, W - 1)
        ye = np.clip(xy_ev[i0:i1, 1].astype(int), 0, H - 1)
        pe = p_ev[i0:i1]
        n_ev = len(xe)

        # Klasyfikacja zdarzeń przez Ground Truth
        is_fast = fast_dil[ye, xe]
        is_slow = slow_dil[ye, xe] & (~is_fast)
        is_static = static_dil[ye, xe] & (~is_fast) & (~is_slow)
        is_bg = (~is_fast) & (~is_slow) & (~is_static)

        n_fast = int(is_fast.sum())
        n_slow = int(is_slow.sum())
        n_static = int(is_static.sum())
        n_bg = int(is_bg.sum())

        # Render panelu 1: Raw events + Contours + Velocity Labels
        img_left = np.ones((H, W, 3), dtype=np.float32)
        # Zdarzenia tła i obiektów z polaryzacją (OFF: niebieski, ON: czerwony)
        off_mask = pe == 0
        on_mask = pe != 0
        img_left[ye[off_mask], xe[off_mask]] = (0.35, 0.45, 0.85)  # OFF niebieski
        img_left[ye[on_mask], xe[on_mask]] = (0.85, 0.35, 0.35)    # ON czerwony

        # Obrysy obiektów na panelu 1
        for oid, s in obj_stats.items():
            omask = (lab == oid)
            cont = get_outline(omask)
            if s["status"] == "IMO":
                color = (0.05, 0.85, 0.15)  # Jaskrawy zielony
            elif s["status"] == "SLOW":
                color = (0.95, 0.55, 0.10)  # Pomarańczowy
            else:
                color = (0.30, 0.30, 0.30)  # Ciemnoszary
            img_left[cont] = color

        # Render panelu 2: Czysty podział Ground Truth zdarzeń
        img_right = np.ones((H, W, 3), dtype=np.float32)
        # Tło
        img_right[ye[is_bg], xe[is_bg]] = (0.78, 0.78, 0.78)
        # Statyczne obiekty
        img_right[ye[is_static], xe[is_static]] = (0.30, 0.50, 0.75)
        # Wolne obiekty
        img_right[ye[is_slow], xe[is_slow]] = (0.95, 0.60, 0.15)
        # CZYSTE TP: Szybkie IMO
        img_right[ye[is_fast], xe[is_fast]] = (0.05, 0.80, 0.15)

        # Dodanie delikatnych obrysów GT na panelu prawym
        if fast_mask.any():
            img_right[get_outline(fast_mask)] = (0.0, 0.0, 0.0)
        if slow_mask.any():
            img_right[get_outline(slow_mask)] = (0.6, 0.1, 0.6)

        # Tworzenie wykresu
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.2))

        ax1.imshow(img_left, interpolation="nearest")
        ax1.set_title(
            f"Strumień zdarzeń ({n_ev:,} zdarzeń w {args.window_ms:.0f} ms)\n"
            f"Niebieski=OFF, Czerwony=ON | Kontury: Zielony=IMO, Pomarańcz=Wolne, Szary=Statyczne",
            fontsize=10,
        )
        ax1.set_xticks([])
        ax1.set_yticks([])

        # Podpisy obiektów na panelu lewym
        for oid, s in obj_stats.items():
            ys, xs = np.where(lab == oid)
            if len(xs) > 0:
                xc, yc = int(xs.mean()), int(ys.mean())
                tag = f"ID{oid}: {s['v']:.2f}m/s"
                c_txt = "darkgreen" if s["status"] == "IMO" else ("darkorange" if s["status"] == "SLOW" else "black")
                ax1.text(
                    xc,
                    yc,
                    tag,
                    color=c_txt,
                    fontsize=8,
                    fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor=c_txt),
                )

        ax2.imshow(img_right, interpolation="nearest")
        ax2.set_title(
            f"Ground Truth Etykiet Zdarzeń (Czysty TP z GT)\n"
            f"Zielony=TP IMO ({n_fast:,}, {100*n_fast/max(1,n_ev):.1f}%) | Pomarańcz=Wolne ({n_slow:,}) | "
            f"Niebieski=Stat ({n_static:,}) | Szary=Tło ({n_bg:,})",
            fontsize=10,
        )
        ax2.set_xticks([])
        ax2.set_yticks([])

        fig.suptitle(
            f"{seq.name} | Klatka {fid:05d} | t = {tsf:.3f} s (okno ±{w_s*500:.0f} ms)\n"
            f"Obiektów w klatce: {len(obj_stats)} (Szybkie IMO: {len(fast_ids)}, Wolne: {len(slow_ids)}, Statyczne: {len(static_ids)})",
            fontsize=11,
            fontweight="bold",
        )

        out_file = out_dir / f"gt_viz_{fid:05d}.png"
        fig.tight_layout()
        fig.savefig(out_file, dpi=120)
        plt.close(fig)
        print(f"  [OK] Zapisano klatkę {fid:05d} -> {out_file.name} (TP IMO: {n_fast:,} zdarzeń)")

    print(f"\nUkończono! Wszystkie pliki zapisane w: {out_dir}")


if __name__ == "__main__":
    main()
