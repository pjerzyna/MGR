#!/usr/bin/env python3
"""Generate dedicated ROC curves and comparative analysis plots specifically for r=3, tau=0.3."""
import json
import csv
import shutil
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

EVAL_DIR = Path("results/v1_rotcomp/benchmark")
ARTIFACT_DIR = Path("/home/pawsooon/.gemini/antigravity/brain/5814172c-e6cf-4394-ae23-449e5e16c752")

SEQUENCES = [
    ("scene6_dyn_train_03_000000", "Spoczynek (czyste IMO)", "Kategoria I"),
    ("scene14_dyn_test_05_000000", "Czysta rotacja (wzorzec)", "Kategoria II"),
    ("scene14_dyn_test_04_000000", "Czysta rotacja (mały ruch)", "Kategoria II"),
    ("scene13_dyn_test_02_000000", "Szybka rotacja (stress-test FPR)", "Kategoria II"),
    ("scene13_dyn_test_01_000000", "Ruch złożony (silna paralaksa)", "Kategoria III"),
    ("scene15_dyn_test_02_000000", "Ruch złożony (duży udział IMO)", "Kategoria III"),
    ("scene15_dyn_test_05_000000", "Ruch złożony (wysokie IMO)", "Kategoria III"),
]

def load_data():
    grid_data = {}
    for seq_name, cat_desc, cat_group in SEQUENCES:
        json_file = EVAL_DIR / f"grid_{seq_name}" / "imo_detect.json"
        if not json_file.exists():
            alt = list(EVAL_DIR.glob(f"*{seq_name}*/*.json"))
            if alt:
                json_file = alt[0]
            else:
                continue
        with open(json_file) as f:
            grid_data[seq_name] = json.load(f)
    return grid_data

def generate_roc_tau03_plot(grid_data):
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, axes = plt.subplots(1, 3, figsize=(19, 6.2), dpi=160)
    
    categories = [
        ("Kategoria I: Spoczynek", ["scene6_dyn_train_03_000000"], axes[0]),
        ("Kategoria II: Dominujący obrót", ["scene14_dyn_test_05_000000", "scene14_dyn_test_04_000000", "scene13_dyn_test_02_000000"], axes[1]),
        ("Kategoria III: Ruch złożony / paralaksa", ["scene13_dyn_test_01_000000", "scene15_dyn_test_02_000000", "scene15_dyn_test_05_000000"], axes[2]),
    ]
    
    colors = ['#1f77b4', '#2ca02c', '#d62728']
    taus = [0.2, 0.3, 0.4, 0.5]
    
    for cat_title, seq_list, ax in categories:
        ax.set_title(cat_title, fontsize=13, fontweight='bold', pad=12)
        ax.set_xlabel("Odsetek błędów tła (FPR) [-]", fontsize=11, labelpad=8)
        ax.set_ylabel("Czułość detekcji (Recall) [-]", fontsize=11, labelpad=8)
        ax.grid(True, linestyle="--", alpha=0.6)
        ax.set_ylim(-0.02, 0.85)
        ax.set_xlim(-0.02, 0.55)
        
        for idx, s in enumerate(seq_list):
            if s not in grid_data:
                continue
            modes = grid_data[s].get("modes", {})
            none_k = [k for k in modes.keys() if "brak" in k.lower()][0]
            rot_k  = [k for k in modes.keys() if "stałe" in k.lower() or "oracle" in k.lower()][0]
            
            c = colors[idx % len(colors)]
            short_name = s.replace("_000000", "").replace("_dyn_test", "").replace("_dyn_train", "")
            
            rec_none, fpr_none = [], []
            rec_rot, fpr_rot = [], []
            
            pt_tau03_none = None
            pt_tau03_rot = None
            
            for tau in taus:
                pk = f"blur3_tau{tau}"
                if pk in modes[none_k]:
                    rn = modes[none_k][pk]["all"]["recall"]
                    fn = modes[none_k][pk]["all"]["fpr"]
                    rec_none.append(rn)
                    fpr_none.append(fn)
                    if tau == 0.3:
                        pt_tau03_none = (fn, rn)
                if pk in modes[rot_k]:
                    rr = modes[rot_k][pk]["all"]["recall"]
                    fr = modes[rot_k][pk]["all"]["fpr"]
                    rec_rot.append(rr)
                    fpr_rot.append(fr)
                    if tau == 0.3:
                        pt_tau03_rot = (fr, rr)
            
            # Plot background full ROC curves
            if fpr_none and rec_none:
                ax.plot(fpr_none, rec_none, '--o', color=c, alpha=0.35, lw=1.2, ms=4, label=f"{short_name} (Brak)")
            if fpr_rot and rec_rot:
                ax.plot(fpr_rot, rec_rot, '-s', color=c, alpha=0.45, lw=1.5, ms=5, label=f"{short_name} (SO3)")
                
            # Highlight target operating point tau=0.3
            if pt_tau03_none and pt_tau03_rot:
                # None marker
                ax.plot(pt_tau03_none[0], pt_tau03_none[1], 'o', color=c, ms=8, mec='black', mew=1.2, alpha=0.85)
                # Rot marker (star / diamond)
                ax.plot(pt_tau03_rot[0], pt_tau03_rot[1], '*', color=c, ms=14, mec='black', mew=1.0)
                
                # Shift arrow from Brak to SO(3)
                dx = pt_tau03_rot[0] - pt_tau03_none[0]
                dy = pt_tau03_rot[1] - pt_tau03_none[1]
                if abs(dx) > 0.015:
                    ax.annotate("", xy=pt_tau03_rot, xytext=pt_tau03_none,
                                arrowprops=dict(arrowstyle="->", color=c, lw=2.0, alpha=0.85, shrinkA=5, shrinkB=6))
                    
                    # Annotate reduction
                    mid_x = (pt_tau03_none[0] + pt_tau03_rot[0]) / 2.0
                    mid_y = (pt_tau03_none[1] + pt_tau03_rot[1]) / 2.0
                    red = pt_tau03_none[0] / pt_tau03_rot[0] if pt_tau03_rot[0] > 0 else 1.0
                    if red > 1.3:
                        ax.text(mid_x, mid_y + 0.035, f"{red:.1f}× mniej\nbłędów", fontsize=8.5,
                                ha='center', va='bottom', color=c, fontweight='bold',
                                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=c, alpha=0.7, lw=0.8))

        # Add a custom legend element explaining marker styles
        ax.legend(fontsize=8.5, loc="lower right", framealpha=0.9)
    
    fig.suptitle("Analiza Punktu Zbalansowanego: Krzywe ROC i przesunięcie punktów pracy dla r=3, τ=0.3\n"
                 "(Gwiazdka (*): Kompensacja SO(3), Kółko (o): Brak kompensacji, Strzałka: redukcja fałszywych alarmów tła)",
                 fontsize=14, fontweight='bold', y=0.99)
    fig.tight_layout()
    out_path = EVAL_DIR / "roc_curves_tau03.png"
    fig.savefig(out_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    print(f"Saved {out_path}")
    shutil.copy(out_path, ARTIFACT_DIR / "roc_curves_tau03.png")

def generate_metrics_summary_tau03(grid_data):
    # Bar chart comparing Brak vs SO3 specifically at r=3, tau=0.3
    qual_csv = EVAL_DIR / "quality_metrics_7seqs.csv"
    if not qual_csv.exists():
        return
        
    rows = []
    with open(qual_csv) as f:
        r = csv.DictReader(f)
        for row in r:
            if "tau=0.3" in row["point"]:
                rows.append(row)
                
    if not rows:
        return
        
    seq_labels = [r["sequence"].replace("_000000", "").replace("_dyn_test", "").replace("_dyn_train", "") for r in rows]
    none_fpr = [float(r["none_fpr"])*100 for r in rows]
    rot_fpr = [float(r["rot_fpr"])*100 for r in rows]
    none_prec = [float(r["none_precision"])*100 for r in rows]
    rot_prec = [float(r["rot_precision"])*100 for r in rows]
    none_rec = [float(r["none_recall"])*100 for r in rows]
    rot_rec = [float(r["rot_recall"])*100 for r in rows]
    reductions = [float(r["fpr_reduction_x"]) for r in rows]
    deltas = [float(r["delta_precision_pp"]) for r in rows]
    
    x = np.arange(len(seq_labels))
    width = 0.36
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17, 5.8), dpi=160)
    
    # 1. Błędy tła (FPR)
    b1 = ax1.bar(x - width/2, none_fpr, width, label="Brak kompensacji", color='#e74c3c', alpha=0.75, edgecolor='black', lw=0.8)
    b2 = ax1.bar(x + width/2, rot_fpr, width, label="Kompensacja SO(3)", color='#27ae60', alpha=0.85, edgecolor='black', lw=0.8)
    ax1.set_title("Odsetek błędów tła (FPR) dla r=3, τ=0.3 (mniej = lepiej)", fontsize=12, fontweight='bold', pad=10)
    ax1.set_ylabel("FPR tła [%]", fontsize=11)
    ax1.set_xticks(x)
    ax1.set_xticklabels(seq_labels, rotation=25, ha='right', fontsize=9.5)
    ax1.set_ylim(0, 55)
    ax1.set_yticks(np.arange(0, 56, 10))
    ax1.grid(True, linestyle="--", alpha=0.6, axis='y')
    ax1.legend(fontsize=9.5, loc="upper right")
    
    # Annotate reductions
    for i, (red, r_fpr) in enumerate(zip(reductions, rot_fpr)):
        if red >= 1.3:
            ax1.annotate(f"{red:.1f}× mniej\nbłędów", (x[i] + width/2, r_fpr + 0.8),
                         ha='center', va='bottom', fontsize=8.5, fontweight='bold', color='#1e824c')
        else:
            ax1.annotate(f"{red:.2f}×", (x[i] + width/2, r_fpr + 0.8),
                         ha='center', va='bottom', fontsize=8.5, color='#555555')
                         
    # 2. Precyzja (Precision)
    b3 = ax2.bar(x - width/2, none_prec, width, label="Brak kompensacji", color='#95a5a6', alpha=0.75, edgecolor='black', lw=0.8)
    b4 = ax2.bar(x + width/2, rot_prec, width, label="Kompensacja SO(3)", color='#2980b9', alpha=0.85, edgecolor='black', lw=0.8)
    ax2.set_title("Precyzja detekcji (Precision) dla r=3, τ=0.3 (więcej = lepiej)", fontsize=12, fontweight='bold', pad=10)
    ax2.set_ylabel("Precyzja [%]", fontsize=11)
    ax2.set_xticks(x)
    ax2.set_xticklabels(seq_labels, rotation=25, ha='right', fontsize=9.5)
    ax2.set_ylim(0, 105)
    ax2.set_yticks(np.arange(0, 106, 20))
    ax2.grid(True, linestyle="--", alpha=0.6, axis='y')
    ax2.legend(fontsize=9.5, loc="upper left")
    
    # Annotate delta precision
    for i, (d, r_p) in enumerate(zip(deltas, rot_prec)):
        if d >= 1.0:
            ax2.annotate(f"+{d:.1f} pp", (x[i] + width/2, r_p + 1.2),
                         ha='center', va='bottom', fontsize=8.5, fontweight='bold', color='#1f4e79')
        elif d <= -0.5:
            ax2.annotate(f"{d:.1f} pp", (x[i] + width/2, r_p + 1.2),
                         ha='center', va='bottom', fontsize=8.5, color='#c0392b')
                         
    fig.suptitle("Ewaluacja Kompensacji Rotacji SO(3) na 7 sekwencjach EVIMO2 (Punkt Zbalansowany: r=3, τ=0.3)",
                 fontsize=14, fontweight='bold', y=0.99)
    fig.tight_layout()
    out_path2 = EVAL_DIR / "operating_point_tau03_summary.png"
    fig.savefig(out_path2, bbox_inches="tight", dpi=180)
    plt.close(fig)
    print(f"Saved {out_path2}")
    shutil.copy(out_path2, ARTIFACT_DIR / "operating_point_tau03_summary.png")

if __name__ == "__main__":
    gd = load_data()
    generate_roc_tau03_plot(gd)
    generate_metrics_summary_tau03(gd)
