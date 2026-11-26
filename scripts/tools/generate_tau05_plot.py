#!/usr/bin/env python3
"""Generate dedicated comparative analysis bar plot specifically for r=3, tau=0.5 (High Precision)."""
import csv
import shutil
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

EVAL_DIR = Path("results/v1_rotcomp/benchmark")
ARTIFACT_DIR = Path("/home/pawsooon/.gemini/antigravity/brain/5814172c-e6cf-4394-ae23-449e5e16c752")

def generate_metrics_summary_tau05():
    qual_csv = EVAL_DIR / "quality_metrics_7seqs.csv"
    if not qual_csv.exists():
        print(f"Error: {qual_csv} not found")
        return
        
    rows = []
    with open(qual_csv) as f:
        r = csv.DictReader(f)
        for row in r:
            if "tau=0.5" in row["point"]:
                rows.append(row)
                
    if not rows:
        print("Error: No tau=0.5 rows found")
        return
        
    seq_labels = [r["sequence"].replace("_000000", "").replace("_dyn_test", "").replace("_dyn_train", "") for r in rows]
    none_fpr = [float(r["none_fpr"])*100 for r in rows]
    rot_fpr = [float(r["rot_fpr"])*100 for r in rows]
    none_prec = [float(r["none_precision"])*100 for r in rows]
    rot_prec = [float(r["rot_precision"])*100 for r in rows]
    none_rec = [float(r["none_recall"])*100 for r in rows]
    rot_rec = [float(r["rot_recall"])*100 for r in rows]
    reductions = [float(r["fpr_reduction_x"]) for r in rows]
    deltas = [round((float(r["rot_precision"]) - float(r["none_precision"]))*100, 1) for r in rows]
    
    x = np.arange(len(seq_labels))
    width = 0.36
    
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17, 5.8), dpi=160)
    
    # 1. Błędy tła (FPR)
    b1 = ax1.bar(x - width/2, none_fpr, width, label="Brak kompensacji", color='#e74c3c', alpha=0.75, edgecolor='black', lw=0.8)
    b2 = ax1.bar(x + width/2, rot_fpr, width, label="Kompensacja SO(3)", color='#27ae60', alpha=0.85, edgecolor='black', lw=0.8)
    ax1.set_title("Odsetek błędów tła (FPR) dla r=3, τ=0.5 (mniej = lepiej)", fontsize=12, fontweight='bold', pad=10)
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
            ax1.annotate(f"{red:.2f}× mniej\nbłędów", (x[i] + width/2, r_fpr + 0.6),
                         ha='center', va='bottom', fontsize=8.5, fontweight='bold', color='#1e824c')
        else:
            ax1.annotate(f"{red:.2f}×", (x[i] + width/2, r_fpr + 0.6),
                         ha='center', va='bottom', fontsize=8.5, color='#555555')
                         
    # 2. Precyzja (Precision)
    b3 = ax2.bar(x - width/2, none_prec, width, label="Brak kompensacji", color='#95a5a6', alpha=0.75, edgecolor='black', lw=0.8)
    b4 = ax2.bar(x + width/2, rot_prec, width, label="Kompensacja SO(3)", color='#2980b9', alpha=0.85, edgecolor='black', lw=0.8)
    ax2.set_title("Precyzja detekcji (Precision) dla r=3, τ=0.5 (więcej = lepiej)", fontsize=12, fontweight='bold', pad=10)
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
        elif d > 0.0:
            ax2.annotate(f"+{d:.1f} pp", (x[i] + width/2, r_p + 1.2),
                         ha='center', va='bottom', fontsize=8.5, color='#1f4e79')
        elif d <= -0.1:
            ax2.annotate(f"{d:.1f} pp", (x[i] + width/2, r_p + 1.2),
                         ha='center', va='bottom', fontsize=8.5, color='#c0392b')
                         
    fig.suptitle("Ewaluacja Kompensacji Rotacji SO(3) na 7 sekwencjach EVIMO2 (Punkt Wysokiej Precyzji: r=3, τ=0.5)",
                 fontsize=14, fontweight='bold', y=0.99)
    fig.tight_layout()
    out_path = EVAL_DIR / "operating_point_tau05_summary.png"
    fig.savefig(out_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    print(f"Saved {out_path}")
    shutil.copy(out_path, ARTIFACT_DIR / "operating_point_tau05_summary.png")
    print(f"Copied to {ARTIFACT_DIR / 'operating_point_tau05_summary.png'}")

if __name__ == "__main__":
    generate_metrics_summary_tau05()
