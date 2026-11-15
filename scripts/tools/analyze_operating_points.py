#!/usr/bin/env python3
"""Detailed mathematical and visual analysis of the operating point selection (r, tau)."""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

EVAL_DIR = Path("results/v1_rotcomp/benchmark")
df = pd.read_csv(EVAL_DIR / "param_sweep_grid.csv")

# 1. Agregacja globalna
grouped = df.groupby(['blur', 'tau']).agg({
    'rot_recall': 'mean',
    'rot_fpr': 'mean',
    'rot_precision': 'mean',
    'rot_iou': 'mean',
}).reset_index()

grouped['youden_J'] = grouped['rot_recall'] - grouped['rot_fpr']
grouped['f1'] = 2 * (grouped['rot_precision'] * grouped['rot_recall']) / (grouped['rot_precision'] + grouped['rot_recall'] + 1e-6)
grouped['dist_ideal'] = np.sqrt(grouped['rot_fpr']**2 + (1 - grouped['rot_recall'])**2)

grouped = grouped.round(4)
grouped.to_csv(EVAL_DIR / "operating_points_aggregate.csv", index=False)
print("Saved operating_points_aggregate.csv")

# 2. Wykres wizualizacyjny Pareto & Heatmapy
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), dpi=150)

# Panel 1: Front Pareto (FPR vs Recall dla różnych r i tau)
markers = {2: '^', 3: 'o', 5: 's'}
colors = {2: '#e67e22', 3: '#2980b9', 5: '#8e44ad'}

for r in [2, 3, 5]:
    sub = grouped[grouped['blur'] == r].sort_values('rot_fpr')
    ax1.plot(sub['rot_fpr'] * 100, sub['rot_recall'] * 100, f'-{markers[r]}',
             color=colors[r], lw=2, ms=8, label=f"Promień r={r} (filtr {2*r+1}x{2*r+1} px)")
    for _, row in sub.iterrows():
        ax1.annotate(f"tau={row['tau']}", (row['rot_fpr']*100 + 0.6, row['rot_recall']*100 - 1.0),
                     fontsize=8.5, color=colors[r])

# Wyróżnienie dwóch kluczowych punktów
pt_balanced = grouped[(grouped['blur'] == 3) & (grouped['tau'] == 0.3)].iloc[0]
pt_high_prec = grouped[(grouped['blur'] == 3) & (grouped['tau'] == 0.5)].iloc[0]

ax1.scatter([pt_balanced['rot_fpr']*100], [pt_balanced['rot_recall']*100],
            s=220, facecolors='none', edgecolors='#27ae60', lw=2.5, zorder=5)
ax1.annotate("Punkt Zbalansowany\n(r=3, tau=0.3)\nRecall=49%, FPR=16.7%",
             xy=(pt_balanced['rot_fpr']*100, pt_balanced['rot_recall']*100),
             xytext=(pt_balanced['rot_fpr']*100 - 9, pt_balanced['rot_recall']*100 + 8),
             arrowprops=dict(facecolor='#27ae60', arrowstyle='->', lw=1.5),
             fontsize=9, fontweight='bold', bbox=dict(boxstyle='round,pad=0.3', facecolor='#eafaf1', edgecolor='#27ae60'))

ax1.scatter([pt_high_prec['rot_fpr']*100], [pt_high_prec['rot_recall']*100],
            s=220, facecolors='none', edgecolors='#c0392b', lw=2.5, zorder=5)
ax1.annotate("Punkt Wysokiej Precyzji\n(r=3, tau=0.5)\nRecall=26%, FPR=6.9%",
             xy=(pt_high_prec['rot_fpr']*100, pt_high_prec['rot_recall']*100),
             xytext=(pt_high_prec['rot_fpr']*100 + 3, pt_high_prec['rot_recall']*100 - 8),
             arrowprops=dict(facecolor='#c0392b', arrowstyle='->', lw=1.5),
             fontsize=9, fontweight='bold', bbox=dict(boxstyle='round,pad=0.3', facecolor='#fdedec', edgecolor='#c0392b'))

ax1.set_title("1. Front Pareto w przestrzeni ROC (FPR vs Recall)", fontsize=12, fontweight='bold')
ax1.set_xlabel("FPR tła [%] (mniej = lepiej)", fontsize=11)
ax1.set_ylabel("Czułość (Recall) [%] (więcej = lepiej)", fontsize=11)
ax1.legend(loc='lower right', fontsize=9.5)
ax1.set_xlim(0, 45)
ax1.set_ylim(10, 85)

# Panel 2: Heatmapa Wskaźnika Youdena J = Recall - FPR
pivot_j = grouped.pivot(index='blur', columns='tau', values='youden_J')
im = ax2.imshow(pivot_j.values, cmap='Blues', aspect='auto')

ax2.set_xticks(range(len(pivot_j.columns)))
ax2.set_xticklabels([f"{t:.1f}" for t in pivot_j.columns], fontsize=10)
ax2.set_yticks(range(len(pivot_j.index)))
ax2.set_yticklabels([f"r={b}" for b in pivot_j.index], fontsize=10)

ax2.set_title("2. Jakość kompromisu: Wskaźnik Youdena (J = Recall - FPR)", fontsize=12, fontweight='bold')
ax2.set_xlabel("Próg niespójności tau [-]", fontsize=11)
ax2.set_ylabel("Promień wygładzania r [px]", fontsize=11)

for i in range(len(pivot_j.index)):
    for j in range(len(pivot_j.columns)):
        val = pivot_j.values[i, j]
        ax2.text(j, i, f"J = {val:.3f}", ha='center', va='center',
                 color='white' if val > 0.3 else 'black', fontweight='bold', fontsize=10)

plt.colorbar(im, ax=ax2, label="Wskaźnik Youdena J [-]")
plt.suptitle("Metodologia i dowód doboru punktów pracy (r, tau) na 7 sekwencjach EVIMO2", fontsize=13, fontweight='bold', y=0.98)
plt.tight_layout()

fig.savefig(EVAL_DIR / "operating_points_analysis.png", bbox_inches='tight', dpi=180)
fig.savefig(Path("/home/pawsooon/.gemini/antigravity/brain/5814172c-e6cf-4394-ae23-449e5e16c752/operating_points_analysis.png"), bbox_inches='tight', dpi=180)
plt.close(fig)
print("Saved operating_points_analysis.png")
