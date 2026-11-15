#!/usr/bin/env python3
"""Generate comprehensive evaluation reports, CSVs and plots from benchmark and grid results."""
import os
import json
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

EVAL_DIR = Path("results/v1_rotcomp/benchmark")

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
            # Try alternate name if any
            alt = list(EVAL_DIR.glob(f"*{seq_name}*/*.json"))
            if alt:
                json_file = alt[0]
            else:
                print(f"Warning: {json_file} does not exist yet.")
                continue
        with open(json_file) as f:
            grid_data[seq_name] = json.load(f)
    return grid_data

def generate_csvs(grid_data):
    # 1. Parameter sweep grid CSV
    grid_rows = []
    for seq_name, cat_desc, cat_group in SEQUENCES:
        if seq_name not in grid_data:
            continue
        d = grid_data[seq_name]
        modes = d.get("modes", {})
        
        # Identify mode keys
        none_key = [k for k in modes.keys() if "brak" in k.lower()]
        rot_key = [k for k in modes.keys() if "stałe" in k.lower() or "oracle" in k.lower()]
        
        m_none = modes[none_key[0]] if none_key else {}
        m_rot = modes[rot_key[0]] if rot_key else {}
        
        for blur in [2, 3, 5]:
            for tau in [0.2, 0.3, 0.4, 0.5]:
                pk = f"blur{blur}_tau{tau}"
                row = {
                    "sequence": seq_name,
                    "category": cat_group,
                    "description": cat_desc,
                    "blur": blur,
                    "tau": tau,
                }
                
                # None mode
                if pk in m_none and "all" in m_none[pk]:
                    v = m_none[pk]["all"]
                    row["none_recall"] = round(v.get("recall", float("nan")), 3)
                    row["none_precision"] = round(v.get("precision", float("nan")), 3)
                    row["none_iou"] = round(v.get("iou", float("nan")), 3)
                    row["none_fpr"] = round(v.get("fpr", float("nan")), 4)
                
                # Rotcomp mode
                if pk in m_rot and "all" in m_rot[pk]:
                    v = m_rot[pk]["all"]
                    row["rot_recall"] = round(v.get("recall", float("nan")), 3)
                    row["rot_precision"] = round(v.get("precision", float("nan")), 3)
                    row["rot_iou"] = round(v.get("iou", float("nan")), 3)
                    row["rot_fpr"] = round(v.get("fpr", float("nan")), 4)
                    
                grid_rows.append(row)
                
    if grid_rows:
        fieldnames = list(grid_rows[0].keys())
        with open(EVAL_DIR / "param_sweep_grid.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(grid_rows)
        print(f"Saved {EVAL_DIR / 'param_sweep_grid.csv'}")

    # 2. Quality metrics at standard operating points (blur=3, tau=0.3 and blur=3, tau=0.5)
    quality_rows = []
    for seq_name, cat_desc, cat_group in SEQUENCES:
        if seq_name not in grid_data:
            continue
        d = grid_data[seq_name]
        modes = d.get("modes", {})
        none_key = [k for k in modes.keys() if "brak" in k.lower()]
        rot_key = [k for k in modes.keys() if "stałe" in k.lower() or "oracle" in k.lower()]
        
        m_none = modes[none_key[0]] if none_key else {}
        m_rot = modes[rot_key[0]] if rot_key else {}
        
        for pk, pt_name in [("blur3_tau0.3", "Zbalansowany (r=3, tau=0.3)"), ("blur3_tau0.5", "Wysoka precyzja (r=3, tau=0.5)")]:
            if pk in m_none and pk in m_rot:
                vn = m_none[pk]["all"]
                vr = m_rot[pk]["all"]
                
                fpr_red = (vn["fpr"] / vr["fpr"]) if vr["fpr"] > 1e-6 else float("nan")
                delta_prec = (vr["precision"] - vn["precision"]) * 100.0 if np.isfinite(vr["precision"]) and np.isfinite(vn["precision"]) else float("nan")
                
                quality_rows.append({
                    "sequence": seq_name,
                    "category": cat_group,
                    "point": pt_name,
                    "none_recall": round(vn.get("recall", 0), 3),
                    "none_fpr": round(vn.get("fpr", 0), 4),
                    "none_precision": round(vn.get("precision", 0), 3),
                    "none_iou": round(vn.get("iou", 0), 3),
                    "rot_recall": round(vr.get("recall", 0), 3),
                    "rot_fpr": round(vr.get("fpr", 0), 4),
                    "rot_precision": round(vr.get("precision", 0), 3),
                    "rot_iou": round(vr.get("iou", 0), 3),
                    "fpr_reduction_x": round(fpr_red, 2),
                    "delta_precision_pp": round(delta_prec, 1)
                })
                
    if quality_rows:
        fieldnames = list(quality_rows[0].keys())
        with open(EVAL_DIR / "quality_metrics_7seqs.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(quality_rows)
        print(f"Saved {EVAL_DIR / 'quality_metrics_7seqs.csv'}")

def generate_plots(grid_data):
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    # 1. ROC-like curves (Recall vs FPR) across 3 categories (for blur=3, tau in [0.2, 0.3, 0.4, 0.5])
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), dpi=150)
    
    categories = [
        ("Kategoria I: Spoczynek", ["scene6_dyn_train_03_000000"], axes[0]),
        ("Kategoria II: Dominujący obrót", ["scene14_dyn_test_05_000000", "scene14_dyn_test_04_000000", "scene13_dyn_test_02_000000"], axes[1]),
        ("Kategoria III: Ruch złożony / paralaksa", ["scene13_dyn_test_01_000000", "scene15_dyn_test_02_000000", "scene15_dyn_test_05_000000"], axes[2]),
    ]
    
    colors = ['#1f77b4', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    taus = [0.2, 0.3, 0.4, 0.5]
    
    for cat_title, seq_list, ax in categories:
        ax.set_title(cat_title, fontsize=12, fontweight='bold', pad=10)
        ax.set_xlabel("FPR tła (odsetek fałszywych alarmów) [-]", fontsize=10)
        ax.set_ylabel("Czułość detekcji (Recall) [-]", fontsize=10)
        ax.grid(True, linestyle="--", alpha=0.6)
        
        for idx, s in enumerate(seq_list):
            if s not in grid_data:
                continue
            modes = grid_data[s].get("modes", {})
            none_k = [k for k in modes.keys() if "brak" in k.lower()][0]
            rot_k  = [k for k in modes.keys() if "stałe" in k.lower() or "oracle" in k.lower()][0]
            
            c = colors[idx % len(colors)]
            short_name = s.replace("_000000", "").replace("_dyn_test", "").replace("_dyn_train", "")
            
            # Extract points for blur=3
            rec_none, fpr_none = [], []
            rec_rot, fpr_rot = [], []
            for tau in taus:
                pk = f"blur3_tau{tau}"
                if pk in modes[none_k]:
                    rec_none.append(modes[none_k][pk]["all"]["recall"])
                    fpr_none.append(modes[none_k][pk]["all"]["fpr"])
                if pk in modes[rot_k]:
                    rec_rot.append(modes[rot_k][pk]["all"]["recall"])
                    fpr_rot.append(modes[rot_k][pk]["all"]["fpr"])
                    
            if fpr_none and rec_none:
                ax.plot(fpr_none, rec_none, '--o', color=c, alpha=0.5, label=f"{short_name} (Brak comp)")
            if fpr_rot and rec_rot:
                ax.plot(fpr_rot, rec_rot, '-s', color=c, lw=2, label=f"{short_name} (SO3 comp)")
                
        ax.legend(fontsize=8, loc="best")
        
    plt.suptitle("Krzywe Czułość vs FPR tła (ROC Trade-Off) dla r=3 na 7 sekwencjach EVIMO2", fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    fig.savefig(EVAL_DIR / "roc_curves_categories.png", bbox_inches="tight", dpi=180)
    plt.close(fig)
    print(f"Saved {EVAL_DIR / 'roc_curves_categories.png'}")
    
    # 2. C++ Throughput and Latency benchmark plot
    bench_csv = EVAL_DIR / "benchmark_cpp_7seqs.csv"
    if bench_csv.exists():
        seq_names = []
        throughputs = []
        latencies = []
        rt_factors = []
        with open(bench_csv) as f:
            r = csv.DictReader(f)
            for row in r:
                sn = row["sequence"].replace("_000000", "").replace("_dyn_test", "").replace("_dyn_train", "")
                seq_names.append(sn)
                throughputs.append(float(row["throughput_mev"]))
                latencies.append(float(row["avg_window_ms"]))
                rt_factors.append(float(row["rt_factor"]))
                
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), dpi=150)
        x = np.arange(len(seq_names))
        
        # Subplot 1: Throughput
        bars1 = ax1.bar(x, throughputs, color='#2980b9', alpha=0.85, width=0.55)
        ax1.axhline(6.0, color='#e74c3c', linestyle='--', lw=2, label='Fizyczny limit Prophesee VGA (~6.0 Mev/s)')
        ax1.set_title("Przepustowość silnika C++ (Throughput)", fontsize=12, fontweight='bold')
        ax1.set_ylabel("Przepustowość [Mev/s]", fontsize=11)
        ax1.set_xticks(x)
        ax1.set_xticklabels(seq_names, rotation=25, ha='right', fontsize=9.5)
        ax1.set_ylim(0, max(max(throughputs)*1.25, 7.0))
        for b in bars1:
            ax1.annotate(f"{b.get_height():.2f}\nMev/s", (b.get_x() + b.get_width()/2, b.get_height() + 0.15),
                         ha='center', va='bottom', fontsize=8.5, fontweight='bold')
        ax1.legend(loc="upper left", fontsize=9)
        
        # Subplot 2: Latency per 40ms window
        bars2 = ax2.bar(x, latencies, color='#27ae60', alpha=0.85, width=0.55)
        ax2.axhline(40.0, color='#e74c3c', linestyle='--', lw=2, label='Granica czasu rzeczywistego (okno 40 ms)')
        ax2.set_title("Średni czas obliczeń okna 40 ms w C++", fontsize=12, fontweight='bold')
        ax2.set_ylabel("Czas per okno [ms] (mniej = lepiej)", fontsize=11)
        ax2.set_xticks(x)
        ax2.set_xticklabels(seq_names, rotation=25, ha='right', fontsize=9.5)
        ax2.set_ylim(0, max(max(latencies)*1.2, 45.0))
        for b, rt in zip(bars2, rt_factors):
            ax2.annotate(f"{b.get_height():.1f} ms\n({rt:.1f}x RT)", (b.get_x() + b.get_width()/2, b.get_height() + 0.8),
                         ha='center', va='bottom', fontsize=8.5, fontweight='bold')
        ax2.legend(loc="upper left", fontsize=9)
        
        plt.suptitle("Wydajność deterministycznego potoku C++ na 7 sekwencjach EVIMO2 (1 rdzeń x86)", fontsize=13, fontweight='bold', y=0.98)
        plt.tight_layout()
        fig.savefig(EVAL_DIR / "throughput_latency_cpp.png", bbox_inches="tight", dpi=180)
        plt.close(fig)
        print(f"Saved {EVAL_DIR / 'throughput_latency_cpp.png'}")

if __name__ == "__main__":
    gd = load_data()
    generate_csvs(gd)
    generate_plots(gd)
