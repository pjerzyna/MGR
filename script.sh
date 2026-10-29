python3 - <<'EOF'
import glob, json, os
import numpy as np

rows = []
for info in glob.glob("data/npz_left_camera_imo/left_camera/imo/**/dataset_info.npz", recursive=True):
    seq = os.path.basename(os.path.dirname(info))
    j = f"out/{seq}/gyro_vs_gt.json"
    if not os.path.exists(j):
        continue
    with np.load(info, allow_pickle=True) as z:
        ros = float(z["meta"].item()["meta"]["ros_time_offset"])
    d = json.load(open(j))
    rows.append((seq, ros, d["lag_s"] * 1e3, d["residual_rel"] * 100))
rows.sort(key=lambda r: r[1])

t0 = rows[0][1]
print(f"{'sekwencja':28s} {'czas nagrania [s]':>18s} {'d [ms]':>8s} {'res [%]':>8s}")
for seq, ros, lag, res in rows:
    print(f"{seq:28s} {ros - t0:18.1f} {lag:8.1f} {res:8.1f}")

print("\nNachylenie d(t) per scena (dla sekwencji z residuum < 15%):")
for scene in sorted({r[0].split('_')[0] for r in rows}):
    pts = [(r[1], r[2]) for r in rows if r[0].startswith(scene + "_") and r[3] < 15]
    if len(pts) < 3:
        print(f"  {scene}: za mało punktów ({len(pts)})"); continue
    t = np.array([p[0] for p in pts]) - pts[0][0]
    d = np.array([p[1] for p in pts])
    slope, icpt = np.polyfit(t, d, 1)
    rms = np.sqrt(np.mean((d - (slope * t + icpt)) ** 2))
    print(f"  {scene}: {slope * 1e3:+.1f} ppm, rozrzut wokół prostej {rms:.2f} ms, n={len(pts)}")
EOF