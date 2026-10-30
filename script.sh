#!/usr/bin/env bash
set -u
cd "$(dirname "$0")"          # zawsze uruchamiaj z katalogu, w którym leży script.sh 
export EVIMO2_OUT=/home/pawsooon/Documents/Vsc/npz_odczyty/evimo2_viewer/out

for s in scene13_dyn_test_02_000000 scene14_dyn_test_05_000000 scene15_dyn_test_03_000000; do
  d=$(find data/npz_left_camera_imo/left_camera/imo -maxdepth 2 -name "$s" -type d | head -1)
  python3 event_lag_check.py "$d" --window-ms 100 --windows 20 --min-omega 0.3
done
python3 event_lag_check.py --summary