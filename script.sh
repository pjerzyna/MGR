#!/usr/bin/env bash
set -u
cd "$(dirname "$0")"
export EVIMO2_OUT=/home/pawsooon/Documents/Vsc/npz_odczyty/evimo2_viewer/out
ROOT=data/npz_left_camera_imo/left_camera/imo

for s in scene6_dyn_train_03_000000 scene14_dyn_test_05_000000 scene14_dyn_test_04_000000 scene13_dyn_test_01_000000; do
  d=$(find $ROOT -maxdepth 2 -name "$s" -type d | head -1)
  python3 imo_detect.py "$d"
done

d=$(find $ROOT -maxdepth 2 -name scene14_dyn_test_05_000000 -type d | head -1)
python3 imo_detect.py "$d" --grid