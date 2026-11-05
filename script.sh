#!/usr/bin/env bash
set -u
cd "$(dirname "$0")"
export EVIMO2_OUT=/home/pawsooon/Documents/Vsc/npz_odczyty/evimo2_viewer/out
ROOT=data/npz_left_camera_imo/left_camera/imo
find_seq() { find $ROOT -maxdepth 2 -name "$1" -type d | head -1; }

# podgląd: gdzie są fałszywe alarmy i brakujące zdarzenia
python3 imo_detect.py "$(find_seq scene6_dyn_train_03_000000)" --viz 6 | sed -n '1,9p'
python3 imo_detect.py "$(find_seq scene14_dyn_test_05_000000)" --viz 6 | sed -n '1,9p'

# siatka na sekwencji treningowej (14_05)
python3 imo_detect.py "$(find_seq scene14_dyn_test_05_000000)" --grid

# test przeniesienia parametrów (rozmycie 5, tau 0,3) na innych sekwencjach
for s in scene14_dyn_test_04_000000 scene13_dyn_test_01_000000 scene15_dyn_test_03_000000 scene6_dyn_train_03_000000; do
  python3 imo_detect.py "$(find_seq $s)" --blur 5 --tau 0.3
done