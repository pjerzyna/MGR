#!/usr/bin/env bash
export EVIMO2_OUT=/home/pawsooon/Documents/Vsc/npz_odczyty/evimo2_viewer/out

python3 imo_labels.py --all data/npz_left_camera_imo/left_camera/imo --validate --csv imo_labels.csv
#python3 imo_labels.py data/npz_left_camera_imo/left_camera/imo/train/scene6_dyn_train_03_000000 --validate --save-masks