#!/usr/bin/env bash
set -u
cd "$(dirname "$0")"

ROOT=data/npz_left_camera_imo/left_camera/imo
find_seq() { find "$ROOT" -maxdepth 2 -name "$1" -type d | head -1; }

# Lista sekwencji do sprawdzenia
SEKWENCJE=(
  "scene6_dyn_train_03_000000"
  "scene14_dyn_test_05_000000"
)

for s in "${SEKWENCJE[@]}"; do
  seq_path="$(find_seq "$s")"
  
  if [ -z "$seq_path" ]; then
    echo " Nie znaleziono katalogu dla: $s"
    continue
  fi

  echo "============================================================"
  echo "Sprawdzanie głębi dla: $s"
  echo "Ścieżka: $seq_path"
  echo "============================================================"
  
  # Wywołanie depth_check.py (opcjonalnie można dodać np. --max-frames 60)
  python3 depth_check.py "$seq_path" --max-frames 60
  echo ""
done