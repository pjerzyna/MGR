#!/usr/bin/env bash
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

ROOT="$ROOT_DIR/data/raw"
find_seq() { find -L "$ROOT" -maxdepth 2 -name "$1" -type d | head -1; }

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
  
  python3 "$SCRIPT_DIR/depth_check.py" "$seq_path" --max-frames 60
  echo ""
done