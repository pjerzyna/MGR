#!/usr/bin/env bash
# Konwertuje obrazy PPM z imo_detect do PNG (ImageMagick lub Python Pillow)
#
# Użycie:
#   ./convert_to_png.sh <katalog_z_obrazami> [--force]
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Użycie: $0 <katalog_z_obrazami> [--force]"
    exit 1
fi

DIR="$1"
FORCE=0
if [[ $# -ge 2 && "$2" == "--force" ]]; then
    FORCE=1
fi

if [[ ! -d "$DIR" ]]; then
    echo "Błąd: Katalog $DIR nie istnieje!" >&2
    exit 1
fi

COUNT=0
CONV_CMD=""

if command -v convert >/dev/null 2>&1; then
    CONV_CMD="convert"
elif command -v magick >/dev/null 2>&1; then
    CONV_CMD="magick"
fi

for f in "$DIR"/*_imo.ppm; do
    [[ -e "$f" ]] || continue
    png="${f%.ppm}.png"
    if [[ -e "$png" && $FORCE -eq 0 ]]; then
        continue
    fi
    if [[ -n "$CONV_CMD" ]]; then
        "$CONV_CMD" "$f" "$png"
    else
        python3 -c "from PIL import Image; Image.open('$f').save('$png')"
    fi
    COUNT=$((COUNT + 1))
done

echo "Skonwertowano $COUNT obrazów PPM do PNG w katalogu: $DIR"
