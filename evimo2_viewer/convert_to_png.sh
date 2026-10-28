#!/usr/bin/env bash
# Konwertuje obrazy .ppm/.pgm z evimo2_viewer do PNG (ImageMagick).
#
# Użycie:
#   ./convert_to_png.sh <folder> [--type pol|ts|all] [--force]
#
#   <folder>   nazwa folderu z obrazami; jeśli taki katalog nie istnieje wprost,
#              szukany jest w out/<folder>  (np. "scene15_06" -> out/scene15_06)
#   --type     pol = tylko frame_*_pol.ppm, ts = tylko frame_*_ts.pgm, all = oba (domyślnie)
#   --force    nadpisuj istniejące PNG (domyślnie istniejące są pomijane)
#
# Wynik: <folder>/png_files/pol/frame_XXXX_pol.png
#        <folder>/png_files/ts/frame_XXXX_ts.png
# Oryginały .ppm/.pgm nie są ruszane.
set -euo pipefail

usage() { sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-1}"; }

folder=""; type="all"; force=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --type)  [[ $# -ge 2 ]] || usage; type="$2"; shift 2 ;;
        --force) force=1; shift ;;
        -h|--help) usage 0 ;;
        -*) echo "Nieznana flaga: $1" >&2; usage ;;
        *)  [[ -z "$folder" ]] || { echo "Podano więcej niż jeden folder" >&2; usage; }
            folder="$1"; shift ;;
    esac
done
[[ -n "$folder" ]] || usage
[[ "$type" == "pol" || "$type" == "ts" || "$type" == "all" ]] || { echo "--type: pol|ts|all" >&2; exit 1; }

# Rozwiąż folder: wprost albo w out/
base="$(cd "$(dirname "$0")" && pwd)"
if   [[ -d "$folder" ]];          then src="$folder"
elif [[ -d "$base/out/$folder" ]]; then src="$base/out/$folder"
else echo "Nie znaleziono katalogu: $folder (ani out/$folder)" >&2; exit 1; fi
src="${src%/}"

# ImageMagick 7 ("magick") albo 6 ("convert")
if   command -v magick  >/dev/null 2>&1; then im=(magick)
elif command -v convert >/dev/null 2>&1; then im=(convert)
else echo "Brak ImageMagick (sudo apt install imagemagick)" >&2; exit 1; fi

files=()
[[ "$type" == "pol" || "$type" == "all" ]] && files+=("$src"/frame_*_pol.ppm)
[[ "$type" == "ts"  || "$type" == "all" ]] && files+=("$src"/frame_*_ts.pgm)

out="$src/png_files"
done_n=0; skip_n=0; found=0
for f in "${files[@]}"; do
    [[ -e "$f" ]] || continue            # wzorzec bez dopasowań zostaje dosłownie
    found=$((found + 1))
    case "$f" in
        *_pol.ppm) sub="pol" ;;
        *_ts.pgm)  sub="ts" ;;
    esac
    png="$out/$sub/$(basename "${f%.*}").png"
    if [[ -e "$png" && $force -eq 0 ]]; then skip_n=$((skip_n + 1)); continue; fi
    mkdir -p "$out/$sub"
    "${im[@]}" "$f" "$png"
    done_n=$((done_n + 1))
done

if [[ $found -eq 0 ]]; then
    echo "Brak plików frame_*_{pol.ppm,ts.pgm} w $src" >&2
    exit 1
fi
echo "Gotowe: skonwertowano $done_n, pominięto $skip_n (istniejące) -> $out/{pol,ts}"