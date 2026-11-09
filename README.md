# evimo2_tools

Narzędzia w Pythonie i C++ do inspekcji, weryfikacji spójności, kompensacji ruchu i detekcji obiektów ruchomych (IMO) dla sekwencji **EVIMO2v2** (Prophesee, kamera lewa).
Projekt stanowi bazę badawczą do pracy magisterskiej: *"Event camera and IMU data fusion for the detection of independently moving objects"*.

Wymagania środowiskowe: Python >= 3.9 (zalecane środowisko conda `ml` lub `embedded-ai`), `numpy`, `matplotlib`, `pytest`.

---

## Struktura projektu

* `data/` – sekwencje danych wejściowych EVIMO2 (w `.gitignore`).
  * `data/raw/train/` oraz `data/raw/eval/` – wygodny symlink do `npz_left_camera_imo/left_camera/imo/`.
* `docs/` – dokumentacja badawcza i notatki:
  * [`docs/IMPORTANT.md`](docs/IMPORTANT.md) – skonsolidowane ustalenia i wnioski eksperymentalne (konwencje poz, parametry kalibracji, jednostki głębi, tablica index, selekcja sekwencji).
  * [`docs/notes.txt`](docs/notes.txt) – robocze notatki z inspekcji sekwencji.
* `results/` – wyniki działania algorytmów (pliki JSON, wizualizacje PNG, maski):
  * `results/summaries/` – tabele zbiorcze (`survey.csv`, `imo_labels.csv`, `object_motion.csv`).
  * `results/<nazwa_sekwencji>/` – wyniki per sekwencja (`gyro_vs_gt.json`, `rotcomp_eval.json`, `imo_detect.json`, etc.).
  * `out/` – symlink do `results/` zapewniający pełną kompatybilność wsteczną.
* `scripts/` – skrypty narzędziowe i algorytmy w Pythonie:
  * `rotcomp.py` / `rotcomp_eval.py` – kompensacja rotacji z IMU ($SO(3)$, RadTan) i jej ewaluacja.
  * `gyro_vs_gt.py` – kalibracja żyroskopu względem GT (ekstrynsyki $R$, bias, opóźnienie).
  * `imo_detect.py` – prototyp detektora IMO v0 (mapa niespójności, klasteryzacja).
  * `gt_motion_viz.py` – dwupanelowa wizualizacja czystego Ground Truth (obiekty z prędkościami i prawdziwe zdarzenia TP IMO).
  * `imo_labels.py` / `object_motion.py` – etykietowanie dynamiki obiektów w czasie i wyznaczanie prędkości w układzie świata.
  * `depth_check.py` / `check_sync.py` / `inspect_evimo2.py` – weryfikacja głębi i synchronizacji.
  * `survey_sequences.py` – przegląd statystyk sekwencji EVIMO2.
* `cpp/` – kod w C++:
  * `cpp/viewer/` – przeglądarka sekwencji EVIMO2 (`main.cpp`, `events.hpp`, `viz.hpp`).
* `tests/` – testy formatu i spójności danych (`pytest tests/`).

---

## Szybki start

### Aktywacja środowiska:
```bash
conda activate ml
```

### Uruchomienie testów spójności formatu sekwencji:
```bash
EVIMO2_SEQ_DIR=data/raw/train/scene14_dyn_test_01_000000 pytest
```

### Przykłady uruchomienia narzędzi:
```bash
# 1. Sprawdzenie głębi i jednostek
python3 scripts/depth_check.py data/raw/train/scene6_dyn_train_03_000000

# 2. Ewaluacja kompensacji rotacji
python3 scripts/rotcomp_eval.py data/raw/eval/scene14_dyn_test_05_000000 --max-windows 10

# 3. Uruchomienie detekcji IMO v0
python3 scripts/imo_detect.py data/raw/eval/scene14_dyn_test_05_000000 --max-windows 10 --viz 5

# 4. Wizualizacja czystego Ground Truth ruchu
python3 scripts/gt_motion_viz.py data/raw/train/scene6_dyn_train_03_000000 --max-viz 5
```
