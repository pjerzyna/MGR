# evimo2_tools: Detekcja Obiektów Ruchomych (IMO) z Kamer Zdarzeniowych i IMU

Narzędzia w Pythonie i deterministyczny silnik C++ do kompensacji ruchu kamery ($SO(3)$) oraz detekcji niezależnie poruszających się obiektów (IMO) dla sekwencji **EVIMO2v2** (czujnik neuromorficzny Prophesee Gen3 VGA).
Projekt stanowi bazę badawczo-inżynierską do pracy magisterskiej: *"Event camera and IMU data fusion for the detection of independently moving objects"*.

Wymagania środowiskowe: Python >= 3.9 (dedykowane środowisko conda `mgr`), `numpy`, `matplotlib`, `pytest`, kompilator C++17 (GCC/Clang).

---

## 1. Struktura Projektu

```text
├── scripts/
│   ├── run_evaluation.py        # Główny punkt wejścia v1: automatyczna ewaluacja 7 sekwencji
│   ├── run_epipolar_foe.py      # Główny punkt wejścia v2: filtr geometrii epipolarnej i FOE
│   └── tools/                   # Narzędzia diagnostyczne, kalibracyjne i przygotowania danych
│       ├── imo_detect.py        # Referencyjny detektor okienkowy w Pythonie
│       ├── generate_eval_reports.py  # Generator raportów CSV/Markdown i wykresów
│       ├── analyze_operating_points.py # Analiza Pareto i indeksu Youdena
│       ├── inspect_evimo2.py    # Narzędzie inspekcji sekwencji
│       └── ...
│
├── python/
│   └── evimo2/                  # Kanoniczna biblioteka algorytmów (Python package)
│       ├── calibration.py       # Kalibracja IMU, opóźnienia czasowe i macierze rotacji
│       ├── common.py            # Ładowanie zdarzeń i metadanych z EVIMO2
│       ├── detector.py          # Obliczanie residual map, SAT box blur i klastrowanie CCL
│       ├── epipolar.py          # Geometria epipolarna, FOE i 2D NCC filtr klastrów
│       ├── geometry.py          # Kwaterniony, RadTan undistort, rzutowania pinhole
│       ├── kinematics.py        # Trajektorie Ground Truth, prędkości obiektów i poz
│       └── rotcomp.py           # Kompensacja rotacji SO(3) ze splattingiem biliniowym IWE
│
├── cpp/
│   ├── include/                 # Nagłówki silnika C++ (rotcomp, detector, geometry, events)
│   ├── apps/
│   │   ├── imo_detect/          # Deterministyczny potok C++ (przepustowość >4 Mev/s)
│   │   └── viewer/              # Szybka przeglądarka sekwencji EVIMO2
│   ├── tests/                   # Rygorystyczne testy numerycznej zgodności z Pythonem
│   └── Makefile                 # Budowanie binariów z flagami -O3 -std=c++17
│
├── results/
│   ├── v1_rotcomp/              # Rozdział I: Kompensacja rotacji SO(3) z żyroskopu
│   │   ├── benchmark/           # Zbiorcze wyniki ewaluacji C++ i Pythona (CSV, wykresy, raport)
│   │   ├── per_sequence/        # Szczegółowe metryki i analizy per sekwencja
│   │   └── visual_samples/      # Wyrenderowane ramki detekcji PNG/PPM
│   │
│   ├── v2_epipolar_foe/         # Rozdział II: Filtr geometrii epipolarnej (Focus of Expansion)
│   │   ├── scene13_01_metrics.csv # Wyniki ewaluacji (redukcja FP o 5.42x, podwojenie F1)
│   │   ├── epipolar_v2_validation_report.md
│   │   └── epipolar_v2_comparison_frame60.png
│   │
│   └── archive/                 # Eksperymenty pomocnicze (np. Oracle GT depth)
│
├── tests/                       # Testy jednostkowe pytest (12/12 testów zaliczonych)
└── scratch/                     # Pliki robocze i jednorazowe eksperymenty (git-ignored)
```

---

## 2. Szybki Start i Przykłady Użycia

### Aktywacja środowiska Conda
```bash
conda activate mgr
```

### Uruchomienie testów jednostkowych
```bash
pytest tests/
./cpp/bin/test_numerics
```

### 1. Zbiorcza ewaluacja detektora v1 na całym benchmarku EVIMO2 (7 sekwencji):
```bash
python scripts/run_evaluation.py            # Pełny potok end-to-end (C++ benchmark + jakość GT + wykresy)
python scripts/run_evaluation.py --plots    # Szybkie wygenerowanie wykresów z istniejących wyników
python scripts/run_evaluation.py --bench    # Tylko benchmark wydajnościowy C++ (throughput, latency)
```

### 2. Uruchomienie ewaluacji detektora v2 (Filtr Epipolarny FOE):
```bash
python scripts/run_evaluation_v2.py         # Pełna ewaluacja v1 vs v2 na 7 sekwencjach (wykresy, benchmark, ROC)
python scripts/run_evaluation_v2.py --bench # Tylko benchmark C++ v2 (throughput, latency, narzut FOE)
python scripts/run_epipolar_foe.py          # Pojedynczy test i wizualizacja na sekwencji z translacją (scene13_01)
```

### 3. Kompilacja i bezpośrednie wywołanie silnika C++:
```bash
make -C cpp

# Detektor bazowy v1 (RotComp):
./cpp/bin/imo_detect data/npz_left_camera_imo/left_camera/imo/train/scene6_dyn_train_03_000000 \
    --window-ms 40 --step-ms 40 --blur 3 --tau 0.5 --lag-ms 14.3

# Detektor v2 z aktywnym filtrem epipolarnym FOE (odrzucanie paralaksy translacji):
./cpp/bin/imo_detect data/npz_left_camera_imo/left_camera/imo/train/scene13_dyn_test_01_000000 \
    --epipolar --epipolar-thresh 35.0 --blur 3 --tau 0.3 --max-windows 50
```
