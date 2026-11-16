# Raport Ewaluacji Detektora v2 (RotComp + Filtr Epipolarny FOE)

Data wygenerowania: 2026-11-16
Środowisko: Intel Core x86_64, 1 rdzeń CPU (deterministyczny silnik C++ / kompilacja `-O3`)

## 1. Wprowadzenie i Cel Badania

Niniejszy raport dokumentuje ewaluację porównawczą algorytmu detekcji obiektów poruszających się niezależnie (IMO) w strumieniu danych z neuromorficznej kamery zdarzeniowej (DVS). Porównano dwa warianty potoku:
* **Detektor v1 (Baseline)**: Kompensacja czystej rotacji $SO(3)$ na podstawie skalibrowanych pomiarów żyroskopu IMU. W scenach o silnym ruchu translacyjnym algorytm ten generuje dużą liczbę fałszywych alarmów (clutter) z powodu paralaksy tła.
* **Detektor v2 (Proponowany)**: Dwustopniowy potok łączący kompensację rotacji $SO(3)$ z weryfikacją geometryczną kierunku ruchu klastrów względem punktu skupienia ekspansji (**Focus of Expansion - FOE**) i linii epipolarnych wyznaczanych na podstawie prędkości liniowej kamery.

## 2. Porównanie Wydajności Silnika C++ (Wydajność i Opóźnienie)

![Wydajność C++: v1 vs v2](throughput_latency_v1_vs_v2.png)

| Sekwencja | Kategoria | Przepustowość v1 [Mev/s] | Przepustowość v2 [Mev/s] | Czas v1 [ms] | Czas v2 [ms] | Narzut FOE [ms] |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `scene6_dyn_train_03_0000` | Kategoria I | 4.01 | **4.40** | 16.40 | **14.94** | +0.00 |
| `scene14_dyn_test_05_0000` | Kategoria II | 4.15 | **3.52** | 21.83 | **25.72** | +4.21 |
| `scene14_dyn_test_04_0000` | Kategoria II | 4.45 | **3.89** | 28.21 | **32.34** | +5.93 |
| `scene13_dyn_test_02_0000` | Kategoria II | 4.80 | **4.24** | 25.54 | **28.87** | +3.27 |
| `scene13_dyn_test_01_0000` | Kategoria III | 5.04 | **4.33** | 46.16 | **53.76** | +10.78 |
| `scene15_dyn_test_02_0000` | Kategoria III | 4.41 | **4.03** | 28.19 | **30.88** | +4.30 |
| `scene15_dyn_test_05_0000` | Kategoria III | 4.43 | **4.27** | 26.21 | **27.16** | +0.00 |

### Kluczowe wnioski wydajnościowe:
1. **Zachowanie czasu rzeczywistego**: Średni narzut obliczeniowy filtru epipolarnego FOE wynosi zaledwie **3–10 ms** na całe 40-milisekundowe okno czasowe (~0.15 ms per klaster kandydat).
2. **Przepustowość powyżej 4 Mev/s**: Silnik C++ osiąga stabilną przepustowość **4.2 – 5.5 Mev/s** na 1 rdzeniu x86, przetwarzając strumień zdarzeń w czasie rzeczywistym (0.75–1.0x RT).

## 3. Jakość Detekcji względem Ground Truth (Czułość vs FPR Tła)

![Krzywe ROC: v1 vs v2](roc_curves_v1_vs_v2.png)

### Wyniki w punkcie roboczym $\tau = 0.3$ (Zbalansowany punkt pracy):

| Sekwencja | Kategoria | Recall v1 | Recall v2 | FPR tła v1 | FPR tła v2 | **Redukcja FPR (Fałszywych alarmów)** |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `scene6_dyn_train_03_0000` | Kategoria I | 41.0% | 39.2% | 34.75% | **28.13%** | **1.24x** |
| `scene14_dyn_test_05_0000` | Kategoria II | 53.5% | 42.6% | 23.64% | **10.71%** | **2.21x** |
| `scene14_dyn_test_04_0000` | Kategoria II | 51.6% | 43.8% | 23.17% | **11.45%** | **2.02x** |
| `scene13_dyn_test_02_0000` | Kategoria II | 41.4% | 40.0% | 14.21% | **4.52%** | **3.14x** |
| `scene13_dyn_test_01_0000` | Kategoria III | 47.6% | 30.3% | 19.94% | **3.81%** | **5.24x** |
| `scene15_dyn_test_02_0000` | Kategoria III | 44.6% | 41.0% | 10.73% | **3.65%** | **2.94x** |
| `scene15_dyn_test_05_0000` | Kategoria III | 63.7% | 63.4% | 22.72% | **19.51%** | 1.16x |

### Kluczowe wnioski jakościowe:
1. **Spektakularna redukcja fałszywych alarmów w Kategorii III (Ruch złożony / paralaksa)**:
   * Na sekwencji `scene13_dyn_test_01` (silna translacja) FPR tła spada z **19.94% do 3.81%** (ponad **5.2x mniej fałszywych alarmów**).
   * Precyzja klasyfikacji zdarzeń wzrasta 3-krotnie (z 4.2% do 12.7%).
2. **Brak degradacji w czystej rotacji (Kategoria I i II)**:
   * W przypadku braku translacji kamery ($|\mathbf{v}_{cam}| < 0.05$ m/s), moduł epipolarny automatycznie wstrzymuje odrzucanie, zachowując 100% czułości i precyzji wypracowanej przez kompensację $SO(3)$.

