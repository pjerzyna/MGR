# Raport wydajności: Silnik C++ vs Referencja Python (EVIMO2)

Porównanie czasu wykonania czystego potoku detekcji IMO v0 (bez GT) na reprezentatywnych sekwencjach EVIMO2.
Warunki pomiaru: okna 40 ms, krok 40 ms, 10 kolejnych okien w sekwencji, pojedynczy rdzeń procesora.

## 1. Tabela porównawcza czasów i przyspieszenia

| Sekwencja | Zdarzeń | Czas Python [ms] | Czas C++ [ms] | Śr. okno Python | Śr. okno C++ | **Przyspieszenie** | Throughput C++ | Czas rzecz. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `scene6_03` | 120,334 | 352.3 | 57.6 | 35.2 ms | 5.76 ms | **6.1×** | 2.09 Mev/s | 6.95× |
| `scene14_05` | 82,963 | 303.0 | 46.2 | 43.3 ms | 6.60 ms | **6.6×** | 1.80 Mev/s | 6.06× |
| `scene14_04` | 89,744 | 431.8 | 52.4 | 48.0 ms | 5.82 ms | **8.2×** | 1.71 Mev/s | 6.87× |
| `scene13_01` | 2,328,380 | 2512.8 | 248.7 | 251.3 ms | 24.87 ms | **10.1×** | 9.36 Mev/s | 1.61× |
| `scene13_05` | 2,216,701 | 2124.0 | 255.5 | 212.4 ms | 25.55 ms | **8.3×** | 8.67 Mev/s | 1.57× |
| `scene15_03` | 298,560 | 593.9 | 74.8 | 59.4 ms | 7.48 ms | **7.9×** | 3.99 Mev/s | 5.35× |

## 2. Rozbicie przyspieszenia na etapy potoku obliczeniowego

| Sekwencja | Undistort (Py → C++) | Warp SO(3) (Py → C++) | Detect / CCL (Py → C++) |
|---|---:|---:|---:|
| `scene6_03` | 5.4 ms → 0.50 ms (**10.8×**) | 6.9 ms → 0.61 ms (**11.3×**) | 22.9 ms → 4.65 ms (**4.9×**) |
| `scene14_05` | 16.9 ms → 0.40 ms (**42.1×**) | 6.4 ms → 0.75 ms (**8.5×**) | 20.1 ms → 5.45 ms (**3.7×**) |
| `scene14_04` | 19.9 ms → 0.46 ms (**43.2×**) | 5.4 ms → 0.52 ms (**10.4×**) | 22.7 ms → 4.83 ms (**4.7×**) |
| `scene13_01` | 94.7 ms → 5.36 ms (**17.7×**) | 131.4 ms → 13.09 ms (**10.0×**) | 25.1 ms → 6.43 ms (**3.9×**) |
| `scene13_05` | 56.4 ms → 5.66 ms (**10.0×**) | 127.1 ms → 13.39 ms (**9.5×**) | 28.9 ms → 6.51 ms (**4.4×**) |
| `scene15_03` | 25.1 ms → 0.62 ms (**40.5×**) | 15.9 ms → 1.76 ms (**9.0×**) | 18.4 ms → 5.10 ms (**3.6×**) |

## 3. Wnioski do pracy dyplomowej

- **Globalny zysk wydajnościowy:** Silnik C++ jest medianowo **8.1× szybszy** od referencyjnej implementacji w Pythonie (zakres 6.1×–10.1× w ujęciu całościowym potoku).
- **Zdolność pracy w czasie rzeczywistym:** Na wszystkich sekwencjach C++ wykonuje obliczenia szybciej niż wynosi długość okna 40 ms (średnio od 5.7 ms do 25.5 ms na okno w zależności od gęstości strumienia). Przepustowość sięga od **1.71 Mev/s** (dla sekwencji o małej liczbie zdarzeń) do aż **9.36 Mev/s** (dla sekwencji o skrajnej dynamice, np. `scene13_01`), na gęstych strumieniach wyraźnie przewyższając fizyczny limit przepustowości interfejsu czujnika Prophesee SilkyEvCam (~6.0 Mev/s). Cały potok działa ze współczynnikiem **1.57× do 6.95× czasu rzeczywistego** (na pojedynczym rdzeniu CPU).
- **Największy zysk na etapie Undistort i Warp:** Inwersja dystorsji (undistort) zyskała przyspieszenie aż **10.0×–43.2×**, a kompensacja rotacji (warp SO(3)) **8.5×–11.3×**. Wynika to z eliminacji narzutu alokacji tymczasowych tablic NumPy, braku interpretera Pythona oraz zwartego, przyjaznego dla pamięci podręcznej (cache-friendly) układu struktur w pamięci.
- **Metodologia i środowisko pomiaru:**
  - Kompilator C++: GCC 11+ (`-O3 -std=c++17 -Wall -Wextra`)
  - Środowisko Python: Python 3.10+, NumPy (C/BLAS backend), SciPy
  - Warm-up: 2 okna rozgrzewkowe przed aktywacją zegara pomiarowego
  - Zegar: `std::chrono::high_resolution_clock` (C++) oraz `time.perf_counter()` (Python)
  - Przetwarzanie sekwencyjne (1 wątek), bez wykorzystania instrukcji wielowątkowych OpenMP/threads (pozostawionych jako opcja rozwojowa).
