# Testy numeryczne C++ (EVIMO2)

Ten katalog zawiera testy numeryczne weryfikujące zgodność implementacji silnika C++ (`cpp/include/`) względem referencyjnego kodu Python (`python/evimo2/`).

## 1. Źródło danych testowych (`cpp/tests/data/`)

Pliki w katalogu `data/` pochodzą z rzeczywistego nagrania z datasetu **EVIMO2**:
* **Sekwencja:** `scene14_dyn_test_05_000000` (skrót: `scene14_05`)
* **Klatka odniesienia:** `frame_id = 85`
* **Środek okna czasowego:** $t_c = 1.433333\text{ s}$
* **Szerokość okna:** $W = 40.0\text{ ms}$ ($[1.413333; 1.453333]\text{ s}$)
* **Liczba zdarzeń w oknie:** 233 144
* **Liczba próbek IMU:** 205 (przedział $\pm 100\text{ ms}$ wokół okna)
* **Kalibracja:** $R_{IMU \to CAM} = \text{diag}(-1, -1, 1)$, $\mathbf{b} = (-0.0295, 0.0154, 0.0180)\text{ rad/s}$, opóźnienie $d = -15.7\text{ ms}$

Dane te wygenerowano automatycznie za pomocą skryptu:
```bash
python python/scripts/export_cpp_reference.py
```

## 2. Zawartość plików referencyjnych

| Plik | Typ / Rozmiar | Opis |
|---|---|---|
| `t_ev.npy` | `float64 (233144,)` | Znaczniki czasu zdarzeń w sekundach |
| `pix_ev.npy` | `float64 (233144, 2)` | Współrzędne pikseli $(x, y)$ |
| `p_ev.npy` | `uint8 (233144,)` | Polaryzacja zdarzeń (0 / 1) |
| `imu_ts.npy` | `float64 (205,)` | Znaczniki czasu próbek IMU |
| `imu_raw_gyro.npy`| `float64 (205, 3)` | Surowe prędkości kątowe z żyroskopu [rad/s] |
| `imu_w.npy` | `float64 (205, 3)` | Skalibrowane prędkości kątowe w układzie kamery |
| `bearing.npy` | `float64 (233144, 3)` | Promienie jednostkowe $(b_x, b_y, 1.0)$ po inwersji dystorsji RadTan |
| `grid.npy` | `float64 (41,)` | Siatka dyskretyzacji całkowania $1\text{ ms}$ |
| `rot_matrices.npy`| `float64 (233144, 3, 3)`| Macierze obrotu $M_i = R(t_{ref})^T R(t_i)$ z całkowania kątowego |
| `warped_xs.npy`, `warped_ys.npy` | `float64 (233144,)` | Pozycje pikseli po rzutowaniu promieni z kompensacją rotacji |
| `iwe_E.npy`, `iwe_L.npy` | `float32 (480, 640)` | Obrazy IWE (Image of Warped Events) dla wczesnej i późnej połówki |
| `D.npy`, `S.npy` | `float64 (480, 640)` | Mapa niespójności $D$ i mapa wsparcia $S$ |
| `pred_imo.npy` | `bool (233144,)` | Ostateczna maska klasyfikacji IMO per zdarzenie wyznaczona przez Pythona |
| `meta.json` | JSON | Metadane okna, intrinsics kamery oraz parametry detektora |

## 3. Uruchamianie testów

```bash
cd cpp
make test
```
Test weryfikuje zgodność każdego etapu obliczeń z dokładnością do $10^{-6}$ (dla algebry $SO(3)$ błąd wynosi poniżej $10^{-15}$).
