# Podsumowanie implementacji prototypu w języku Python (`evimo2`)

Dokument stanowi techniczny opis funkcjonalności, modeli matematycznych oraz algorytmów zaimplementowanych i zweryfikowanych w referencyjnym prototypie Pythona (pakiet `python/evimo2` oraz narzędzia badawcze `python/scripts`).

---

## 1. Architektura modułowa pakietu `evimo2`

Kod biblioteczny został podzielony na wyspecjalizowane moduły funkcyjne:

* **`geometry.py`**: Algebra kwaternionów, macierze rotacji $SO(3)$, wzór Rodriguesa, model kamery Pinhole oraz inwersja dystorsji RadTan.
* **`calibration.py`**: Model kalibracyjny żyroskopu IMU, kompensacja biasu, estymacja ekstrynsyk $R_{imu \to cam}$ (algorytm Kabscha/SVD) oraz estymacja opóźnienia czasowego $d$.
* **`rotcomp.py`**: Przygotowanie okien zdarzeń, dyskretne całkowanie prędkości kątowej, warping rotacyjny promieni $SO(3)$ oraz akumulacja biliniowa do obrazu IWE.
* **`detector.py`**: Dwustronna mapa niespójności czasowej (Early vs. Late), filtracja pudełkowa (SAT), progowanie gęstościowe, redukcja do komórek i klasteryzacja składowych spójnych (CCL).
* **`kinematics.py`**: Składanie transformacji obiektów w świecie ($T_{w \leftarrow o} = T_{w \leftarrow c} \cdot T_{c \leftarrow o}$), wyliczanie prędkości różnicami skończonymi, generowanie masek Ground Truth.
* **`common.py`**: I/O danych EVIMO2 (pliki `.npz`, `.npy`), indeksowanie zdarzeń, ekstrakcja strumieni czasowych.

---

## 2. Zastosowane koncepty matematyczne i geometryczne

### 2.1. Formuła Rodriguesa dla $SO(3)$ (`rodrigues_batch`)
Konwersja wektora obrotu (reprezentacja axis-angle $\mathbf{\theta} = \boldsymbol{\omega} \Delta t$) na macierz rotacji $3 \times 3$:
$$\mathbf{R}(\mathbf{\theta}) = \mathbf{I} + \frac{\sin \theta}{\theta} [\mathbf{\theta}]_\times + \frac{1 - \cos \theta}{\theta^2} [\mathbf{\theta}]_\times^2$$
gdzie $[\mathbf{\theta}]_\times$ to macierz antysymetryczna (ang. *skew-symmetric matrix*).
* **Numeryczna stabilność:** Dla kątów bliskich zeru ($\theta < 10^{-6}\text{ rad}$) zastosowano rozwinięcie w szereg Taylora:
  $$\frac{\sin \theta}{\theta} \approx 1 - \frac{\theta^2}{6}, \quad \frac{1 - \cos \theta}{\theta^2} \approx \frac{1}{2} - \frac{\theta^2}{24}$$
  zapobiegając dzieleniu przez zero i utracie precyzji zmiennoprzecinkowej.

### 2.2. Algebra kwaternionów Hamiltona (`geometry.py`)
* **Iloczyn Hamiltona ($\mathbf{a} \otimes \mathbf{b}$):** Wykorzystywany do sekwencyjnego składania obrotów bez zjawiska blokady kardanowej (*gimbal lock*).
* **Obrót wektora:** $\mathbf{v}' = \mathbf{q} \otimes (0, \mathbf{v}) \otimes \mathbf{q}^*$.
* **Unifikacja półsfery (`unify_hemisphere`):** Ponieważ kwaterniony $\mathbf{q}$ oraz $-\mathbf{q}$ reprezentują ten sam fizyczny obrót w $SO(3)$, w szeregach czasowych narzucono ciągłość znaku ($\mathbf{q}_i \cdot \mathbf{q}_{i+1} \ge 0$).
* **NLERP (Normalized Linear Interpolation):** Efektywna interpolacja orientacji w dowolnych chwilach wystąpienia zdarzeń $t_e$:
  $$\mathbf{q}(t) = \frac{(1 - \alpha) \mathbf{q}_i + \alpha \mathbf{q}_{i+1}}{\|(1 - \alpha) \mathbf{q}_i + \alpha \mathbf{q}_{i+1}\|}, \quad \alpha = \frac{t - t_i}{t_{i+1} - t_i}$$

### 2.3. Model kamery i inwersja dystorsji RadTan (`undistort`)
Model dystorsji radialno-tangencjalnej Browna-Conrady'ego:
$$x_d = x (1 + k_1 r^2 + k_2 r^4 + k_3 r^6) + 2 p_1 x y + p_2 (r^2 + 2x^2)$$
$$y_d = y (1 + k_1 r^2 + k_2 r^4 + k_3 r^6) + p_1 (r^2 + 2y^2) + 2 p_2 x y$$
* **Odwracanie metodą punktu stałego (Fixed-Point Iteration):**
  Numeryczne wyznaczanie znormalizowanego promienia $\mathbf{b} = (x, y, 1)^T$ z piksela $(u, v)$ przy 8 iteracjach punktu stałego:
  $$x^{(k+1)} = \frac{x_d - \Delta x_{\text{tang}}(x^{(k)}, y^{(k)})}{1 + k_1 (r^{(k)})^2 + k_2 (r^{(k)})^4 + k_3 (r^{(k)})^6}$$

---

## 3. Fuzja danych sensorycznych (DVS + IMU)

### 3.1. Równanie kalibracji żyroskopu (`calibrated_gyro`)
Przekształcenie surowego odczytu prędkości kątowej z żyroskopu IMU do układu współrzędnych kamery DVS:
$$\boldsymbol{\omega}_{cam}(t) = \mathbf{R}_{imu \to cam} (\boldsymbol{\omega}_{raw}(t + d) - \mathbf{b}_{gyro})$$
* **Ekstrynsyki:** $\mathbf{R}_{imu \to cam} \approx \text{diag}(-1, -1, 1)$ (odwrócony montaż układu scalonego IMU względem matrycy DVS na płytce PCB).
* **Bias żyroskopu:** $\mathbf{b}_{gyro} \approx (-0.0295, 0.0154, 0.0180)\text{ rad/s}$ (stały offset temperaturowy/sprzętowy).
* **Opóźnienie sprzętowe (Hardware Lag $d$):** Uwzględnienie asynchroniczności zegarów o przesunięcie $d \approx -15.7\text{ ms}$.

### 3.2. Dopasowanie Kabscha / SVD (`calibration.py`, `gyro_vs_gt.py`)
Analityczne wyznaczenie macierzy rotacji i wektora przesunięcia między trajektorią żyroskopu a pozycjami referencyjnymi Ground Truth z OptiTrack:
1. Centrowanie chmur wektorów: $\mathbf{X}_c = \mathbf{X} - \bar{\mathbf{X}}$, $\mathbf{Y}_c = \mathbf{Y} - \bar{\mathbf{Y}}$.
2. Rozkład według wartości osobliwych kowariancji: $\mathbf{U} \mathbf{\Sigma} \mathbf{V}^T = \text{SVD}(\mathbf{X}_c^T \mathbf{Y}_c)$.
3. Macierz rotacji z wymuszeniem $\det(\mathbf{R}) = +1$ (eliminacja niepożądanych odbić lustrzanych):
   $$\mathbf{R} = \mathbf{V} \operatorname{diag}\left(1, 1, \operatorname{sgn}(\det(\mathbf{V} \mathbf{U}^T))\right) \mathbf{U}^T$$

---

## 4. Kompensacja ruchu rotacyjnego $SO(3)$ (Motion Warping)

### 4.1. Dyskretyzacja i schodkowe całkowanie prędkości kątowej (`rot_matrices`)
Dla okna czasowego o szerokości $W = 40\text{ ms}$ wokół chwili centralnej $t_c$:
1. Zdefiniowano siatkę całkowania o stałym kroku $\Delta t = 1.0\text{ ms}$ ($nst = 40$ kroków).
2. Wyznaczono obroty węzłowe metodą Eulera-Rodriguesa:
   $$\mathbf{R}_g[k+1] = \mathbf{R}_g[k] \cdot \mathbf{R}(\boldsymbol{\omega}_k \Delta t)$$
3. **Interpolacja podkrokowa (Sub-bin Interpolation):** Dla każdego zdarzenia $e_i = (x_i, y_i, t_i)$ wyznaczany jest ułamek kroku $\delta t = t_i - g[k_i]$, dający pełną macierz orientacji w chwili zdarzenia:
   $$\mathbf{R}_{wc}(t_i) = \mathbf{R}_g[k_i] \cdot \mathbf{R}(\boldsymbol{\omega}_{k_i} \delta t)$$
4. **Względna transformacja rotacyjna:** Przekształcenie promienia do chwili referencyjnej środka okna:
   $$\mathbf{M}_i = \mathbf{R}_{wc}(t_c)^T \mathbf{R}_{wc}(t_i)$$

### 4.2. Warping współrzędnych i rzut perspektywiczny (`warp`)
Transformacja znormalizowanego wektora kierunkowego $\mathbf{b}_i$:
$$\mathbf{b}'_i = \mathbf{M}_i \mathbf{b}_i$$
$$x_s = f_x \frac{b'_{i,x}}{b'_{i,z}} + c_x, \quad y_s = f_y \frac{b'_{i,y}}{b'_{i,z}} + c_y$$

### 4.3. Akumulacja biliniowa do obrazu IWE (`iwe`)
Generowanie zrekonstruowanego obrazu skompensowanych zdarzeń (*Image of Warped Events*) z podziałem wagi zdarzenia na 4 sąsiednie piksele siatki (ang. *bilinear splatting*):
$$\Delta x = x_s - \lfloor x_s \rfloor, \quad \Delta y = y_s - \lfloor y_s \rfloor$$
$$w_{00} = (1 - \Delta x)(1 - \Delta y), \quad w_{10} = \Delta x (1 - \Delta y), \quad w_{01} = (1 - \Delta x)\Delta y, \quad w_{11} = \Delta x \Delta y$$
Zrealizowane wektorowo przez `np.bincount`, co eliminuje aliasing i utratę energii krawędzi.

---

## 5. Algorytm detekcji obiektów niezależnie poruszających się (IMO v0 Pipeline)

Algorytm opiera się na hipotezie, że po idealnej kompensacji rotacji $SO(3)$ statyczne tło pozostaje niezmienne w czasie, podczas gdy obiekty o niezależnym ruchu wywołują rozbieżność czasową między wczesną a późną połówką okna.

```
Strumień zdarzeń (okno W)
       │
       ▼
 [Warping SO(3)]
       │
       ├─────────────────────────┐
       ▼                         ▼
Połówka Early (t < tc)   Połówka Late (t >= tc)
       │                         │
  IWE: E (norm)             IWE: L (norm)
       │                         │
  BoxBlur 2x (SAT)          BoxBlur 2x (SAT)
       │                         │
       └───────────┬─────────────┘
                   ▼
       Mapa niespójności D = |E - L| / (E + L + eps)
       Mapa gęstości     S = E + L
                   │
                   ▼
       Progowanie podwójne: D >= tau  AND  S >= support * mean(S)
                   │
                   ▼
       Downsampling na siatkę bloków 4x4 px (próg aktywności >= 30%)
                   │
                   ▼
       Klasteryzacja spójnych składowych (CCL, 8-sąsiedztwo)
                   │
                   ▼
       Filtracja małych klastrów (< min_cells) + Fill Holes
                   │
                   ▼
       Maska logiczna zdarzeń IMO (True/False per zdarzenie)
```

### 5.1. Szybki filtr pudełkowy przez tablicę sum skumulowanych (`box_blur`)
Odpowiednik zoptymalizowanego filtru uśredniającego 2D w złożoności $O(1)$ na piksel:
1. Replikacja brzegów (`mode="edge"`).
2. Dwuosiowa tablica sum skumulowanych (SAT): `c[1:, 1:] = p.cumsum(0).cumsum(1)`.
3. Wyznaczenie sumy w oknie $(2r+1) \times (2r+1)$ operacją 4 odczytów.
Filtr aplikowany jest dwukrotnie ($2 \times$), co odpowiada przybliżeniu wygładzania filtrem trójkątnym/Gaussa.

### 5.2. Normalizowana mapa niespójności (`residual_map`)
$$D(x, y) = \frac{|E(x, y) - L(x, y)|}{E(x, y) + L(x, y) + \epsilon}, \quad S(x, y) = E(x, y) + L(x, y)$$
* $E, L$: znormalizowane rozkłady gęstości ($E / \sum E$, $L / \sum L$).
* $\epsilon = 10^{-12}$: składnik regularyzacyjny zabezpieczający przed dzieleniem przez zero w obszarach braku zdarzeń.

### 5.3. Siatka komórek i klasteryzacja (`label_cc`, `fill_holes`)
1. **Downsampling $4 \times 4$ px:** Zmniejsza szum pojedynczych pikseli. Blok komórki zostaje uznany za aktywny, gdy co najmniej 5 z 16 pikseli ($\ge 31.25\%$) spełnia warunki progowe.
2. **CCL (Connected Component Labeling):** Algorytm przeszukiwania wszerz (BFS) w 8-sąsiedztwie. Klastry o liczbie komórek mniejszej niż `min_cells` (domyślnie 8 komórek, tj. $128\text{ px}^2$) są odrzucane jako szum.
3. **Wypełnianie dziur (`fill_holes`):** Flood-fill w 4-sąsiedztwie od krawędzi kadru zapobiega powstawaniu pustych wnętrz wewnątrz zamkniętych obrysów poruszających się brył.

---

## 6. Kinematyka obiektów i weryfikacja Ground Truth

### 6.1. Złożenie transformacji obiektów w świecie (`kinematics.py`)
W zbiorze EVIMO2 trajektorie obiektów rejestrowane są w układzie współrzędnych kamery:
$$\mathbf{T}_{c \leftarrow o} = \begin{bmatrix} \mathbf{R}_{c \leftarrow o} & \mathbf{t}_{c \leftarrow o} \\ \mathbf{0} & 1 \end{bmatrix}$$
Aby wyznaczyć fizyczną prędkość obiektu w bezwładnościowym układzie świata, konieczne jest złożenie transformacji z pozą kamery $\mathbf{T}_{w \leftarrow c}$:
$$\mathbf{T}_{w \leftarrow o}(t) = \mathbf{T}_{w \leftarrow c}(t) \cdot \mathbf{T}_{c \leftarrow o}(t)$$
$$\mathbf{p}_{w, o}(t) = \mathbf{R}_{w \leftarrow c}(t) \mathbf{t}_{c \leftarrow o}(t) + \mathbf{t}_{w \leftarrow c}(t)$$
Prędkość liniowa obiektu wyznaczana jest centralną różnicą skończoną:
$$\mathbf{v}_o(t) = \frac{\mathbf{p}_{w, o}(t + \Delta t) - \mathbf{p}_{w, o}(t - \Delta t)}{2 \Delta t}, \quad v_o = \|\mathbf{v}_o\|_2$$

### 6.2. Kompensacja 3D Oracle z mapą głębi LiDAR (`warp_translation`)
Eksperyment referencyjny wyznaczający teoretyczną górną granicę (*upper bound*) kompensacji ego-ruchu:
$$\mathbf{P}_{ref} = \mathbf{M}_i \mathbf{P}_i + \mathbf{R}_{wc}(t_{ref})^T \left(\mathbf{p}_c(t_i) - \mathbf{p}_c(t_{ref})\right)$$
gdzie $\mathbf{P}_i = Z_i (b_{i,x}, b_{i,y}, 1)^T$, a $Z_i$ pochodzi bezpośrednio ze skojarzonej mapy głębi LiDAR/RGB-D (`dataset_depth.npz`).

---

## 7. Weryfikacja testami jednostkowymi i sanity-checkami

Poprawność prototypu została potwierdzona zestawem testów automatycznych (`pytest`):

| Test | Plik | Cel weryfikacji i warunek zaliczenia | Wynik |
|---|---|---|:---:|
| `test_package_exports` | `test_evimo2_pkg.py` | Kompletność API biblioteki `evimo2` | **PASSED** |
| `test_geometry_quaternions` | `test_evimo2_pkg.py` | Iloczyn Hamiltona, koniugacja, obrót wektora $[1,0,0]^T \to [0,1,0]^T$ przy $90^\circ$ | **PASSED** |
| `test_rodrigues_batch` | `test_evimo2_pkg.py` | Macierz jednostkowa dla obrotu $0^\circ$, analityczna macierz rotacji dla $90^\circ$ | **PASSED** |
| `test_undistort_pinhole` | `test_evimo2_pkg.py` | Odwracanie dystorsji dla $k_i=p_i=0$ (zgodność z rzutem perspektywicznym) | **PASSED** |
| `test_calibrated_gyro` | `test_evimo2_pkg.py` | Odejście biasu sprzętowego i przekształcenie osi $\mathbf{R}_{imu \to cam}$ | **PASSED** |
| `test_iwe_accumulation` | `test_evimo2_pkg.py` | Zachowanie masy biliniowej: $\sum IWE = N_{\text{events}}$ | **PASSED** |
| `test_detector_box_blur` | `test_evimo2_pkg.py` | Odpowiedź filtru pudełkowego na impuls Diraca: wartość centralna $1/(2r+1)^2$ | **PASSED** |
| `test_warp_translation` | `test_evimo2_pkg.py` | Zgodność analityczna rzutu 3D z translacją kamery i głębią $Z=2\text{ m}$ | **PASSED** |
| `test_sequence_format` | `test_sequence_format.py` | Walidacja integralności datasetu EVIMO2 (formaty `.npy`, zakresy czasowe, brak NaN) | **PASSED** |
