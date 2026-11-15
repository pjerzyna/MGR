# Metodologia Doboru Punktów Pracy (r, tau) Detektora IMO

Niniejszy dokument formalizuje matematyczny i empiryczny proces doboru optymalnych hiperparametrów detektora IMO:
* **$r$** – promień 2-krotnego filtru pudełkowego (Box Blur przez tablicę sum skumulowanych SAT), określający rozmiar jądra przestrzennego $(2r+1) \times (2r+1)\text{ px}$,
* **$\tau$** – próg niespójności znormalizowanej mapy różnicowej $D = \frac{|E - L|}{E + L + \epsilon} \in [0, 1]$.

Wykres analityczny: [`operating_points_analysis.png`](operating_points_analysis.png)  
Surowe dane zagregowane: [`operating_points_aggregate.csv`](operating_points_aggregate.csv)  
Pełna siatka pomiarowa (84 wiersze): [`param_sweep_grid.csv`](param_sweep_grid.csv)

---

## 1. Algorytm i kroki wyznaczania punktów pracy

W problemie detekcji obiektów w ruchu własnym (IMO) mamy do czynienia z **optymalizacją wielokryterialną** o sprzecznych celach:
1. **Maksymalizacja czułości (Recall):** chcemy wykryć jak największą część zdarzeń generowanych przez poruszający się obiekt ($TP / P \to 1$).
2. **Minimalizacja fałszywych alarmów (FPR tła):** chcemy uniknąć błędnego oznaczania skompensowanego tła jako obiekt ($FP / N \to 0$).

Wybór hiperparametrów zrealizowano w 4 krokach:

### Krok 1: Dyskretyzacja przestrzeni parametrów (Grid Search)
Zbadano siatkę $3 \times 4 = 12$ konfiguracji parametrów na wszystkich 7 sekwencjach referencyjnych EVIMO2:
$$r \in \{2, 3, 5\} \quad (\text{jądra } 5\times 5, 7\times 7, 11\times 11\text{ px}), \qquad \tau \in \{0.2, 0.3, 0.4, 0.5\}$$

### Krok 2: Konstrukcja Frontu Pareto w przestrzeni ROC
Dla każdej wartości promienia $r$ wyznaczono krzywą kompromisu $(\text{FPR}(\tau), \text{Recall}(\tau))$. 
Punkty niebędące zdominowane (czyli takie, dla których nie istnieje inny punkt o jednocześnie wyższym Recall i niższym FPR) tworzą **Front Pareto**.

### Krok 3: Wybór promienia wygładzania $r$ (Fizyka sensora DVS)
* **$r = 2$ ($5 \times 5\text{ px}$): ZBYT MAŁE ROZMYCIE.** Dyskretne zdarzenia tła i mikroskopijny jitter sensora nie nakładają się w pełni między połówką wczesną $E$ i późną $L$. Efekt: wysoki szum resztkowy i drastyczny wzrost błędów tła ($\text{FPR} = 11.2\%\text{--}38.7\%$).
* **$r = 5$ ($11 \times 11\text{ px}$): ZBYT DUŻE ROZMYCIE.** Filtr wygładza nie tylko szum tła, ale także drobne krawędzie małych obiektów. Efekt: sygnał obiektu ulega stłumieniu, a czułość gwałtownie spada ($\text{Recall}$ maleje poniżej $17\%\text{--}36\%$).
* **$r = 3$ ($7 \times 7\text{ px}$): PUNKTY OPTYMALNE.** Zapewnia optymalną równowagę — skutecznie skleja krawędzie tła o grubości 1–2 px, nie zacierając konturów poruszających się obiektów.

### Krok 4: Wybór progu niespójności $\tau$ według kryteriów decyzyjnych

Zależnie od przeznaczenia algorytmu w systemie robotycznym zdefiniowano dwa rozłączne kryteria:

#### A. Punkt Zbalansowany (Maksymalizacja wskaźnika Youdena $J$):
Wskaźnik Youdena $J$ definiuje punkt na krzywej ROC najbardziej oddalony od losowej klasyfikacji:
$$J(r, \tau) = \text{Recall}(r, \tau) - \text{FPR}(r, \tau)$$
* **Wynik:** Maksimum wskaźnika Youdena $J = 0.323\text{--}0.374$ oraz wskaźnika $F_1 = 0.44\text{--}0.48$ przypada na **$r = 3, \tau = 0.3$**.
* **Charakterystyka:** $\text{Recall} \approx 49.0\%$, $\text{FPR} \approx 16.7\%$. Punkt ten maksymalizuje liczbę poprawnie wykrytych zdarzeń obiektu, akceptując umiarkowany poziom szumu.

#### B. Punkt Wysokiej Precyzji (Ograniczenie Neymana-Pearsona na błędy tła):
Dla autonomicznego robota fałszywa detekcja przeszkody (ghost detection) może spowodować niepotrzebne gwałtowne hamowanie. Wprowadzamy twarde ograniczenie inżynierskie:
$$\max \text{Recall} \quad \text{przy warunku} \quad \text{FPR} \le 0.08 \quad (\le 8\%)$$
* **Wynik:** Jedynym punktem na froncie Pareto spełniającym ten warunek przy zachowaniu sensownej czułości jest **$r = 3, \tau = 0.5$**.
* **Charakterystyka:** $\text{FPR} = 6.90\%$ (w sekwencjach czystej rotacji $\text{FPR} \le 4.6\%\text{--}7.5\%$, a w `scene15_02` $\text{FPR} = 1.76\%$), $\text{Precision} \approx 43.8\%$, $\text{Recall} \approx 25.9\%$.

---

## 2. Zagregowana tabela decyzyjna (Średnie ze wszystkich 7 sekwencji)

| Promień $r$ | Próg $\tau$ | Czułość (Recall) | Błędy tła (FPR) | Precyzja (Precision) | Jakość (IoU) | **Wskaźnik Youdena $J$** | **$F_1$-score** | Ocena inżynierska |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| 2 | 0.2 | 74.7% | 38.7% | 35.2% | 30.1% | 0.360 | 0.478 | Zbyt zaszumione tło |
| 2 | 0.3 | 58.8% | 24.8% | 37.0% | 26.7% | 0.339 | 0.454 | Zbyt wysoki FPR |
| 2 | 0.4 | 45.0% | 16.5% | 38.5% | 22.9% | 0.286 | 0.415 | Gorszy niż $r=3$ |
| 2 | 0.5 | 34.1% | 11.2% | 39.9% | 19.1% | 0.228 | 0.368 | Gorszy niż $r=3$ |
| **3** | **0.2** | 65.6% | 28.2% | 37.5% | 29.0% | **0.374** | 0.477 | Bardzo wysoki Recall, ale spory szum |
| **3** | **0.3** | **49.0%** | **16.7%** | **40.2%** | **25.0%** | **0.323** | **0.441** | 🟢 **PUNKT ZBALANSOWANY (Standard)** |
| 3 | 0.4 | 35.5% | 10.6% | 42.0% | 20.6% | 0.249 | 0.385 | Przejściowy |
| **3** | **0.5** | **25.9%** | **6.9%** | **43.8%** | **16.6%** | **0.190** | **0.325** | 🟢 **PUNKT WYSOKIEJ PRECYZJI (Low-FPR)** |
| 5 | 0.2 | 51.8% | 17.6% | 41.5% | 26.9% | 0.342 | 0.461 | Zbyt duży narzut obliczeniowy |
| 5 | 0.3 | 36.0% | 9.6% | 45.4% | 22.3% | 0.263 | 0.401 | Tłumienie małych obiektów |
| 5 | 0.4 | 25.1% | 5.8% | 47.7% | 17.8% | 0.194 | 0.329 | Silne tłumienie krawędzi |
| 5 | 0.5 | 17.2% | 3.5% | 49.6% | 13.7% | 0.137 | 0.256 | Obiekty niemal niewykrywalne |
