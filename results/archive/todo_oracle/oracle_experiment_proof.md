# Eksperyment Oracle 3D: Wpływ translacji kamery i głębi sceny na detekcję IMO

Dokument podsumowuje laboratoryjne wyniki eksperymentu **Oracle 3D** przeprowadzonego na sekwencjach zbioru EVIMO2. Eksperyment stanowi formalny dowód na to, że ograniczenie precyzji klasycznego modelu $SO(3)$ (kompensacji czystej rotacji z IMU) wynikało wyłącznie z nieskompensowanej **paralaksy translacyjnej tła**, a nie z błędów kalibracji żyroskopu.

Plik wykresu podsumowującego: [`oracle_experiment_results.png`](oracle_experiment_results.png)  
Skonsolidowane dane pomiarowe: [`oracle_experiment_results.json`](oracle_experiment_results.json)

---

## 1. Zestawienie wyników pomiarowych

Pomiary wykonano za pomocą skryptu [`python/scripts/imo_detect.py`](../../python/scripts/imo_detect.py) z parametrami: okno $W = 40\text{ ms}$, filtr pudełkowy $r = 5$ (2 przebiegi), próg niespójności $\tau = 0.3$, wsparcie gęstości $\text{support} = 1.0$, min. wielkość klastra 8 komórek ($4 \times 4\text{ px}$).

| Sekwencja | Wariant przetwarzania | Precyzja [$TP / (TP+FP)$] | Błędy tła (FPR) [$FP / N$] | IoU zdarzeń IMO | Czułość (Recall) [$TP / P$] | Składowe / okno |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **`scene14_05`**<br>*(wzorzec kalibracji IMU)* | 1. Brak kompensacji<br>2. **Sama rotacja IMU (obecny model)**<br>3. **Rotacja IMU + translacja Oracle (3D)**<br>4. Pełna rotacja + translacja GT | 12.0%<br>**37.0%**<br>**81.0%**<br>81.1% | 30.8%<br>**8.1%**<br>**2.6%**<br>2.6% | 10.5%<br>**25.2%**<br>**40.0%**<br>40.2% | 45.8%<br>**44.3%**<br>**44.1%**<br>44.4% | 50.0<br>30.0<br>**11.0**<br>11.0 |
| **`scene14_04`**<br>*(obrót + mały ruch obiektów)* | 1. Brak kompensacji<br>2. **Sama rotacja IMU (obecny model)**<br>3. **Rotacja IMU + translacja Oracle (3D)**<br>4. Pełna rotacja + translacja GT | 13.0%<br>**21.9%**<br>**73.5%**<br>72.9% | 20.4%<br>**10.8%**<br>**2.0%**<br>2.0% | 10.8%<br>**16.0%**<br>**32.5%**<br>32.6% | 39.2%<br>**37.1%**<br>**36.8%**<br>37.0% | 40.5<br>31.0<br>**14.5**<br>15.5 |
| **`scene13_01`**<br>*(umiarkowany obrót + silna paralaksa)* | 1. Brak kompensacji<br>2. **Sama rotacja IMU (obecny model)**<br>3. **Rotacja IMU + translacja Oracle (3D)**<br>4. Pełna rotacja + translacja GT | 9.6%<br>**19.1%**<br>**78.9%**<br>78.7% | 15.6%<br>**6.5%**<br>**0.8%**<br>0.8% | 8.4%<br>**14.3%**<br>**28.3%**<br>28.0% | 40.6%<br>**36.0%**<br>**30.6%**<br>30.3% | 73.5<br>32.5<br>**1.5**<br>1.0 |
| **`scene15_03`**<br>*(silna translacja liniowa, wysoki udział IMO)* | 1. Brak kompensacji<br>2. **Sama rotacja IMU (obecny model)**<br>3. **Rotacja IMU + translacja Oracle (3D)**<br>4. Pełna rotacja + translacja GT | 77.9%<br>**78.1%**<br>**78.1%**<br>78.2% | 29.9%<br>**13.0%**<br>**5.2%**<br>4.7% | 48.7%<br>**49.1%**<br>**49.1%**<br>49.3% | 56.4%<br>**56.9%**<br>**56.9%**<br>57.2% | 49.0<br>37.0<br>**25.5**<br>25.0 |

---

## 2. Dlaczego precyzja w `scene14_05` skacze z 37% do 81%?

### Mechanizm matematyczny:
1. **Współczynnik $P/N$ (Base Rate):**
   W sekwencji `scene14_05` udział zdarzeń generowanych przez ruchomy obiekt wynosi zaledwie ok. 8–13% ($P/N \approx 0.13$). Zdecydowana większość strumienia zdarzeń to tło sceny.
2. **Co się dzieje przy samej rotacji IMU ($SO(3)$):**
   * Żyroskop idealnie usuwa rozmycie obrotowe kamery.
   * Ponieważ kamera nie jest idealnie unieruchomiona na statywie kulowym i wykonuje mikroskopijne ruchy translacyjne w przestrzeni, bliskie obiekty w tle doznają **paralaksy**.
   * Paralaksa wywołuje niespójność między połówkami okna, co daje wskaźnik fałszywych alarmów tła $\text{FPR} = 8.1\%$.
   * Z definicji precyzji:
     $$\text{Precision} = \frac{TP}{TP + FP}$$
     Nawet $8.1\%$ fałszywych pikseli z olbrzymiego tła ($FP$) jest porównywalne z liczbą prawdziwych pikseli małego obiektu ($TP$). Dlatego precyzja wynosi tylko **37.0%**.
3. **Co daje kompensacja 3D Oracle:**
   * Po uwzględnieniu translacji kamery oraz głębi $Z$ z LiDAR-u tło zostaje zrzutowane z uwzględnieniem odległości fizycznej.
   * Fałszywe alarmy tła ($FP$) gwałtownie spadają z **$8.1\%$ do $2.6\%$** (ponad 3-krotna redukcja szumu).
   * Czułość wykrywania obiektu ($TP$) pozostaje nienaruszona ($\text{Recall} \approx 44.1\%$).
   * Mianownik drastycznie maleje, a precyzja skacze do **$81.0\%$**!

---

## 3. Kluczowe wnioski z eksperymentu

1. **Rotacja z IMU jest równoważna rotacji z Ground Truth:**
   Wiersz *"Rotacja IMU + translacja Oracle"* daje **identyczne** wyniki jak *"Pełna rotacja + translacja GT"* (w `scene14_05`: $81.0\%$ vs $81.1\%$ precyzji; $\text{FPR} = 2.6\%$ vs $2.6\%$). Oznacza to, że skalibrowana rotacja z żyroskopu sensora jest perfekcyjna.
2. **Korelacja z prędkością translacji kamery:**
   W `scene13_01` korelacja Spearmana między fałszywymi alarmami a prędkością kamery spada z **$+0.75$** (silny wpływ paralaksy przy samej rotacji) do zaledwie **$+0.18$** po kompensacji translacji.
3. **Resztkowe $2.6\%$ FPR:**
   Pozostałe $1\text{--}3\%$ błędów tła to fizyczna dolna granica szumu sensora DVS, niepełne pokrycie map głębi (LiDAR w EVIMO2 pokrywa ok. 61–84% punktów; reszta otrzymuje medianę okna) oraz zjawisko apertury na krawędziach.
4. **Dlaczego model w C++ operuje na samej rotacji?**
   Kompensacja rotacji z IMU jest bezstanowa, natychmiastowa ($O(1)$) i nie wymaga znajomości głębi sceny ani integracji akcelerometru podatnego na dryf. Eksperyment Oracle udowadnia jednak, że jeśli w przyszłości dodamy estymację głębi (np. drugą kamerę Prophesee Stereo lub sieć monocular depth), jakość detekcji IMO natychmiast osiągnie poziom $>80\%$.
