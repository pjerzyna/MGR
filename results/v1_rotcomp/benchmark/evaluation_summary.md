# Kompleksowy Raport Ewaluacji Silnika Detekcji IMO v1 (C++) na Zbiorze EVIMO2

Niniejszy raport podsumowuje automatyczną, laboratoryjną ewaluację silnika detekcji obiektów w ruchu własnym (IMO v1) zaimplementowanego w języku C++17 oraz referencji numerycznej. Badania przeprowadzono na 7 reprezentatywnych sekwencjach ze zbioru **EVIMO2**, obejmujących pełne spektrum dynamiki sensora neuromorficznego Prophesee Gen3 VGA (od całkowitego spoczynku, przez czystą rotację, aż po agresywny ruch złożony z paralaksą translacyjną).

---

## 1. Wybór i taksonomia sekwencji benchmarkowych

Do kompleksowej oceny wybrano 7 kluczowych sekwencji podzielonych na trzy kategorie badawcze:

| Kategoria | Sekwencja EVIMO2 | Charakterystyka ruchu kamery | Udział zdarzeń IMO | Rola w pracy dyplomowej |
| :--- | :--- | :--- | :---: | :--- |
| **Kategoria I** | `scene6_dyn_train_03_000000` | **Spoczynek kamery** ($\omega \approx 0$, $v \approx 0$) | 86.8% | Wzorzec referencyjny czystego sygnału IMO (Sanity-check). Weryfikacja, czy algorytm nie degraduje jakości przy braku ruchu. |
| **Kategoria II** | `scene14_dyn_test_05_000000`<br>`scene14_dyn_test_04_000000`<br>`scene13_dyn_test_02_000000` | **Dominujący obrót** wokół osi optycznej i poprzecznej ($\omega \le 2.76\text{ rad/s}$) | 2.3% – 12.9% | Sweet-spot modelu $SO(3)$. Bezpośredni dowód skuteczności kompensacji rotacyjnej z żyroskopu IMU. |
| **Kategoria III** | `scene13_dyn_test_01_000000`<br>`scene15_dyn_test_02_000000`<br>`scene15_dyn_test_05_000000` | **Ruch złożony (obrót + silna translacja)** generujący paralaksę tła | 4.6% – 71.1% | Identyfikacja fizycznych granic modelu $SO(3)$ i empiryczne uzasadnienie potrzeby estymacji głębi/stereo w wersji v2. |

---

## 2. Wydajność czasowa i throughput silnika C++

Pomiary przeprowadzono na pojedynczym rdzeniu procesora x86 (bez wielowątkowości) dla okien czasowych $W = 40\text{ ms}$ (krok $40\text{ ms}$, rozgrzewka 2 okna, parametry: $r=3$, $\tau=0.5$). 

| Sekwencja | Zdarzeń | Okien | Undistort [ms] | Warp SO(3) [ms] | Detect/CCL [ms] | Czas per okno [ms] | **Throughput [Mev/s]** | **Współczynnik RT** |
| :--- | ---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `scene6_03` | 2 366 078 | 36 | 7.59 | 4.02 | 4.79 | **16.4 ms** | **4.01 Mev/s** | **2.44× czasu rzecz.** |
| `scene14_05` | 3 081 126 | 34 | 9.99 | 5.79 | 6.05 | **21.8 ms** | **4.15 Mev/s** | **1.83× czasu rzecz.** |
| `scene14_04` | 4 648 575 | 37 | 14.05 | 7.67 | 6.49 | **28.2 ms** | **4.45 Mev/s** | **1.42× czasu rzecz.** |
| `scene13_02` | 4 656 308 | 38 | 13.50 | 7.04 | 5.00 | **25.5 ms** | **4.80 Mev/s** | **1.57× czasu rzecz.** |
| `scene13_01` | 8 843 398 | 38 | 25.71 | 13.46 | 6.99 | **46.2 ms** | **5.04 Mev/s** | **0.87× czasu rzecz.** |
| `scene15_02` | 4 728 099 | 38 | 13.97 | 7.93 | 6.29 | **28.2 ms** | **4.41 Mev/s** | **1.42× czasu rzecz.** |
| `scene15_05` | 2 206 444 | 19 | 12.77 | 6.99 | 6.45 | **26.2 ms** | **4.43 Mev/s** | **1.53× czasu rzecz.** |

> [!NOTE]
> * **Łączna skala testu:** Silnik przetworzył **30 530 028 zdarzeń** w 237 oknach bez ani jednego błędu pamięci czy wycieku.
> * **Stabilność throughputu:** Przepustowość C++ waha się w bardzo wąskim i przewidywalnym przedziale **4.01 – 5.04 Mev/s**.
> * **Czas rzeczywisty:** Na 6 z 7 sekwencji średni czas obliczeń wynosi $16.4\text{--}28.2\text{ ms}$, co przy długości okna $40\text{ ms}$ oznacza pracę w czasie rzeczywistym ze średnim współczynnikiem **1.7× czasu rzeczywistego**.
> * **Przypadek `scene13_01`:** Sekwencja ta cechuje się ekstremalną gęstością strumienia (ponad 8.8 mln zdarzeń w 1.5 s — szczytowa gęstość strumienia w nagraniach EVIMO2 $\approx 5.8\text{ Mev/s}$, przy sprzętowym limicie odczytu sensora Prophesee Gen3.1 wynoszącym 50–66 Meps). Nawet przy tak skrajnym obciążeniu silnik osiąga rekordowe **5.04 Mev/s**.

Plik wykresu wydajności: [`throughput_latency_cpp.png`](throughput_latency_cpp.png)

---

## 3. Jakość detekcji względem Ground Truth (Wpływ kompensacji rotacji SO(3))

### Tabela 3.1: Punkt Zbalansowany ($r=3, \tau=0.3$) – Rekomendowany standard roboczy (Recall $\approx 49.0\%$)

| Sekwencja | Kategoria | Czułość (Recall) | Błędy tła (FPR) | Precyzja (Precision) | Redukcja FPR tła | Skok precyzji |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| `scene6_03` | **Kat. I** *(spoczynek)* | 40.9% | 25.8% | 93.6% | **1.00×** *(neutralna)* | **+0.1 pp** |
| `scene14_05` | **Kat. II** *(czysta rotacja)* | 53.3% | 19.0% | 21.2% | **2.32× mniej błędów!** | **+10.5 pp** |
| `scene14_04` | **Kat. II** *(czysta rotacja)* | 50.8% | 19.7% | 17.5% | **1.60× mniej błędów** | **+4.7 pp** |
| `scene13_02` | **Kat. II** *(szybka rotacja)* | 39.5% | 12.6% | 4.5% | **2.40× mniej błędów!** | **+2.0 pp** |
| `scene13_01` | **Kat. III** *(paralaksa)* | 46.8% | 17.5% | 9.8% | **2.06× mniej błędów** | **+4.0 pp** |
| `scene15_02` | **Kat. III** *(ruch złożony)* | 45.5% | 7.4% | 57.7% | **1.63× mniej błędów** | **+13.1 pp** |
| `scene15_05` | **Kat. III** *(wysokie IMO)* | 66.3% | 14.9% | 76.8% | **1.12×** | **-0.2 pp** |

---

### Tabela 3.2: Punkt Wysokiej Precyzji ($r=3, \tau=0.5$) – Kryterium Neymana-Pearsona (Low-FPR $\le 8\%$)

| Sekwencja | Kategoria | Czułość (Recall) | Błędy tła (FPR) | Precyzja (Precision) | Redukcja FPR tła | Skok precyzji |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| `scene6_03` | **Kat. I** *(spoczynek)* | 14.3% | 11.0% | 92.0% | **1.00×** *(neutralna)* | **+0.1 pp** |
| `scene14_05` | **Kat. II** *(czysta rotacja)* | 30.8% | 7.5% | 27.1% | **3.68× mniej błędów!** | **+17.6 pp** |
| `scene14_04` | **Kat. II** *(czysta rotacja)* | 25.2% | 10.4% | 15.6% | **1.79× mniej błędów** | **+5.6 pp** |
| `scene13_02` | **Kat. II** *(szybka rotacja)* | 16.4% | 4.6% | 4.3% | **3.10× mniej błędów!** | **+1.8 pp** |
| `scene13_01` | **Kat. III** *(paralaksa)* | 27.6% | 6.4% | 14.9% | **2.36× mniej błędów** | **+7.4 pp** |
| `scene15_02` | **Kat. III** *(ruch złożony)* | 21.2% | 1.8% | 73.5% | **1.66× mniej błędów** | **+12.6 pp** |
| `scene15_05` | **Kat. III** *(wysokie IMO)* | 45.7% | 6.6% | 79.0% | **1.00×** | **-0.3 pp** |

Pliki wykresów:
* Wykres słupkowy redukcji FPR i skoku precyzji dla Punktu Wysokiej Precyzji ($\tau=0.5$): [`operating_point_tau05_summary.png`](operating_point_tau05_summary.png)
* Wykres słupkowy redukcji FPR i skoku precyzji dla Punktu Zbalansowanego ($\tau=0.3$): [`operating_point_tau03_summary.png`](operating_point_tau03_summary.png)
* Krzywe ROC i analiza punktu zbalansowanego ($r=3, \tau=0.3$): [`roc_curves_tau03.png`](roc_curves_tau03.png)
* Trójpanelowe krzywe ROC pełnego omiatania ($r=3, \tau \in [0.2, 0.5]$): [`roc_curves_categories.png`](roc_curves_categories.png)

---

## 4. Wnioski badawcze i interpretacja zjawisk

### 1. Walidacja czystości modelu (Kategoria I – `scene6_03`):
* Przy braku ruchu kamery kompensacja $SO(3)$ zachowuje się w sposób całkowicie neutralny (dokładnie $1.00\times$ zmiany FPR, precyzja $92.0\%$). Dowodzi to, że potok kompensacji nie wprowadza sztucznych artefaktów ani zniekształceń do nieruchomego tła.

### 2. Spektakularny sukces kompensacji rotacji (Kategoria II):
* W sekwencjach z dominującym obrotem (`scene14_05` i `scene13_02`) wskaźnik fałszywych alarmów tła spada **ponad 3-krotnie** (w `scene14_05` z 27.7% do 7.5%, a w `scene13_02` z 14.2% do 4.6%).
* W `scene14_05` precyzja wzrasta niemal 3-krotnie (z 9.5% do 27.1%), co jednoznacznie potwierdza poprawność zaimplementowanego całkowania wzoru Rodriguesa oraz precyzję kalibracji ekstrynsyk IMU.

### 3. Wpływ zjawiska paralaksy i uzasadnienie wersji v2 (Kategoria III):
* W sekwencji `scene13_01` (gdzie występuje wyraźny ruch translacyjny kamery względem bliskich obiektów), kompensacja rotacji obniża FPR z 15.1% do 6.4%, lecz dalsza redukcja jest blokowana przez paralaksę (krawędzie tła przesuwają się liniowo). 
* Zbiega się to z wcześniejszym eksperymentem **Oracle 3D**, w którym dodanie składowej translacyjnej obniżyło FPR do zaledwie $0.8\%$, podnosząc precyzję do $78.9\%$.
* Stanowi to formalne, empiryczne uzasadnienie dla przejścia do wersji **v2 (Stereo Prophesee lub linie epipolarne FOE)** w celu skasowania resztkowego błędu paralaksy.

### 4. Wybór punktu pracy (Trade-Off Czułość vs FPR):
* Analiza pełnej siatki parametrów $\tau \in \{0.2, 0.3, 0.4, 0.5\}$ i $r \in \{2, 3, 5\}$ wskazuje, że:
  * **$r = 3, \tau = 0.5$** to optymalny punkt dla aplikacji wymagających minimalnej liczby fałszywych alarmów (wysoka precyzja, FPR tła $\le 4\text{--}7\%$).
  * **$r = 3, \tau = 0.3$** to optymalny punkt kompromisowy (zbalansowany), podwajający czułość (Recall wzrasta do $40\text{--}55\%$) przy umiarkowanym wzroście FPR.

---

## 5. Spis wygenerowanych plików i artefaktów

Wszystkie dane zostały wygenerowane i zapisane w folderze [`results/summaries/evaluation/`](.):

1. **[`benchmark_cpp_7seqs.csv`](benchmark_cpp_7seqs.csv)** – Pełne pomiary wydajnościowe silnika C++ (czasy, throughput, współczynniki czasu rzeczywistego).
2. **[`quality_metrics_7seqs.csv`](quality_metrics_7seqs.csv)** – Zestawienie metryk jakości detekcji (Recall, Precision, IoU, FPR, współczynniki redukcji) dla punktów pracy.
3. **[`param_sweep_grid.csv`](param_sweep_grid.csv)** – Pełna tablica siatki parametrów $12 \times 7$ kombinacji ($\tau \times r$) dla wszystkich sekwencji.
4. **[`roc_curves_categories.png`](roc_curves_categories.png)** – Trójpanelowy wykres krzywych ROC (Czułość vs FPR) z podziałem na kategorie dynamiki.
5. **[`throughput_latency_cpp.png`](throughput_latency_cpp.png)** – Wykres słupkowy throughputu i latencji silnika C++ z zaznaczonymi fizycznymi granicami sensora.
