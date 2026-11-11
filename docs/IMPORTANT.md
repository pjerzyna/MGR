## Status pytań badawczych i ustaleń [PODSUMOWANIE]
- Zakresy czasu zdarzeń, IMU i klatek GT pokrywają się (start ≈ 0,017 s); parametr `ros_time_offset` nie jest potrzebny, ale znaczniki IMU wykazują opóźnienie względem GT rzędu −31…+14 ms.
- Konwencja poz kamery: ROZSTRZYGNIĘTE → Konwencja A ($T_{świat \leftarrow kamera}$, $t$ = pozycja kamery w świecie; potwierdzona w 19/19 sekwencji).
- Ekstrynsyka IMU→kamera: ROZSTRZYGNIĘTE → $R_{imu \to kamera} \approx \text{diag}(-1, -1, 1)$ (obrót o 180° wokół osi z, odchylenia $< 0,8^\circ$).
- Bias żyroskopu: ROZSTRZYGNIĘTE → stały bias sprzętowy $\mathbf{b} \approx (-0,0295; 0,0154; 0,0180)$ rad/s ($|\mathbf{b}| \approx 0,038$ rad/s, stabilność $\pm 0,001$ rad/s).
- Konwencja poz obiektów: ROZSTRZYGNIĘTE → $T_{kamera \leftarrow obiekt}$ (klasa $R*S$, potwierdzona w 11/11 sekwencji z obrotem $\ge 10^\circ$). Do obliczenia prędkości w świecie konieczne jest złożenie $T_{w \leftarrow o} = T_{w \leftarrow c} \cdot T_{c \leftarrow o}$.
- Jednostki głębi w `dataset_depth.npz`: ROZSTRZYGNIĘTE → milimetry (pozy w metrach; stosunek mediany głębi pod maską do współrzędnej $z$ z pozy wynosi 974–989). Typowe głębie obiektów: 0,5–1,3 m.
- Tablica `index` w `dataset_info.npz` (654 elementy): ROZSTRZYGNIĘTE → prekomputowany indeks czasowy zdarzeń o stałym kroku `discretization` = 0,01 s (10 ms). Wpis `index[k]` wskazuje dokładny indeks w `dataset_events_t.npy` dla chwili $t_0 + k \cdot 0,01$ s, co umożliwia natychmiastowy dostęp czasowy $O(1)$ bez przeszukiwania binarnego $O(\log N)$ milionów zdarzeń.
- Opóźnienie IMU względem GT: ROZSTRZYGNIĘTE → stałe w obrębie sekwencji ($\le 13$ s, dryf $\le 0,15$ ms); wykazuje liniowy dryf zegara w sesji ok. +11 ppm; estymowane per nagranie.
- Opóźnienie IMU względem zdarzeń z maksymalizacji kontrastu: NIEROZSTRZYGNIĘTE → czułość w oknach 30–100 ms jest rzędu ułamka piksela, co uniemożliwia stabilną estymację online; do ewaluacji algorytmu należy przyjąć skalibrowane opóźnienie z GT.
- Wpływ paralaksy i model translacji 3D: ROZSTRZYGNIĘTE (Oracle) → Dodanie translacji kamery i głębi GT per zdarzenie redukuje FPR tła 2,5–8× (do 0,008–0,026) i podnosi precyzję do 0,74–0,81 (IoU do 0,28–0,40). Rotacja z IMU jest w pełni równoważna rotacji z GT. Cały nieskompensowany błąd modelu rotacyjnego wynikał z paralaksy.

## Przegląd 19 sekwencji (survey_sequences.py)
- Prędkość kątowa $\omega$ kamery z GT: maks. od 0,04 do 2,76 rad/s; czas trwania sekwencji: 1,4–13 s.
- Przepustowość maksymalna w oknie 10 ms wynosi 5,7–6,04 Mev/s we wszystkich 19 sekwencjach (sprzętowy sufit interfejsu kamery). Do wymagań determinizmu czasu rzeczywistego w C++/ARM należy przyjąć 6 Mev/s jako szczytowe obciążenie strumienia.
- Globalna kolumna metadanych „ruchome*” jest bezużyteczna (prawie wszystkie obiekty oznaczono jako ruchome) — nie stosować do doboru sekwencji.
- Wiarygodne dopasowanie żyroskopu do GT (residuum $< 12\%$, bootstrap p95 $< 1,2^\circ$): 10 z 19 sekwencji (`scene13_01`, `02`, `03`, `05`; `scene14_00`, `01`, `03`, `04`, `05`; `scene15_03`). Pozostałe sekwencje posiadają zbyt małą dynamikę ruchu (residuum 20–55%).

## Kalibracja IMU względem kamery i zegary (gyro_vs_gt.py)
- Konwencja poz A ($T_{świat \leftarrow kamera}$): stosunek residuum do konwencji B $\ge 2,4$ we wszystkich sekwencjach; konwencja B daje ujemny wyznacznik i niefizyczne obroty.
- Macierz ekstrynsyk $R_{imu \to kamera} \approx \text{diag}(-1, -1, 1)$ odpowiada odwróconemu montażowi układu IMU na płytce PCB kamery. Odchylenia między sekwencjami wynoszą 0,3–0,8°.
- Pomiary opóźnienia IMU vs GT: `scene13`: −30,9…−21,3 ms; `scene14`: −17,8…−15,6 ms; `scene15`: −5,6…+0,1 ms; `scene6_03`: +14,3 ms (inna sesja, 22 dni wcześniej).
- Liniowy wzrost opóźnienia w czasie sesji wynosi +11,0 ppm (rozrzut 0,44 ms, maks. 1,28 ms), co odpowiada typowemu dryfowi niesynchronizowanych oscylatorów kwarcowych (~40 ms/h).
- W obrębie sekwencji ($\le 13$ s) dryf wynosi $\le 0,15$ ms, dlatego opóźnienie jest stałe podczas jednego nagrania. Ignorowanie opóźnienia generuje do ok. 28 px błędu rotacji przy $\omega = 2$ rad/s.

## Kompensacja rotacji zdarzeń z IMU (rotcomp.py, rotcomp_eval.py)
- Kompensacja z wykorzystaniem $R$, biasu i opóźnienia usuwa 90–97% pozornego przemieszczenia od ego-rotacji w oknach z ruchem $\ge 2$ px (błąd p95 spada z 2,09–4,18 px do 0,11–0,38 px).
- Wariancja obrazu zdarzeń (IWE) po kompensacji rośnie 2,7–3,2× w oknach 30 ms i 4,2× w oknie 100 ms (`scene13_02`), co potwierdza poprawność modelu odwracania dystorsji RadTan i całkowania $SO(3)$.
- Wrażliwość błędu na odchyłkę opóźnienia w oknach 30 ms jest mała: błąd $\pm 25$ ms daje jedynie 0,56–0,58 px błędu p95.
- Zastosowanie stałego opóźnienia dla zbioru (mediana leave-one-out z innych sekwencji: ok. −16 ms) podnosi błąd p95 zaledwie o 0,25–0,34 px względem wyroczni GT i kosztuje $\le 0,6$ pkt proc. FPR detektora.
- Zysk ostrości (wariancja IWE) przy małych przemieszczeniach (30 ms) bywa bliski 1,0 nawet dla orientacji z GT — nie nadaje się jako miara ewaluacyjna; miarą wiodącą jest błąd pikselowy względem GT.

## Opóźnienie IMU względem zdarzeń (event_lag_check.py)
- Maksymalizacja wariancji IWE w oknach 30 i 100 ms daje rozrzucone optima (std 22–32 ms); zaledwie 26% okien trafia w granicach 10 ms od $d_{GT}$ (mediana błędu okna 24,6 ms).
- Analiza czułości ($f \cdot \dot{\omega} \cdot \Delta d \cdot \tau$): dla typowego $\dot{\omega} \approx 3$ rad/s² błąd opóźnienia 5 ms przesuwa zdarzenie na krańcu okna 30 ms zaledwie o ok. 0,12 px (0,46 px w 100 ms), co leży poniżej rozdzielczości metody i ginie w szumie dyskretyzacji oraz paralaksie.
- Decyzja projektowa: do dalszych prac i oceny detekcji należy stosować opóźnienie z GT per sekwencja (lub stałą medianę leave-one-out) jako skalibrowany parametr wejściowy.

## Pozy obiektów, ruch w czasie i etykiety IMO (object_motion.py, imo_labels.py)
- Pozy obiektów w metadanych są rejestrowane w układzie kamery ($T_{kamera \leftarrow obiekt}$, klasa $R*S$). Obiekty statyczne w świecie wykazują przemieszczenie p95 $< 0,5$ cm i szum rotacji do 7,3°.
- Prędkości obiektów w świecie należy wyznaczać różnicą centralną pozycji po złożeniu $T_{w \leftarrow o} = T_{w \leftarrow c} \cdot T_{c \leftarrow o}$ z progiem prędkości liniowej ($v > 0,05$ m/s) lub kątowej ($\omega > 1,0$ rad/s).
- Duża część obiektów porusza się z małymi prędkościami ($v \in [0,05; 0,3]$ m/s, np. epizodyczne przemieszczenia dużych elementów o ID 9 i 22 o 1–6 cm); obiekty te zajmują do 30–45% kadru, lecz generują znikomy strumień zdarzeń. Prawdziwe IMO generujące wyraźny sygnał mają $v > 0,3$ m/s. Wyniki detekcji należy raportować osobno dla szybkich IMO ($v > 0,3$ m/s) oraz wszystkich ruchomych ($v > 0,05$ m/s).
- Sekwencja `scene6_03` stanowi w pełni zwalidowany przypadek referencyjny czystego IMO przy nieruchomej kamerze: gęstość zdarzeń w masce IMO jest 27,7× wyższa niż w tle (87% zdarzeń w 19% powierzchni). Ruch ma charakter impulsowy (początek w 0,65 s).
- W całym zbiorze EVIMO2 dynamika kamery i obecność IMO są silnie ujemnie skorelowane (Spearman −0,78); brak sekwencji łączących jednocześnie gwałtowny obrót kamery i duży udział IMO.

## Wybór sekwencji benchmarkowych do pracy magisterskiej
- `scene6_03`: nieruchoma kamera, impulsowy ruch wielu obiektów — wzorzec czystej segmentacji bez konieczności kompensacji ego-ruchu (walidacja 27,7×).
- `scene14_05` oraz `scene14_04`: dominujący obrót, wiarygodna kalibracja (5,4% i 9,1%) — sekwencje referencyjne do oceny redukcji FPR dzięki kompensacji rotacji (IMO 0,9–1,0%).
- `scene13_01`: umiarkowany obrót i translacja — weryfikacja wpływu paralaksy na fałszywe alarmy (korelacja FPR z prędkością kamery $+0,75$).
- `scene15_03`: duża translacja liniowa — badanie ograniczeń modelu czystej rotacji i konieczności uwzględnienia głębi sceny.
- `scene15_05` oraz `scene15_02`: mały ruch kamery, wysoki udział IMO — sekwencje do oceny czułości wykrywania obiektów.
- `scene13_02` oraz `scene14_00`: długie sekwencje z silnym ruchem kamery i znikomym IMO ($< 0,3\%$) — test odporności na fałszywe alarmy (FPR) oraz test obciążeniowy throughputu (do 6 Mev/s).
- `scene13_05`: odrzucona z benchmarku z powodu zaledwie 13 dostępnych masek GT.

## Detektor IMO v0 i fuzja danych (imo_detect.py)
- Algorytm v0: kompensacja rotacji IMU $\to$ podział okna na połówki wczesną i późną $\to$ wygładzona mapa niespójności $D = |E - L| / (E + L)$ $\to$ progowanie gęstościowe $S \ge \text{support} \cdot \bar{S}$ $\to$ klasteryzacja składowych (CCL na komórkach $4 \times 4$ px).
- Kompensacja rotacji zmniejsza wskaźnik fałszywych alarmów (FPR) 1,8–3,9× (np. w `scene14_05` z 0,277 do 0,075; w `scene13_01` z 0,151 do 0,064).
- Bezwzględne IoU zdarzeń wynosi 0,10–0,25. Niski poziom wynika ze specyfiki kamer neuromorficznych (obiekt generuje zdarzenia wyłącznie na krawędziach, wnętrza pozostają puste — zjawisko apertury) oraz ze znoszenia się powolnych obiektów przy różnicowaniu połówek okna.
- Paralaksa od translacji kamery stanowi główne źródło nieskompensowanych fałszywych alarmów (potwierdzone silną korelacją FPR z prędkością kamery w `scene13_01`: $+0,75$). Wskazuje to na konieczność rozbudowy modelu o fuzję z informacją o głębi/translacji.
- Etykietowanie temporalne (`--gt-mode temporal`): przypisywanie zdarzeniom etykiety z najbliższej w czasie klatki z maską ($\le 12$ ms) jest niezbędne, aby wyeliminować pozorne błędy FPR wynikające z przesunięcia szybkiego obiektu w czasie trwania okna.
- Wypełnianie wnętrz (`--fill`) na siatce komórek wymaga domkniętych obrysów; domknięcie morfologiczne (`--close 1/2`) podnosi czułość kosztem FPR wzdłuż tej samej krzywej ROC.

## Wizualizacja Ground Truth i narzędzia (gt_motion_viz.py) [ZROBIONE]
- Zaimplementowano dedykowany skrypt `gt_motion_viz.py` generujący dwupanelowe zestawienia czystego Ground Truth bez nakładania predykcji detektora:
  * Panel lewy: surowy strumień zdarzeń (polaryzacja ON/OFF) z nałożonymi obrysami obiektów z GT pokolorowanymi według stanu dynamiki (zielony = szybkie IMO $v > 0,3$ m/s, pomarańczowy = obiekty wolne $0,05 < v \le 0,3$ m/s, ciemnoszary = obiekty statyczne) oraz etykietami prędkości liniowej $v$ w m/s.
  * Panel prawy: czyste etykiety zdarzeń z GT (jaskrawy zielony = prawdziwe zdarzenia IMO [czysty TP z GT], pomarańczowy = zdarzenia obiektów wolnych [strefa niejednoznaczna], niebieski = zdarzenia obiektów statycznych, jasnoszary = tło).

## Eksperyment Oracle translacji i głębi (imo_detect.py --translation)
- Model geometryczny kompensacji 3D: $\mathbf{P}_{ref} = M_i \mathbf{P}_i + R_{wc}(t_{ref})^T (\mathbf{p}_c(t_i) - \mathbf{p}_c(t_{ref}))$, gdzie $\mathbf{P}_i = Z_i (u_i, v_i, 1)^T$ ($Z_i$ z mapy głębi GT w mm przeliczana na metry), a $M_i = R_{wc}(t_{ref})^T R_{wc}(t_i)$ pochodzi z całkowania żyroskopu IMU lub z kwaternionów GT. Następnie punkty rzutowane są perspektywicznie na płaszczyznę obrazu.
- Wyniki na 4 sekwencjach referencyjnych (`13_01`, `14_05`, `14_04`, `15_03`): dodanie translacji i głębi GT obniża FPR tła 2,5–8× względem samej rotacji (`13_01`: 0,065 $\to$ 0,008; `14_05`: 0,081 $\to$ 0,026; `14_04`: 0,108 $\to$ 0,020; `15_03`: 0,130 $\to$ 0,052), a łącznie względem braku kompensacji aż 6–20×.
- Precyzja detekcji przy małym udziale IMO ($P/N \approx 0,05$–$0,15$) rośnie z 0,19–0,37 do 0,74–0,81 (IoU zdarzeń rośnie z 0,14–0,25 do 0,28–0,40), co dowodzi, że niska precyzja modelu czystej rotacji wynikała z fałszywych alarmów paralaksy tła.
- Rotacja z IMU jest w pełni równoważna rotacji z GT: wariant "rotacja IMU + translacja" daje identyczny FPR i precyzję jak "rotacja+translacja GT" (np. w `14_05`: FPR 0,026 vs 0,026; w `14_04`: 0,020 vs 0,020; w `13_01`: 0,008 vs 0,008). Kalibracja żyroskopu całkowicie eliminuje błąd składowej rotacyjnej.
- Korelacja Spearmana FPR z prędkością kamery w `scene13_01` spada z +0,75 do +0,18, co ostatecznie potwierdza, że źródłem nieskompensowanego sygnału w tle była paralaksa translacyjna.
- W sekwencjach `14_05`, `14_04`, `15_03` po kompensacji translacji korelacja z prędkością kamery staje się ujemna (−0,44…−0,66): przy szybszym ruchu rośnie liczba zdarzeń generowanych przez tło ($N$), co przy stabilnym poziomie szumu obniża wskaźnik $FP/N$.
- Resztkowy poziom FPR (1–5%) po pełnej kompensacji stanowi dolną granicę szumową (wynika ze zjawiska apertury, szumu sensora DVS, niepełnego pokrycia map głębi oraz tolerancji krawędziowej masek GT).
- Pokrycie map głębi w EVIMO2: głębia GT z LiDAR/RGB-D dostępna jest dla 61–84% zdarzeń w oknie (reszta pikseli przyjmuje medianę okna), co oznacza, że zysk z idealnej głębi ciągłej mógłby być jeszcze wyższy.
- Decyzja projektowa i metodyczna: wariant Oracle (GT w pętli) stanowi w pracy magisterskiej teoretyczną górną granicę możliwości kompensacji ego-ruchu (upper bound) i formalne rozstrzygnięcie wpływu paralaksy. Samodzielny moduł translacji bez GT pozostaje odłożony ze względu na złożoność i dryf całkowania akcelerometru, a biblioteka `evimo2` zachowuje czysty interfejs do opcjonalnej fuzji głębi i pozycji.

## Implementacja C++ i walidacja numeryczna (cpp/) [ZROBIONE]
- Pełny port potoku detekcji runtime z Pythona do C++17 (`camera.hpp`, `calibration.hpp`, `rotcomp.hpp`, `detector.hpp`) bez zewnętrznych zależności (brak Eigen i OpenCV, zero-copy mmap przez `npy.hpp`).
- Weryfikacja numeryczna (`test_numerics.cpp`): 100% zgodności z referencją pythonową na 233 144 zdarzeniach (maks. błąd promieni $1,1 \cdot 10^{-16}$, macierzy $1,0 \cdot 10^{-15}$, współrzędnych $3,0 \cdot 10^{-5}$ px, 0 / 233 144 błędów klasyfikacji IMO).
- Wydajność i throughput (`bin/imo_detect`): przetworzenie 3,08 miliona zdarzeń w 465 ms (średnio 13,7 ms na okno 40 ms; przepustowość 6,61 Mev/s). Przekracza szczytową przepustowość czujnika Prophesee (~6 Mev/s) i osiąga współczynnik czasu rzeczywistego 2,92× na pojedynczym rdzeniu CPU.


