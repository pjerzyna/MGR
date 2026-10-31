## Korekty istniejących wpisów
- Wpis „Zdarzenia, IMU i klatki GT są w jednej osi czasu”: zakresy czasu zdarzeń, IMU i klatek GT pokrywają się (start ≈ 0,017 s) i ros_time_offset nie jest potrzebny, ale znaczniki IMU są przesunięte względem GT o ok. −31…+14 ms (patrz niżej).
- Z listy „Otwarte pytania”: konwencja poz → ROZSTRZYGNIĘTE (A); ekstrynsyka IMU→kamera → ROZSTRZYGNIĘTE (niżej); opóźnienie IMU → zmierzone względem GT, względem zdarzeń NIEROZSTRZYGNIĘTE. Nadal otwarte: tablica `index` (654 elementy), jednostki głębi (wartości ok. 550 – mm?).

## Przegląd 19 sekwencji (survey_sequences.py)
- ω z GT maks. od 0,04 do 2,76 rad/s; czas 1,4–13 s. Przepustowość maks. w oknie 10 ms: 5,7–6,04 Mev/s, niemal identyczna we wszystkich 19. Hipoteza (niepotwierdzona): sufit sprzętowy. Do wymagań czasu rzeczywistego przyjąć ok. 6 Mev/s jako szczyt.
- Kolumna „ruchome*” z przeglądu jest bezwartościowa (prawie wszystkie obiekty oznaczone jako ruchome) – nie używać do wyboru sekwencji.
- Wiarygodne dopasowanie żyroskopu do GT (residuum < 12%, bootstrap p95 < 1,2°): scene13_01, 02, 03, 05; scene14_00, 01, 03, 04, 05; scene15_03 (10 z 19). Pozostałe mają za mało ruchu (residuum 20–55%).

## Kalibracja IMU względem kamery (gyro_vs_gt.py) [USTALONE]
- Konwencja poz A (T_świat←kamera, t = pozycja kamery): 19/19 sekwencji, stosunek residuum ≥ 2,4. Konwencja B daje za każdym razem absurdalny obrót i ujemny wyznacznik.
- R_imu→kamera ≈ obrót o 180° wokół osi z (x→−x, y→−y, z→z); odchylenia 0,3–0,8° między sekwencjami (bootstrap 0,05–0,2° to niedoszacowanie, realnie ok. 0,5°).
- Bias żyroskopu ≈ (−0,0295; 0,0154; 0,0180) rad/s (|b| ≈ 0,038), stały na 4 sekwencjach (±0,001).
- Pozy obiektów są prawdopodobnie w układzie kamery (ω p95 pięciu obiektów w scene14_05 równe ω kamery) [HIPOTEZA]. Tabela ruchu obiektów w gyro_vs_gt.py zakłada układ świata, więc jest niewiarygodna; wymaga złożenia z pozą kamery.

## Opóźnienie IMU względem GT [USTALONE dla GT]
- scene13: −30,9…−21,3 ms; scene14: −17,8…−15,6; scene15: −5,6…+0,1 (rośnie z indeksem = czasem nagrania). scene6_03: +14,3 (inna sesja, słabe dopasowanie).
- Liniowy wzrost w sesji ok. +11 ppm (scena 13: +10,0; scena 14: +11,8; 18 sekwencji scen 13–15 razem: +11,0 ppm, rozrzut 0,44 ms, max 1,28 ms). Sesja ≈ 45 min; między sesjami reset (scene6 jest 22 dni wcześniej). Który zegar dryfuje – nieznane [HIPOTEZA oparta na 3 znanych ros_time_offset i monotoniczności].
- W obrębie sekwencji (≤13 s) dryf ≤ 0,15 ms, więc opóźnienie jest stałe w sekwencji i trzeba je estymować osobno dla każdej (wartość początkowa −31…0 ms). Szacunek: ignorowanie opóźnienia = do ok. 28 px błędu rotacji przy ω = 2 rad/s (nie zmierzone).

## Kompensacja rotacji zdarzeń z IMU [USTALONE]
- Z R, biasem i opóźnieniem z GT wariancja obrazu zdarzeń rośnie 2,7–3,2× (okno 30 ms; scene13_02, 14_05, 15_03) i 4,2× (okno 100 ms, scene13_02); krawędzie wyraźnie ostre. To pierwsza weryfikacja R, biasu i modelu kamery (w tym radtan) na prawdziwych zdarzeniach.

## Opóźnienie względem zdarzeń [NIEROZSTRZYGNIĘTE]
- Maksymalizacja kontrastu (model czystego obrotu), okna 30 i 100 ms: optima per okno rozrzucone po całym zakresie (std 22–32 ms), 1–4 okna na sekwencję na brzegu siatki. Tylko 8 z 31 okien (26%) ma optimum w granicach 10 ms od d_GT; mediana błędu okna 24,6 ms.
- Średnia krzywa: d_ev − d_GT = −0,5 / +12,5 / −1,8 ms przy błędach ok. ±6 / ±12 / ±16 ms; chi² 1,1 (d_ev = d_GT, 3 st. sw.) vs 5,0 (d_ev = const, 2 st. sw.); --summary: NIEJEDNOZNACZNY. Wniosek „nachylenie bliskie 1” ze starej wersji skryptu był fałszywie pewny – NIE cytować.
- Diagnostyka (31 okien, Spearman, próg ok. ±0,36): błąd okna nie koreluje z |ω| (−0,17), ω̇ (−0,10) ani z |v| kamery (+0,05). Hipoteza translacji jako głównej przyczyny NIE znajduje poparcia; hipoteza „opóźnienie widać tylko przy dużym ω̇” również nie (okna z ω̇ ≥ 5 rad/s²: 4/10 dobrych, pozostałe 4/21; Fisher p = 0,21). Ograniczenia: n = 31; |v| to średnia prędkość kamery, a paralaksa zależy od głębi poszczególnych pikseli; efekty umiarkowane nie są wykluczone.
- Czułość (szacunek liniowy f·ω̇·Δd·τ, ω̇ w oknach 0,5–10,5 rad/s², mediana ok. 3): błąd opóźnienia 5 ms zmienia kompensację rotacji na krańcach okna 100 ms o 0,46 px (mediana; maks. 1,5 px); błąd 25 ms o ok. 2,3 px (mediana). Rozdzielczość metody rzędu 20–30 ms jest zgodna z małą czułością pojedynczych okien.
- Nieprzetestowane źródła rozrzutu: ruchome obiekty (test: tylko zdarzenia z pikseli tła według masek GT), błąd skali żyroskopu, zależność od głębi pikseli.
- Uwaga: okna 30 ms leżą wewnątrz okien 100 ms o tych samych środkach (wspólne zdarzenia), więc ich zgodność (r = 0,96) NIE jest niezależnym dowodem błędu systematycznego.

## Decyzje robocze (do potwierdzenia)
- Do rozwoju i oceny używać opóźnienia z GT per sekwencja jako skalibrowanego parametru; estymacja ze zdarzeń jako osobne zagadnienie.
- Konfiguracja: wyniki w $EVIMO2_OUT/sceneXX_NN/.