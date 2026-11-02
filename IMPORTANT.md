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

## Pozy obiektów
- Konwencja poz obiektów (T_kamera←obiekt, klasa R*S): w 4 sekwencjach z JSON i obrotem kamery >= 10° stałych jest odpowiednio 4/8, 7/9, 4/8, 3/8 obiektów w tej klasie i 0 w każdej innej [USTALONE]. Do potwierdzenia na pełnych 19 sekwencjach.

- Konwencja poz obiektów T_kamera←obiekt (klasa R*S): 11/11 sekwencji z obrotem kamery >= 10°, w każdej 3–7 obiektów stałych w tej klasie i 0 w każdej innej [USTALONE].
- Statyczne obiekty mają przemieszczenie p95 < 0,5 cm (100 z 163 obiektów; luka 0,5–1 cm), szum rotacji statycznych sięga 7,3°. Ruchome nie mają wyraźnej luki od statycznych poniżej 10 cm (16 obiektów 1–10 cm). Globalna flaga „ruchomy” nie nadaje się jako ground truth: id 9 (ok. 30% obrazu) jest oznaczony w 9 sekwencjach przy przemieszczeniu 1,2–6,3 cm. Potrzebne etykiety w czasie z progiem prędkości (statyczne v p95 <= 1,7 cm/s, ruchome ok. 1,1–1,3 m/s).
- scene6_03 (klatki zdarzeń co 100 ms): zdarzenia skupione na kilku strukturach przesuwających się po obrazie (duży prostokątny element z uchwytem, okrągły obiekt, panel przy krawędzi), tło w większości bez zdarzeń, klatki 0000 i 0050 prawie puste. Zgodne z ruchomą sceną i nieruchomą kamerą; ruch jest impulsowy, więc potrzebne etykiety zależne od czasu. Dowód jakościowy (nie ilościowy) - do potwierdzenia walidacją gęstości zdarzeń w masce IMO.
- Etykiety ruchu w czasie: imo_labels.py (v > 5 cm/s lub ω > 1 rad/s, progi wstępne); żadne etykiety ani sekwencje nie są wykluczane, jakość oznaczana walidacją.
- scene6_03 ZWALIDOWANA: gęstość zdarzeń w masce IMO / poza = 27,7× (87% zdarzeń w 19% powierzchni), nieruchoma kamera. Ruch obiektów impulsowy: brak ruchu do 0,65 s; klatki 39–54 12–13 z 13 obiektów naraz; luka masek dla klatek 61–89 (1,02–1,5 s). Czysty przypadek IMO (bez ruchu własnego).
- W zbiorze ruch kamery i obecność IMO są odwrotnie skorelowane (Spearman −0,78, 19 sekwencji): ω_gt >= 0,8 rad/s: IMO średnio 1,3%; ω_gt < 0,3 rad/s: 7,7%. Nie ma sekwencji z silnym ruchem kamery i dużym IMO.
- Walidacja zdarzeniowa etykiet IMO: wiarygodna tylko przy słabym ruchu kamery (stosunek gęstości 2,7–13,8×, scene6_03 27,7×); przy silnym ruchu ≈ 1 (tło generuje zdarzenia), nie rozstrzyga.
- Podejrzenie: p95 udziału IMO 33–46% (14_02, 14_03, 15_02, 15_04, 15_05, 15_06) - flagowanie dużych obiektów (9, 22) w krótkich epizodach; do rozstrzygnięcia rozbiciem „tylko ω” vs „v” w imo_labels.py.
- Próg prędkości: płytkie minimum w histogramie przy 0,05–0,1 m/s (449 par obiekt–klatka); v >= 0,05: 8,6% par, v >= 0,1: 7,3%. Wybór progu do potwierdzenia sweepem (0,03/0,05/0,1/0,2) i walidacją.
- scene13_05 ma tylko 13 masek - nie nadaje się do oceny.
- Rozbicie flag ruchu: oznaczenia wyłącznie przez ω (v <= 5 cm/s) to margines (<= 0,5 pkt proc. IMO w najgorszych sekwencjach: 15_05 5,79→5,30%, 15_06 7,92→7,50%; w pozostałych <= 0,2 pkt). Wysokie p95 IMO (33–46%) w 14_02, 14_03, 15_02, 15_04, 15_05, 15_06 pochodzi z prędkości liniowej dużych obiektów (9, 22) > 5 cm/s, nie z szumu rotacji. Hipoteza: powolne ruchy (5–10 cm/s, ok. 0,3 s) dużych obiektów, zgodne z ich całkowitym przemieszczeniem 1–6 cm; wspierane walidacją zdarzeniową (7,3×, 9,0×, 13,8× przy słabym ruchu kamery).
- Konsekwencja dla oceny: obiekty powolne (v rzędu 5–30 cm/s) generują bardzo mało zdarzeń, a mogą zajmować 30–45% obrazu, więc dominowałyby w IoU. Raportować wyniki osobno dla IMO szybkich (v > 0,3 m/s) i wszystkich (v > 0,05 m/s); żadne etykiety nie są usuwane.
- scene14_00: praktycznie brak IMO (średnio 0,08%, walidacja 0,9×) - test fałszywych alarmów, nie ocena detekcji.

## Przyszły wybór sewkencji???
- scene6_03	nieruchoma kamera, dużo IMO	zwalidowana (27,7×), czysta segmentacja bez kompensacji
- scene15_05	słaby ruch, dużo IMO	295 masek, IMO 5,8%, walidacja 9,0×
- scene15_00	słaby ruch, dużo IMO	393 maski, IMO 8,7%, baseline z wcześniejszej listy
- scene14_05	obrót dominujący	kalibracja wiarygodna (5,4%), IMO 0,9%, 252 maski
- scene14_04	umiarkowany ruch	kalibracja 9,1%, IMO 1,0%, 285 masek
- scene14_03	umiarkowany ruch, większe IMO	kalibracja 10,9%, IMO 2,7% (do sprawdzenia flagowanie dużych obiektów)
- scene15_03	duża translacja	kalibracja 4,0%, IMO 2,0%, walidacja 9,4×
- scene13_01	umiarkowany ruch	kalibracja 7,5%, IMO 1,2%, 392 maski
- scene13_02	skrajny ruch, test fałszywych alarmów	13 s, 112°, 750 masek, IMO 0,2%
- scene13_03 lub 14_00	długi, silny ruch, test przepustowości	12,8 s lub 8 s; IMO < 0,3%

## Kompensacja rotacji pierwsze wnioski
- Udział IMO w funkcji progu prędkości: serie 15, 14_02 i 14_03 mają głównie obiekty wolne (< 0,3 m/s): scene15_00 8,71% (v>0,05) → 0,08% (v>0,3). Obiekty szybkie (v>0,3): scene6_03 3,57%, 13_00 3,66%, 13_01 0,83%, 14_04 0,85%, 14_05 0,79%, 15_01 1,41%, 15_02 1,36%, 15_05 1,30%, 15_06 1,02%. Propozycja sekwencji do oceny detekcji: 6_03, 13_00, 13_01, 14_03, 14_04, 14_05, 15_05, 15_02; scene15_00 tylko jako test fałszywych alarmów.
- Kompensacja rotacji v0 (rotcomp.py, stałe kalibracyjne R = diag(-1,-1,1), bias; bez GT w przetwarzaniu): błąd względem orientacji z GT, mediana po 36–40 oknach 30 ms (p50/p95 [px]): scene13_02: bez kompensacji 1,00/2,09, oracle 0,08/0,18, stałe −15,7 ms 0,15/0,38, zero 0,24/0,62; scene14_05: 0,31/0,73 → 0,08/0,19 (stałe = oracle, bo d_GT = −15,6); scene15_03: 0,16/0,35 → 0,05/0,11 (wszystkie warianty ≈ 0,05 = dolna granica pomiaru); scene14_04: 0,25/0,51 → stałe 0,06/0,20, zero 0,09/0,29.
- Wrażliwość na opóźnienie przy oknach 30 ms jest mała: ±25 ms daje 0,56–0,58 px p95 w scene13_02 (zgodne z f·ω̇·Δd·τ ~ 0,6–0,8 px); błąd rośnie ∝ długości okna (prognoza: ~3× przy 100 ms, niezmierzona).
- Mediana po oknach jest zdominowana przez spokojne okna (typowa ego-rotacja w 30 ms to 0,16–0,31 px w 3 z 4 sekwencji); rozkład i tabela dla okien z p95 ≥ 2 px w rotcomp_eval.py (do uzupełnienia po ponownym uruchomieniu).
- Zysk ostrości obrazu zdarzeń jest ≈1,0 nawet dla orientacji z GT (14_05: 1,04; 15_03: 0,99; 14_04: 1,02; 13_02: 1,35) - nie nadaje się jako miara przy oknach 30 ms; podstawowa miara to błąd vs GT. Zysk IMO ≈1,0, ale w 13_02 tylko 1 okno z szybkimi IMO.
- Okna z ego-rotacją >= 2 px (p95 bez kompensacji): scene14_05 39% okien, scene15_03 22%, scene13_02 52%, scene14_04 26%; >= 5 px: 22%, 18%, 22%, 10%. W tych oknach kompensacja z IMU usuwa 90–97% przemieszczenia (błąd p95 w pikselach: oracle 0,19–0,39; stałe 0,21–0,58; zero 0,25–0,86; ±25 ms 0,66–1,79). Tabele z okien o słabej rotacji zbiegają do poziomu szumu (ok. 0,05 px). Wpływ translacji nie był mierzony.
- Wariant „stałe d” w rotcomp_eval dla 14_04 i 14_05 nie jest niezależny od oracle (mediana 19 sekwencji = d_GT scene14_04). Od teraz stałe d = mediana d_GT z INNYCH sekwencji (leave-one-out, >= 3); uczciwy koszt braku synchronizacji: scene13_02 +0,25 px, scene15_03 +0,34 px (p95).
- Detektor IMO v0 (imo_detect.py: kompensacja rotacji -> niespójność dwóch połówek okna -> klasteryzacja składowych): na syntetyku brak kompensacji daje FPR 27,6%, z kompensacją 0,3%; czułość niska (0,25–0,53) z powodu znoszenia się wnętrza obiektu przy różnicowaniu; rozdzielczość detektora (rozmycie) jest sprzężona z dokładnością kompensacji (blur 1: FPR do 30% przy błędzie opóźnienia 28 ms, blur 3: 0,5%). Wyniki na danych rzeczywistych: do uzupełnienia.