- Średnia przepustowość w tym oknie (4,3 Mev/s) jest ponad dwa razy wyższa od średniej z całej sekwencji (1,93 Mev/s), a szczyt z 10 ms-owego okna dla całej sekwencji to 5,86 Mev/s. Wymagania czasu rzeczywistego trzeba więc liczyć na szczyt, a nie na średnią. Dla pełnego obrazu przydałaby się statystyka z całej sekwencji w oknach o stałej długości.
- Format danych w tej sekwencji to t jako float32 (N,), a nie (N,1) jak w dokumentacji datasetu. Precyzja czasu przy końcu nagrania to ok. 0,5 µs.
- Zdarzenia, IMU i klatki GT są w jednej osi czasu, a ros_time_offset nie jest potrzebny.
- Ruch kamery w sekwencji scene15_dyn_test_00 jest mały (ω maks. ok. 0,17 rad/s), więc sekwencja jest słaba do testu kompensacji rotacji.
- Maski zawierają wszystkie obiekty, w tym statyczne. Etykietę „niezależnie poruszający się” trzeba wyliczyć z poz.
- Otwarte pytania: konwencja poz (T_świat←kamera czy odwrotnie), ekstrynsyka IMU→kamera, prawdziwe opóźnienie IMU (heurystyka dała −15 ms, ale na słabym ruchu), znaczenie tablicy index (654 elementy) i jednostki głębi.
- Sekwencja scene6_dyn_train_03 ma 1,717 s i 2,53 mln zdarzeń. Średnia to 1,47 Mev/s, a maksimum w oknie 10 ms to 5,70 Mev/s.
- Wynik jest bardzo zbliżony do scene15_00 (średnio 1,93, szczyt 5,86 Mev/s), więc na razie na dwóch sekwencjach rząd wielkości szczytu to ok. 5,7–5,9 Mev/s. To dopiero dwa punkty danych - nie traktujmy tego jako ogólnej reguły.
- Mediana (0,80) jest dwa razy niższa od średniej (1,47), a p95 to 4,85. Strumień jest więc bardzo nierównomierny: większość czasu jest spokojna, a krótkie okresy ruchu dają wielokrotnie wyższe przepustowości. To argument za projektowaniem pod szczyt, a nie pod średnią.
- Konwencja poz: A (T_świat←kamera), potwierdzona na 19/19 sekwencji, stosunek residuum co najmniej 2,4.
- Rotacja IMU→kamera ≈ obrót o 180° wokół osi z kamery, odchylenia poniżej 1°. Bias żyroskopu ≈ (−0,0295; 0,0154; 0,0180) rad/s.
- Opóźnienie IMU względem GT zmienia się liniowo z czasem nagrania o ok. +11 ppm w obrębie sesji, a w obrębie jednej sekwencji jest stałe. Trzeba je estymować osobno dla każdej sekwencji (wartość początkowa od −31 do 0 ms). Zmierzone względem GT, a nie względem zdarzeń --> Do sprawdzenia.

# Shortlista z różnymi reżimami ruchu i różnymi scenami
- Bardzo dynamiczna sekwencja to scene13_dyn_test_02, która ma 13 s nagrania, predkosc kątowa do 2,76 rad/s (ok. 15 px na 10 ms od samej rotacji), obrót kamery o 112° i v do 1,2 m/s.
- Duża translacja: scene15_03
- Obrót dominujący: scene14_05 [jest to z folderu eval, narazie nie tranformowałem tego]
- Umiarkowany ruch: scene14_00 (w 1,4, 8 s)
- Baseline (prawie statyczna kamera, porusza się tylko scena): scene15_00
- Test przepustowości: scene13_05, która ma (4,89 Mev/s, 2,4 s) [jest to z folderu eval, narazie nie tranformowałem tego]