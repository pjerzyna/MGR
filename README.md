# evimo2_tools

Narzędzia w Pythonie do inspekcji i weryfikacji spójności sekwencji EVIMO2v2
(Prophesee, kamera lewa). Wymagania: Python >= 3.9, numpy.

| Skrypt | Cel |
|---|---|
| `inspect_evimo2.py <seq_dir>` | zwięzła inspekcja plików, meta, statystyk zdarzeń, masek i głębi |
| `check_sync.py <seq_dir>` | spójność czasowa zdarzeń / IMU / klatek GT, jakość IMU, mapowanie GT |
| `evimo2_common.py` | wspólne funkcje (wczytywanie, przeszukiwanie struktury meta) |


Bez EVIMO2_SEQ_DIR testy są pomijane. Plik tests/test_sequence_format.py zawiera trzy testy: format zdarzeń (kształty, typy, zakresy, monotoniczność czasu), struktura dataset_info.npz oraz spójność IMU i ground truth (pokrycie czasu, brak luk, numery masek/głębi vs klatki meta).