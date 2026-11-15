"""Testy jednostkowe modułu geometrii epipolarnej i FOE (evimo2.epipolar)."""
import numpy as np
import pytest

from evimo2.epipolar import epipolar_direction, estimate_cluster_shift, filter_clusters_epipolar


def test_epipolar_direction_pure_y_translation():
    """Weryfikacja kierunku wektora epipolarnego przy czystej translacji wzdłuż osi Y."""
    intr = {"fx": 500.0, "fy": 500.0, "cx": 320.0, "cy": 240.0}
    # Kamera porusza się w dół (-Y w układzie kamery: vy = -0.2 m/s)
    v_cam = np.array([0.0, -0.2, 0.0])
    
    # Dla dowolnego punktu, pozorny ruch tła u_dot = 0, v_dot > 0 (kierunek [0, 1])
    u_epi, norm = epipolar_direction(320.0, 240.0, v_cam, intr)
    assert norm > 0
    assert np.allclose(u_epi, [0.0, 1.0], atol=1e-5)


def test_epipolar_direction_pure_z_expansion():
    """Weryfikacja Focus of Expansion (FOE) przy ruchu kamery do przodu (+Z)."""
    intr = {"fx": 500.0, "fy": 500.0, "cx": 320.0, "cy": 240.0}
    # Kamera porusza się do przodu (vz = +1.0 m/s)
    v_cam = np.array([0.0, 0.0, 1.0])
    
    # W centrum obrazu (320, 240), pozorna prędkość wynosi 0 (FOE w centrum)
    u_epi_center, norm_center = epipolar_direction(320.0, 240.0, v_cam, intr)
    assert norm_center < 1e-5
    
    # W prawym górnym rogu (x > cx, y < cy), wektory rozchodzą się promieniście na zewnątrz
    u_epi_tr, norm_tr = epipolar_direction(420.0, 140.0, v_cam, intr)
    assert u_epi_tr[0] > 0  # Ruch w prawo
    assert u_epi_tr[1] < 0  # Ruch w górę


def test_estimate_cluster_shift_synthetic():
    """Weryfikacja estymacji przesunięcia wycinka przez korelację wzajemną."""
    H, W = 100, 100
    Es = np.zeros((H, W), dtype=np.float32)
    Ls = np.zeros((H, W), dtype=np.float32)
    
    # Obiekt kwadratowy w Es w pozycji (40, 40)
    Es[40:50, 40:50] = 1.0
    # W Ls przesunięty o dx = +2, dy = -3 -> pozycja (42, 37)
    Ls[37:47, 42:52] = 1.0
    
    shift, corr = estimate_cluster_shift(Es, Ls, 35, 55, 35, 55, max_shift=5)
    assert shift == (2, -3)
    assert corr > 0.99


def test_filter_clusters_epipolar_separation():
    """Weryfikacja odrzucenia klastra tła (paralaksy) i zachowania obiektu ruchomego."""
    intr = {"fx": 500.0, "fy": 500.0, "cx": 320.0, "cy": 240.0}
    v_cam = np.array([0.0, -0.2, 0.0])  # Wektor epipolarny w [0, 1]
    
    H, W = 480, 640
    Es = np.zeros((H, W), dtype=np.float32)
    Ls = np.zeros((H, W), dtype=np.float32)
    
    # Klaster 1: Tło / paralaksa (przesunięcie w dół: dx = 0, dy = +2)
    Es[100:120, 200:220] = 1.0
    Ls[102:122, 200:220] = 1.0
    
    # Klaster 2: IMO (przesunięcie w lewo: dx = -4, dy = 0)
    Es[300:320, 400:420] = 1.0
    Ls[300:320, 396:416] = 1.0
    
    clusters = [
        {"cid": 1, "bbox": (195, 225, 95, 125)},
        {"cid": 2, "bbox": (390, 425, 295, 325)},
    ]
    
    kept = filter_clusters_epipolar(clusters, Es, Ls, v_cam, intr, min_speed=0.05, angle_thresh_deg=30.0)
    
    # Klaster 1 powinien zostać odrzucony jako paralaksa (angle ~ 0 deg)
    # Klaster 2 powinien zostać zachowany jako IMO (angle ~ 90 deg)
    assert len(kept) == 1
    assert kept[0]["cid"] == 2
    assert kept[0]["angle_deg"] > 80.0
