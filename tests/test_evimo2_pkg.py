"""Testy jednostkowe modułów biblioteki evimo2 (bez konieczności posiadania datasetu)."""
import numpy as np
import pytest

import evimo2
from evimo2 import calibration, common, detector, geometry, kinematics, rotcomp


def test_package_exports():
    """Weryfikacja dostępności kluczowych funkcji w głównym namespace evimo2."""
    assert hasattr(evimo2, "load_info")
    assert hasattr(evimo2, "load_events")
    assert hasattr(evimo2, "warp")
    assert hasattr(evimo2, "iwe")
    assert hasattr(evimo2, "calibrated_gyro")
    assert hasattr(evimo2, "CLASSES")
    assert hasattr(evimo2, "detect_imo")


def test_geometry_quaternions():
    """Poprawność iloczynu Hamiltona, koniugacji i obrotu wektora."""
    # q dla obrotu o 90 deg wokół osi Z: cos(45) + k*sin(45)
    s = np.sin(np.pi / 4.0)
    c = np.cos(np.pi / 4.0)
    qz = np.array([[c, 0.0, 0.0, s]])
    v = np.array([[1.0, 0.0, 0.0]])

    # Obrót wektora [1, 0, 0] o 90 deg wokół Z powinien dać [0, 1, 0]
    v_rot = geometry.qrot(qz, v)
    assert np.allclose(v_rot, [[0.0, 1.0, 0.0]], atol=1e-6)

    # Iloczyn kwaternionu z jego sprzężeniem powinien dać kwaternion jednostkowy [1, 0, 0, 0]
    q_inv = geometry.qconj(qz)
    q_id = geometry.qmul(qz, q_inv)
    assert np.allclose(q_id, [[1.0, 0.0, 0.0, 0.0]], atol=1e-6)


def test_rodrigues_batch():
    """Weryfikacja wzoru Rodrigues dla obrotu zerowego i 90 stopni."""
    # Obrót zerowy -> macierz jednostkowa I3
    rv_zero = np.array([[0.0, 0.0, 0.0]])
    R_zero = geometry.rodrigues_batch(rv_zero)
    assert np.allclose(R_zero[0], np.eye(3))

    # Obrót o 90 deg wokół osi Z
    rv_z90 = np.array([[0.0, 0.0, np.pi / 2.0]])
    R_z90 = geometry.rodrigues_batch(rv_z90)
    expected = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    assert np.allclose(R_z90[0], expected, atol=1e-6)


def test_undistort_pinhole():
    """Gdy współczynniki dystorsji są 0, undistort zwraca znormalizowane współrzędne pinhole."""
    intr = {"fx": 500.0, "fy": 500.0, "cx": 320.0, "cy": 240.0, "k1": 0.0, "k2": 0.0, "p1": 0.0, "p2": 0.0}
    px = np.array([320.0, 570.0])
    py = np.array([240.0, 490.0])
    x, y = geometry.undistort(px, py, intr)
    assert np.allclose(x, [0.0, 0.5])
    assert np.allclose(y, [0.0, 0.5])


def test_calibrated_gyro():
    """Kompensacja biasu i rotacji IMU->kamera: R (omega - b)."""
    R = np.diag([-1.0, -1.0, 1.0])
    bias = np.array([0.1, -0.2, 0.3])
    raw_gyro = np.array([[0.1, -0.2, 0.3], [1.1, -0.2, 0.3]])

    cal = calibration.calibrated_gyro(raw_gyro, R=R, bias=bias)
    # Pierwsza próbka to dokładnie bias -> po odjęciu 0
    assert np.allclose(cal[0], [0.0, 0.0, 0.0])
    # Druga próbka: delta x = 1.0 -> po R diag(-1, -1, 1) powinno dać [-1.0, 0.0, 0.0]
    assert np.allclose(cal[1], [-1.0, 0.0, 0.0])


def test_iwe_accumulation():
    """Weryfikacja tworzenia obrazu zdarzeń IWE (Image of Warped Events)."""
    xs = np.array([100.0, 100.5])
    ys = np.array([50.0, 50.0])
    img = rotcomp.iwe(xs, ys, W=640, H=480)
    assert img.shape == (480, 640)
    # Suma wag w obrazie powinna równać się liczbie zdarzeń
    assert np.isclose(img.sum(), 2.0)
    # Punkt [100.0, 50.0] powinien mieć wagę >= 1.0
    assert img[50, 100] > 0.9


def test_detector_box_blur():
    """Weryfikacja szybkiego filtru pudełkowego."""
    img = np.zeros((20, 20), dtype=float)
    img[10, 10] = 1.0
    blurred = detector.box_blur(img, r=1)
    assert blurred.shape == (20, 20)
    # Promień 1 to okno 3x3 (9 komórek)
    assert np.isclose(blurred.sum(), 1.0)
    assert np.isclose(blurred[10, 10], 1.0 / 9.0)
