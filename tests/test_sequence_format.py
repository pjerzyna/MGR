"""Testy formatu, spójności i ustalonej kalibracji prawdziwej sekwencji EVIMO2v2 (kamera lewa, Prophesee).

Uruchomienie:
    EVIMO2_SEQ_DIR=/ścieżka/do/scene15_dyn_test_00_000000 pytest
    EVIMO2_SEQ_DIR=/home/pawsooon/Documents/Vsc/npz_odczyty/data/npz_left_camera_imo/left_camera/imo/train/scene15_dyn_test_00_000000 pytest
    
Bez zmiennej EVIMO2_SEQ_DIR testy są pomijane (np. po sklonowaniu repo bez danych).
Testy niczego nie modyfikują - tylko czytają wskazaną sekwencję.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

import evimo2_common as ec
import gyro_vs_gt

LEFT_IMU = "/prophesee/left/imu"


@pytest.fixture(scope="module")
def seq_dir():
    path = os.environ.get("EVIMO2_SEQ_DIR")
    if not path:
        pytest.skip("Ustaw EVIMO2_SEQ_DIR na katalog sekwencji EVIMO2, aby uruchomić te testy")
    p = Path(path)
    if not p.is_dir():
        pytest.fail(f"EVIMO2_SEQ_DIR wskazuje na nieistniejący katalog: {p}")
    return p


@pytest.fixture(scope="module")
def info(seq_dir):
    data = ec.load_info(seq_dir)
    assert data is not None, "brak pliku dataset_info.npz"
    return data


def test_events_format(seq_dir, info):
    """Pliki zdarzeń mają kształty, typy i zakresy, na których opiera się reszta potoku."""
    t = np.load(seq_dir / "dataset_events_t.npy", mmap_mode="r")
    xy = np.load(seq_dir / "dataset_events_xy.npy", mmap_mode="r")
    p = np.load(seq_dir / "dataset_events_p.npy", mmap_mode="r")
    n = len(t)

    assert n > 0
    assert t.shape == (n,) and xy.shape == (n, 2) and p.shape == (n,)
    assert t.dtype.kind == "f", f"czas ma typ {t.dtype}"
    assert xy.dtype.kind == "u", f"xy ma typ {xy.dtype}"
    assert p.dtype.kind in "ui", f"polaryzacja ma typ {p.dtype}"
    assert np.isin(np.unique(p), [0, 1]).all(), "polaryzacja inna niż 0/1"

    intr = ec.find_dict_with(info, ["res_x", "res_y"])
    assert intr is not None, "brak res_x/res_y w meta"
    assert xy[:, 0].max() < intr["res_x"] and xy[:, 1].max() < intr["res_y"]
    assert np.all(np.diff(t) >= 0), "czas zdarzeń nie jest niemalejący"


def test_info_structure(info):
    """dataset_info.npz zawiera klucze i pola, z których korzystamy dalej."""
    assert {"index", "discretization", "K", "D", "meta"} <= set(info)
    assert info["K"].shape == (3, 3)
    assert info["D"].shape == (4,)

    meta = info["meta"]
    assert {"frames", "full_trajectory", "imu", "meta"} <= set(meta)
    assert {"fx", "fy", "cx", "cy", "res_x", "res_y"} <= set(meta["meta"])

    assert LEFT_IMU in meta["imu"], f"brak IMU lewej kamery, są: {list(meta['imu'])}"
    sample = meta["imu"][LEFT_IMU][0]
    assert {"angular_velocity", "linear_acceleration", "ts"} <= set(sample)
    assert set("xyz") <= set(sample["angular_velocity"])
    assert set("xyz") <= set(sample["linear_acceleration"])

    assert {"id", "ts", "gt_frame", "cam"} <= set(meta["frames"][0])


def test_imu_and_ground_truth_consistency(seq_dir, info):
    """IMU pokrywa zdarzenia bez luk, a maski/głębia wskazują istniejące klatki z meta."""
    t = np.load(seq_dir / "dataset_events_t.npy", mmap_mode="r")
    e0, e1 = float(t[0]), float(t[-1])

    imu_ts = np.array([s["ts"] for s in info["meta"]["imu"][LEFT_IMU]])
    dt = np.diff(imu_ts)
    overlap = (min(e1, imu_ts[-1]) - max(e0, imu_ts[0])) / (e1 - e0)
    assert overlap >= 0.95, f"IMU pokrywa tylko {overlap:.0%} czasu zdarzeń"
    assert dt.min() > 0, "ts IMU nie jest ściśle rosnący"
    assert dt.max() <= 3 * np.median(dt), f"luka w IMU: {dt.max() * 1e3:.2f} ms"

    frame_ids = {f["id"] for f in info["meta"]["frames"]}
    with np.load(seq_dir / "dataset_mask.npz") as z:
        mask_ids = {int(k.split("_")[-1]) for k in z.files}
        for k in z.files[:10]:                       # próbka, żeby test był szybki
            m = z[k]
            assert (m[m > 0] % 1000 == 0).all(), f"{k}: wartości maski nie są wielokrotnością 1000"
    with np.load(seq_dir / "dataset_depth.npz") as z:
        depth_ids = {int(k.split("_")[-1]) for k in z.files}
    # Maski/głębia istnieją tylko dla części klatek (np. 75 z 104), więc wymagamy
    # zgodności masek z głębią i tego, że każda maska wskazuje istniejącą klatkę meta.
    assert mask_ids == depth_ids, "numery masek i głębi się różnią"
    assert mask_ids <= frame_ids, "maski dla klatek, których nie ma w meta"


def test_gyro_vs_gt_calibration(seq_dir, tmp_path, monkeypatch):
    """Ustalona kalibracja zestawu: konwencja poz A, R_imu->cam ~ obrót o 180 deg wokół osi z, bias.

    Pilnuje wyników, które potok C++ ma przyjąć na stałe. Wymaga sekwencji z wyraźnym ruchem
    (residuum dopasowania < 15%); dla sekwencji z małym ruchem test jest pomijany, bo kalibracja
    z takich danych jest niewiarygodna (np. scene15_00, scene6_03).
    """
    monkeypatch.setattr(sys, "argv", ["gyro_vs_gt.py", str(seq_dir), "--out", str(tmp_path), "--no-plot"])
    gyro_vs_gt.main()
    res = json.loads((tmp_path / "gyro_vs_gt.json").read_text())
    if res["residual_rel"] > 0.15:
        pytest.skip(f"za mało ruchu w tej sekwencji (residuum {res['residual_rel']:.0%}); wybierz np. scene14_05 lub scene13_02")

    assert res["convention"] == "A"
    R = np.array(res["R_imu_to_cam"])
    cos = (np.trace(R @ np.diag([-1.0, -1.0, 1.0])) - 1) / 2          # kąt między R a diag(-1, -1, 1)
    assert np.degrees(np.arccos(np.clip(cos, -1, 1))) < 2.0
    assert np.allclose(res["gyro_bias_imu_rad_s"], [-0.0295, 0.0154, 0.0180], atol=0.01)
    assert -0.05 < res["lag_s"] < 0.02, "opóźnienie IMU poza zmierzonym zakresem -31..+0 ms"