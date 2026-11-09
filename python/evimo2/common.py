"""Moduł operacji wejścia/wyjścia, metadanych i ścieżek dla sekwencji EVIMO2."""
import glob
import json
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np


def unwrap(a: Any) -> Any:
    """0-wymiarowe tablice obiektowe (np. słownik meta) -> zwykły obiekt Pythona."""
    if isinstance(a, np.ndarray) and a.shape == () and a.dtype == object:
        return a.item()
    return a


def load_info(seq_dir: Union[str, Path]) -> Optional[Dict[str, Any]]:
    """Wczytuje dataset_info.npz do słownika {klucz: obiekt}."""
    path = os.path.join(seq_dir, "dataset_info.npz")
    if not os.path.isfile(path):
        return None
    with np.load(path, allow_pickle=True) as z:
        return {k: unwrap(z[k]) for k in z.files}


def load_events(seq_dir: Union[str, Path]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Zwraca (t, xy, p) jako widoki mmap o kształtach (N,), (N,2), (N,)."""
    t = np.load(os.path.join(seq_dir, "dataset_events_t.npy"), mmap_mode="r").reshape(-1)
    xy = np.load(os.path.join(seq_dir, "dataset_events_xy.npy"), mmap_mode="r")
    p = np.load(os.path.join(seq_dir, "dataset_events_p.npy"), mmap_mode="r").reshape(-1)
    return t, xy, p


def event_index_at_time(info: Dict[str, Any], t_target: float, t0: float, total_events: int) -> int:
    """Szybki dostęp O(1) do indeksu zdarzenia w dataset_events_t za pomocą tablicy index.

    Zwraca przybliżony indeks w zdarzeniach dla zadanego czasu t_target.
    """
    if "index" in info and "discretization" in info:
        disc = float(info["discretization"])
        idx_arr = info["index"]
        bin_idx = int(np.floor((t_target - t0) / disc))
        if bin_idx <= 0:
            return 0
        if bin_idx >= len(idx_arr):
            return int(idx_arr[-1])
        return int(idx_arr[bin_idx])
    return 0


def walk(obj: Any, path: str = "meta"):
    """Rekurencyjnie przechodzi po słownikach/listach (dla list tylko po 1. elemencie)."""
    yield path, obj
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from walk(v, f"{path}/{k}")
    elif isinstance(obj, (list, tuple)) and len(obj) > 0:
        yield from walk(obj[0], f"{path}[0]")


def find_series(info: Optional[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    """Listy słowników z kluczem 'ts' -> {ścieżka: lista} (IMU, klatki GT, trajektorie...)."""
    out = {}
    for name, root in (info or {}).items():
        for path, o in walk(root, name):
            if isinstance(o, list) and o and isinstance(o[0], dict) and "ts" in o[0]:
                out[path] = o
    return out


def find_dict_with(info: Optional[Dict[str, Any]], keys: List[str]) -> Optional[Dict[str, Any]]:
    """Pierwszy słownik zawierający wszystkie podane klucze (np. intrinsics)."""
    for name, root in (info or {}).items():
        for _, o in walk(root, name):
            if isinstance(o, dict) and all(k in o for k in keys):
                return o
    return None


def is_imu_series(lst: List[Dict[str, Any]]) -> bool:
    """Sprawdza, czy lista słowników reprezentuje próbki IMU."""
    return bool(lst and "angular_velocity" in lst[0] and "linear_acceleration" in lst[0])


def imu_arrays(lst: List[Dict[str, Any]]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Ekstrahuje szeregi czasowe IMU jako (ts, gyro, acc)."""
    ts = np.array([d["ts"] for d in lst], dtype=np.float64)
    gyro = np.array([[d["angular_velocity"][c] for c in "xyz"] for d in lst], dtype=np.float64)
    acc = np.array([[d["linear_acceleration"][c] for c in "xyz"] for d in lst], dtype=np.float64)
    return ts, gyro, acc


def short(x: Any, n: int = 90) -> str:
    s = repr(x)
    return s if len(s) <= n else s[: n - 3] + "..."


def describe(obj: Any, name: str = "root", depth: int = 0, max_depth: int = 3, out=print):
    """Zwięzłe drzewo struktury (dla list: długość + 1. element)."""
    pad = "  " * depth
    if isinstance(obj, dict):
        out(f"{pad}{name}: dict[{len(obj)}] klucze={list(obj)[:14]}{'...' if len(obj) > 14 else ''}")
        if depth < max_depth:
            for k, v in obj.items():
                describe(v, str(k), depth + 1, max_depth, out)
    elif isinstance(obj, (list, tuple)):
        extra = ""
        if obj and isinstance(obj[0], dict) and "ts" in obj[0]:
            try:
                extra = f"  ts: {obj[0]['ts']:.6f} .. {obj[-1]['ts']:.6f}"
            except Exception:
                pass
        out(f"{pad}{name}: list[{len(obj)}]{extra}")
        if obj and depth < max_depth:
            describe(obj[0], "[0]", depth + 1, max_depth, out)
    elif isinstance(obj, np.ndarray):
        s = f"{pad}{name}: ndarray shape={obj.shape} dtype={obj.dtype}"
        if obj.size and obj.dtype.kind in "iufb" and obj.size <= 50_000_000:
            s += f" min={obj.min()} max={obj.max()}"
        out(s)
    else:
        out(f"{pad}{name}: {type(obj).__name__} = {short(obj)}")


# ---------- lokalizacja wyników (pliki json, png, maski) ----------
def seq_label(name: str) -> str:
    """scene13_dyn_test_02_000000 -> scene13_02 (inne nazwy bez zmian)."""
    m = re.match(r"^(scene\d+)_dyn_(?:test|train)_(\d+)_\d+$", name)
    return f"{m.group(1)}_{m.group(2)}" if m else name


def default_out_root() -> Path:
    """Katalog główny wyników: zmienna EVIMO2_OUT albo results/sequences."""
    if "EVIMO2_OUT" in os.environ:
        return Path(os.environ["EVIMO2_OUT"])
    for cand in (
        "results/sequences",
        "../results/sequences",
        "../../results/sequences",
        "results",
        "../results",
        "out/sequences",
        "out",
    ):
        p = Path(cand)
        if p.is_dir():
            return p
    return Path("results/sequences")


def result_dir_candidates(seq_name: str, out_root: Optional[Union[str, Path]] = None) -> List[Path]:
    root = Path(out_root) if out_root else default_out_root()
    return [root / seq_label(seq_name), root / seq_name]


def result_dir(seq_name: str, out_root: Optional[Union[str, Path]] = None) -> Path:
    """Istniejący katalog wyników sekwencji (krótka nazwa ma pierwszeństwo) albo nowy."""
    cands = result_dir_candidates(seq_name, out_root)
    return next((c for c in cands if c.is_dir()), cands[0])


def find_result_file(seq_name: str, filename: str, out_root: Optional[Union[str, Path]] = None) -> Optional[Path]:
    """Pierwszy istniejący plik <katalog wyników>/<filename> albo None."""
    for c in result_dir_candidates(seq_name, out_root):
        cand_file = c / filename
        if cand_file.is_file():
            return cand_file
    return None


def lag_const_loo(out_root: Optional[Union[str, Path]], seq_name: str, fallback: float) -> Tuple[float, int]:
    """Stałe opóźnienie IMU: mediana d_GT z INNYCH sekwencji o residuum < 15% (leave-one-out)."""
    root = Path(out_root) if out_root else default_out_root()
    lags = []
    for f in sorted(glob.glob(str(root / "*" / "gyro_vs_gt.json"))):
        try:
            with open(f) as fp:
                j = json.load(fp)
            if j.get("sequence") != seq_name and j.get("residual_rel", 1.0) < 0.15:
                lags.append(float(j["lag_s"]))
        except Exception:
            continue
    return (float(np.median(lags)), len(lags)) if len(lags) >= 3 else (fallback, 0)
