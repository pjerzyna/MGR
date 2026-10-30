"""Wspólne funkcje pomocnicze dla narzędzi EVIMO2 (wczytywanie, przeszukiwanie meta)."""
import os
import re
from pathlib import Path

import numpy as np


def unwrap(a):
    """0-wymiarowe tablice obiektowe (np. słownik meta) -> zwykły obiekt Pythona."""
    if isinstance(a, np.ndarray) and a.shape == () and a.dtype == object:
        return a.item()
    return a


def load_info(seq_dir):
    """Wczytuje dataset_info.npz do zwykłego słownika {klucz: obiekt}."""
    path = os.path.join(seq_dir, "dataset_info.npz")
    if not os.path.isfile(path):
        return None
    with np.load(path, allow_pickle=True) as z:
        return {k: unwrap(z[k]) for k in z.files}


def load_events(seq_dir):
    """Zwraca (t, xy, p) jako widoki mmap o kształtach (N,), (N,2), (N,)."""
    t = np.load(os.path.join(seq_dir, "dataset_events_t.npy"), mmap_mode="r").reshape(-1)
    xy = np.load(os.path.join(seq_dir, "dataset_events_xy.npy"), mmap_mode="r")
    p = np.load(os.path.join(seq_dir, "dataset_events_p.npy"), mmap_mode="r").reshape(-1)
    return t, xy, p


def walk(obj, path="meta"):
    """Rekurencyjnie przechodzi po słownikach/listach (dla list tylko po 1. elemencie)."""
    yield path, obj
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from walk(v, f"{path}/{k}")
    elif isinstance(obj, (list, tuple)) and len(obj) > 0:
        yield from walk(obj[0], f"{path}[0]")


def find_series(info):
    """Listy słowników z kluczem 'ts' -> {ścieżka: lista} (IMU, klatki GT, trajektorie...)."""
    out = {}
    for name, root in (info or {}).items():
        for path, o in walk(root, name):
            if isinstance(o, list) and o and isinstance(o[0], dict) and "ts" in o[0]:
                out[path] = o
    return out


def find_dict_with(info, keys):
    """Pierwszy słownik zawierający wszystkie podane klucze (np. intrinsics)."""
    for name, root in (info or {}).items():
        for _, o in walk(root, name):
            if isinstance(o, dict) and all(k in o for k in keys):
                return o
    return None


def is_imu_series(lst):
    return "angular_velocity" in lst[0] and "linear_acceleration" in lst[0]


def imu_arrays(lst):
    ts = np.array([d["ts"] for d in lst], dtype=np.float64)
    gyro = np.array([[d["angular_velocity"][c] for c in "xyz"] for d in lst])
    acc = np.array([[d["linear_acceleration"][c] for c in "xyz"] for d in lst])
    return ts, gyro, acc


def short(x, n=90):
    s = repr(x)
    return s if len(s) <= n else s[: n - 3] + "..."


def describe(obj, name="root", depth=0, max_depth=3, out=print):
    """Zwięzłe drzewo struktury (listy: tylko długość + 1. element)."""
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


# ---------- lokalizacja wyników (gyro_vs_gt.json, event_lag_check.json, ...) ----------
def seq_label(name):
    """scene13_dyn_test_02_000000 -> scene13_02 (inne nazwy bez zmian)."""
    m = re.match(r"^(scene\d+)_dyn_(?:test|train)_(\d+)_\d+$", name)
    return f"{m.group(1)}_{m.group(2)}" if m else name


def default_out_root():
    """Katalog główny wyników: zmienna środowiskowa EVIMO2_OUT albo ./out."""
    return Path(os.environ.get("EVIMO2_OUT", "out"))


def result_dir_candidates(seq_name, out_root=None):
    root = Path(out_root) if out_root else default_out_root()
    return [root / seq_label(seq_name), root / seq_name]


def result_dir(seq_name, out_root=None):
    """Istniejący katalog wyników sekwencji (krótka nazwa ma pierwszeństwo) albo nowy (krótka nazwa)."""
    cands = result_dir_candidates(seq_name, out_root)
    return next((c for c in cands if c.is_dir()), cands[0])


def find_result_file(seq_name, filename, out_root=None):
    """Pierwszy istniejący plik <katalog wyników>/<filename> albo None."""
    for c in result_dir_candidates(seq_name, out_root):
        if (c / filename).is_file():
            return c / filename
    return None