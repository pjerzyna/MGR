"""Moduł geometrii epipolarnej i weryfikacji Focus of Expansion (FOE) dla detektora IMO v2.

Eliminuje fałszywe alarmy pochodzące od paralaksy translacyjnej po uprzedniej kompensacji rotacji SO(3).
"""
from typing import Dict, Tuple, Optional
import numpy as np

from .detector import box_blur


def epipolar_direction(
    u: float,
    v: float,
    v_cam: np.ndarray,
    intr: Dict[str, float]
) -> Tuple[np.ndarray, float]:
    """Wyznacza teoretyczny kierunek wektora pozornego ruchu tła na matrycy w punkcie (u, v).
    
    Wskutek translacji kamery z prędkością liniową v_cam = (vx, vy, vz) w układzie kamery,
    dowolny punkt nieruchomego tła o głębi Z > 0 porusza się zgodnie z wektorem:
        u_dot = fx * (-vx + x * vz) / Z
        v_dot = fy * (-vy + y * vz) / Z
        gdzie x = (u - cx)/fx, y = (v - cy)/fy.
        
    Zwraca:
        (u_epi_unit, epi_norm): znormalizowany wektor 2D [e_u, e_v] oraz jego długość bez uwzględnienia Z.
    """
    fx = float(intr["fx"])
    fy = float(intr["fy"])
    cx = float(intr["cx"])
    cy = float(intr["cy"])
    
    x_norm = (u - cx) / fx
    y_norm = (v - cy) / fy
    
    # Pozorny ruch krawędzi tła (współrzędne pikselowe na sekundę, pomnożone przez Z):
    u_dot = fx * (-v_cam[0] + x_norm * v_cam[2])
    v_dot = fy * (-v_cam[1] + y_norm * v_cam[2])
    
    epi_norm = float(np.hypot(u_dot, v_dot))
    if epi_norm < 1e-7:
        return np.zeros(2, dtype=np.float32), 0.0
    
    return np.array([u_dot / epi_norm, v_dot / epi_norm], dtype=np.float32), epi_norm


def estimate_cluster_shift(
    Es: np.ndarray,
    Ls: np.ndarray,
    x0: int,
    x1: int,
    y0: int,
    y1: int,
    max_shift: int = 8
) -> Tuple[Tuple[int, int], float]:
    """Wyznacza dyskretne przesunięcie (dx, dy) między wczesnym obrazem Es a późnym Ls w wycinku klastra.
    
    Maksymalizuje znormalizowaną korelację wzajemną (NCC):
        corr = sum(E * L_shifted) / (||E|| * ||L_shifted||)
        
    Zwraca:
        ((best_dx, best_dy), best_corr)
    """
    H, W = Es.shape
    crop_E = Es[y0:y1, x0:x1]
    norm_E = float(np.linalg.norm(crop_E))
    if norm_E < 1e-4:
        return (0, 0), 0.0
        
    best_corr = -1.0
    best_shift = (0, 0)
    
    for dy in range(-max_shift, max_shift + 1):
        y0_l = y0 + dy
        y1_l = y1 + dy
        if y0_l < 0 or y1_l > H:
            continue
        for dx in range(-max_shift, max_shift + 1):
            x0_l = x0 + dx
            x1_l = x1 + dx
            if x0_l < 0 or x1_l > W:
                continue
                
            crop_L = Ls[y0_l:y1_l, x0_l:x1_l]
            norm_L = float(np.linalg.norm(crop_L))
            if norm_L < 1e-4:
                continue
                
            corr = float(np.sum(crop_E * crop_L) / (norm_E * norm_L))
            if corr > best_corr:
                best_corr = corr
                best_shift = (dx, dy)
                
    return best_shift, best_corr


def filter_clusters_epipolar(
    clusters: list,
    Es: np.ndarray,
    Ls: np.ndarray,
    v_cam: np.ndarray,
    intr: Dict[str, float],
    min_speed: float = 0.05,
    angle_thresh_deg: float = 35.0,
    max_shift: int = 8
) -> list:
    """Weryfikuje klastry kandydatów pod kątem zgodności z geometrią epipolarną.
    
    Dla każdego klastra sprawdza, czy jego wektor przemieszczenia leży w stożku
    linii epipolarnych rozchodzących się od punktu FOE.
    
    Zwraca listę klastrów zaklasyfikowanych jako prawdziwe IMO (odrzucając paralaksę tła).
    """
    speed = float(np.linalg.norm(v_cam))
    kept_clusters = []
    
    for c in clusters:
        x0, x1, y0, y1 = c["bbox"]
        shift, corr = estimate_cluster_shift(Es, Ls, x0, x1, y0, y1, max_shift=max_shift)
        dx, dy = shift
        shift_len = float(np.hypot(dx, dy))
        
        c["shift"] = (dx, dy)
        c["shift_len"] = shift_len
        c["corr"] = corr
        
        xm = 0.5 * (x0 + x1)
        ym = 0.5 * (y0 + y1)
        u_epi, epi_norm = epipolar_direction(xm, ym, v_cam, intr)
        
        if epi_norm > 1e-4 and shift_len > 0 and speed > min_speed:
            u_s = np.array([dx, dy], dtype=np.float32) / shift_len
            cos_th = float(np.clip(np.dot(u_s, u_epi), -1.0, 1.0))
            ang_deg = float(np.degrees(np.arccos(cos_th)))
            c["angle_deg"] = ang_deg
            
            # Warunek IMO:
            # Ruch odchyla się od linii epipolarnej o co najmniej angle_thresh_deg LUB
            # jest skierowany przeciwnie do wektora epipolarnego (cos_th < 0)
            is_imo = (ang_deg >= angle_thresh_deg)
        else:
            # Kamera nie porusza się translacyjnie lub klaster nie wykazuje ruchu dyskretnego
            c["angle_deg"] = None
            is_imo = True  # Pozostawiamy decyzję z detektora rotacji v1
            
        c["is_imo"] = is_imo
        if is_imo:
            kept_clusters.append(c)
            
    return kept_clusters
