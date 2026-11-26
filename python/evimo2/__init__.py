"""EVIMO2 Processing Pipeline: Fuzja danych z kamer zdarzeniowych i IMU do detekcji IMO."""

from .calibration import GYRO_BIAS, LAG_CONST_S, R_IMU2CAM, calibrated_gyro
from .common import event_index_at_time, load_events, load_info
from .detector import depth_for_events, detect_imo, residual_map, regions_for_frame
from .epipolar import epipolar_direction, estimate_cluster_shift, filter_clusters_epipolar
from .geometry import UndistortLUT, nlerp, rodrigues_batch, undistort
from .kinematics import CLASSES, kinematics, world_pose, load_tracks, pose_series
from .rotcomp import (
    AdaptiveWindowConfig,
    generate_adaptive_windows,
    generate_fixed_windows,
    gt_rot_matrices,
    gt_warp,
    iwe,
    prep_window,
    rot_matrices,
    warp,
    warp_translation,
)

__version__ = "0.1.0"
__all__ = [
    "load_info",
    "load_events",
    "event_index_at_time",
    "R_IMU2CAM",
    "GYRO_BIAS",
    "LAG_CONST_S",
    "calibrated_gyro",
    "undistort",
    "UndistortLUT",
    "rodrigues_batch",
    "nlerp",
    "prep_window",
    "AdaptiveWindowConfig",
    "generate_adaptive_windows",
    "generate_fixed_windows",
    "rot_matrices",
    "warp",
    "iwe",
    "gt_rot_matrices",
    "gt_warp",
    "warp_translation",
    "CLASSES",
    "world_pose",
    "load_tracks",
    "pose_series",
    "kinematics",
    "detect_imo",
    "residual_map",
    "regions_for_frame",
    "depth_for_events",
    "epipolar_direction",
    "estimate_cluster_shift",
    "filter_clusters_epipolar",
]
