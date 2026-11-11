"""Shim dla kompatybilności wstecznej - moduł evimo2_common został przeniesiony do evimo2.common."""
import sys
from pathlib import Path

_pkg_root = str(Path(__file__).resolve().parents[1])
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)

from evimo2.common import *
