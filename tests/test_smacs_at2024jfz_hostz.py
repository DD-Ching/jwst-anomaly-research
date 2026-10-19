import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import smacs_at2024jfz_hostz as hz  # noqa: E402


def test_classify_windows():
    assert hz.classify(None, None) == "untested"
    assert hz.classify(0.37, 0.40) == "member_or_fg"
    assert hz.classify(0.30, 0.50) == "untested"
    assert hz.classify(1.36, 1.48) == "tested_window"
    assert hz.classify(1.80, 2.40) == "tested_window"
    assert hz.classify(1.40, 1.60) == "undetermined_window"
    assert hz.classify(1.60, 1.90) == "undetermined_window"
    assert hz.classify(1.50, 1.54) == "untested"
    assert hz.classify(2.60, 3.10) == "untested"
    assert hz.classify(np.nan, 1.0) == "untested"


def test_pick_host_normalised_distance():
    sep = np.array([0.3, 0.6, 1.2])
    rh = np.array([0.05, 0.5, 2.0])
    # 0.3 / 0.1 (floor) = 3.0 vs 0.6 / 0.5 = 1.2; the third is outside R_MAX
    assert hz.pick_host(sep, rh) == 1
    assert hz.pick_host(np.array([1.5]), np.array([1.0])) is None
