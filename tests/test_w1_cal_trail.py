import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w1_cal_trail as t  # noqa: E402


def test_classify():
    assert t.classify([]).startswith("undecided")
    assert t.classify([True, True]).startswith("in every exposure")
    assert t.classify([False, False]).startswith("in no exposure")
    assert t.classify([True, False, False]).startswith("in 1 of 3")


def test_contrast_and_track_on_synthetic_line():
    rng = np.random.default_rng(1)
    stamp = rng.normal(0.0, 1.0, (201, 600))
    rows = np.round(100 + 0.01 * (np.arange(600) - 300)).astype(int)
    stamp[rows, np.arange(600)] += 30.0
    prof = t.row_profile(stamp)
    c = t.contrast(prof, 100, scale=0.03)
    assert c["contrast"] > t.IN_SIGMA
    tr = t.line_track(stamp, 100)
    assert abs(tr["slope_px_per_px"] - 0.01) < 0.003
    flat = t.contrast(t.row_profile(rng.normal(0.0, 1.0, (201, 600))), 100, scale=0.03)
    assert flat["contrast"] < t.IN_SIGMA
