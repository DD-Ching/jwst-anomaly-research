import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import smacs_zphot_nofit as s  # noqa: E402

FIT_NUSE = np.arange(5, 25).astype(float)  # 5th percentile 5.95
FIT_R = np.linspace(2.0, 10.0, 101)  # 10th percentile 2.8


def test_coverage_wins_when_most_nofit_rows_lack_filters():
    nofit = np.array([2, 3, 3, 4, 20], float)  # 4 / 5 = 0.8 below 5.95
    assert s.decide(nofit, FIT_NUSE, np.full(5, 1.0), FIT_R, 4) == "H-cov: coverage flag"


def test_star_needs_small_median_and_compact_top_rows():
    nofit = np.full(5, 20.0)
    small = np.full(5, 2.0)
    assert s.decide(nofit, FIT_NUSE, small, FIT_R, 2) == "H-star supported"
    assert s.decide(nofit, FIT_NUSE, small, FIT_R, 1) == "H-mix / undetermined"
    assert s.decide(nofit, FIT_NUSE, small, FIT_R, None) == "H-mix / undetermined"
    assert s.decide(nofit, FIT_NUSE, np.full(5, 5.0), FIT_R, 4) == "H-mix / undetermined"


def test_coverage_below_threshold_is_not_coverage():
    nofit = np.array([2, 3, 3, 20, 20], float)  # 0.6
    assert s.decide(nofit, FIT_NUSE, np.full(5, 5.0), FIT_R, 0) == "H-mix / undetermined"


def test_by_mag_bins_groups():
    mag = np.array([23.0, 23.5, 27.5, 30.0])
    rad = np.array([1.0, 3.0, 2.0, 5.0])
    g = {"a": np.array([True, True, False, True]), "b": np.array([False, False, True, False])}
    rows = s.by_mag(mag, rad, g)
    assert rows[0]["a"] == {"n": 2, "median": 2.0}
    assert rows[3]["b"] == {"n": 1, "median": 2.0}
    assert rows[5]["a"]["n"] == 1 and rows[5]["b"] == {"n": 0, "median": None}
