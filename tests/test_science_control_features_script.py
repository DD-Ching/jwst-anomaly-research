import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import science_control_features as s  # noqa: E402


def test_verdict_follows_plan_thresholds():
    assert s.verdict(0.05) == "comparable"
    assert s.verdict(0.1) == "comparable"
    assert s.verdict(0.15) == "intermediate"
    assert s.verdict(0.2) == "intermediate"
    assert s.verdict(0.21) == "differ strongly"
    assert s.verdict(float("nan")) == "not judged"


def test_compare_small_bins_not_judged_and_shift_detected():
    rng = np.random.default_rng(0)
    small = s.compare(rng.normal(size=10), rng.normal(size=500))
    assert np.isnan(small["ks_d"])
    same = s.compare(rng.normal(size=2000), rng.normal(size=2000))
    shifted = s.compare(rng.normal(size=2000), rng.normal(1.0, size=2000))
    assert same["ks_d"] < s.D_OK < s.D_STRONG < shifted["ks_d"]
