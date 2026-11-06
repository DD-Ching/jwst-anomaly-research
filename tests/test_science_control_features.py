"""Frozen rules of scripts/science_control_features.py (verdict thresholds, shared features)."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

_SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location(
    "science_control_features", _SCRIPTS / "science_control_features.py"
)
scf = importlib.util.module_from_spec(_spec)
sys.modules["science_control_features"] = scf
_spec.loader.exec_module(scf)


def test_verdict_thresholds():
    assert scf.verdict([]) == "untested"
    assert scf.verdict([0.05, 0.1]) == "comparable"
    assert scf.verdict([0.05, 0.15]) == "intermediate"
    assert scf.verdict([0.2]) == "intermediate"
    assert scf.verdict([0.05, 0.21]) == "differ strongly"


def test_compare_identical_and_shifted():
    rng = np.random.default_rng(0)
    x = rng.normal(size=500)
    same = scf.compare(x, x)
    assert same["ks_d"] == 0.0
    shifted = scf.compare(x, x + 1.0)
    assert shifted["ks_d"] > scf.D_STRONG
    small = scf.compare(x[:5], x[:5])
    assert small["ks_d"] is None and small["n_science"] == 5


def test_shared_features_skip_ids_and_text():
    a = Table({"source_uid": ["a"], "ra": [1.0], "ref_mag": [20.0], "color_f090w_f150w": [0.1]})
    b = Table({"source_uid": ["b"], "ra": [1.0], "ref_mag": [21.0], "color_f115w_f150w": [0.2]})
    assert scf.shared_features(a, b) == ["ref_mag"]


def test_compact_blue_mask_needs_both_and_finite():
    t = Table(
        {
            "ref_log_ci_50_30": [0.20, 0.20, 0.30, np.nan],
            "color_f277w_f356w": [-0.5, 0.0, -0.5, -0.5],
        }
    )
    assert scf.compact_blue_mask(t).tolist() == [True, False, False, False]
