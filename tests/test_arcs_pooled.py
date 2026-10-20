"""Tests for scripts/arcs_pooled.py: frozen outcome rule, pooling, post hoc S1063 comparison."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ap = _load("arcs_pooled")


@pytest.mark.parametrize(
    ("pa", "pb", "want"),
    [
        (0.01, 0.04, "pooled excess"),
        (0.05, 0.2, "no pooled excess"),
        (0.01, 0.05, "mixed"),
        (0.2, 0.01, "mixed"),
    ],
)
def test_outcome_rule(pa, pb, want):
    assert ap.outcome(pa, pb) == want


def test_pooled_sums_and_poisson():
    stats = {"a": {"anti": 2, "expected": 0.5}, "b": {"anti": 0, "expected": 0.5}}
    r = ap.pooled(stats, ("a", "b"))
    assert r["anti"] == 2 and r["expected"] == 1.0
    assert r["p_ge_observed"] == pytest.approx(1 - 2 / 2.718281828459045, rel=1e-6)


def test_variant_b_drops_only_macs1423():
    assert set(ap.VARIANTS["A"]) - set(ap.VARIANTS["B"]) == {"MACS1423"}


def test_s1063_vs_pooled_scales_and_conditions():
    r = ap.s1063_vs_pooled({"anti": 1, "expected": 1.0}, {"anti": 2, "expected": 1.0})
    assert r["pool_ratio"] == 2.0 and r["scaled_expected"] == 2.0
    assert r["p_conditional"] == pytest.approx(1 - 0.5**3)  # P(>= 1 of 3 at p = 0.5)


def test_background_keeps_strong_shear_background_only(tmp_path):
    from astropy.table import Table

    t = Table(
        {
            "strong_shear": [True, True, False, True],
            "population": ["background", "member_or_foreground", "background", "background"],
            "g_min": [0.3, 0.5, 0.1, 1.2],
            "orientation_class": ["anti", "anti", "anti", "aligned"],
        }
    )
    path = tmp_path / "arcs.ecsv"
    t.write(path)
    b = ap.background(path)
    assert list(b["g_min"]) == [0.3, 1.2]


def test_stored_check_tolerates_only_rounding():
    stored = {"n": 4, "anti": 1, "expected": 0.095}
    assert ap.stored_check({"n": 4, "anti": 1, "expected": 0.09512}, stored)["matches"]
    assert not ap.stored_check({"n": 4, "anti": 1, "expected": 0.0962}, stored)["matches"]
    assert not ap.stored_check({"n": 5, "anti": 1, "expected": 0.095}, stored)["matches"]


def test_ratio_interval_garwood():
    lo, hi = ap.ratio_interval(0, 1.0)
    assert lo == 0.0 and hi == pytest.approx(3.6889, abs=1e-3)
    lo, hi = ap.ratio_interval(6, 2.0)
    assert lo == pytest.approx(2.2019 / 2, abs=1e-3) and hi == pytest.approx(13.059 / 2, abs=1e-3)


def test_frozen_fields_point_at_stored_runs():
    for rel in (*ap.FIELDS.values(), ap.S1063):
        path = ap.ROOT / rel
        assert path.exists() and (path.parent / "arcs.json").exists()


def test_sigma_crossing_finds_root_or_none():
    assert ap.sigma_crossing(lambda s: s, 0.0, 0.1) == pytest.approx(0.05, abs=1e-3)
    assert ap.sigma_crossing(lambda s: 0.5, 0.0, 0.1) is None


def test_inside_matches_anti_probability_boundary():
    from astropy.table import Table

    t = Table({"g_min": [0.5, 0.996, 1.004, 1.2]})
    assert list(ap.inside(t)) == [False, False, False, True]


def test_field_stats_rejects_empty_rows():
    from astropy.table import Table

    t = Table({"g_min": np.array([], float), "orientation_class": np.array([], str)})
    with pytest.raises(ValueError):
        ap.field_stats(t, 0.26, 0.5, 60.0)
