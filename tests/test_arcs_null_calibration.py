"""Unit tests for scripts/arcs_null_calibration.py (frozen outcome rule and likelihood helpers)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"


def _load(name: str):
    sys.path.insert(0, str(_DIR))
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


nc = _load("arcs_null_calibration")


@pytest.mark.parametrize(
    "best, upper, want",
    [
        (0.05, 0.9, "calibration explains the pooled excess"),
        (0.2, 0.01, "calibration explains the pooled excess"),
        (0.049, 0.049, "null calibrated; excess stands"),
        (0.01, 0.05, "weak rows lack power"),
    ],
)
def test_outcome_rule(best, upper, want):
    assert nc.outcome(best, upper) == want


def test_with_iso_mixes_toward_one_third():
    p = np.array([0.0, 0.3, 1.0])
    assert np.allclose(nc.with_iso(p, 0.0), p)
    assert np.allclose(nc.with_iso(p, 1.0), 1 / 3)
    assert np.allclose(nc.with_iso(p, 0.5), 0.5 * p + 1 / 6)


def test_neg2lnl_bernoulli():
    p = np.array([0.2, 0.5])
    a = np.array([True, False])
    assert nc.neg2lnl(p, a) == pytest.approx(-2 * (np.log(0.2) + np.log(0.5)))


def test_profile_interval_and_edges():
    grid = np.arange(0.0, 1.01, 0.1)
    curve = (grid - 0.3) ** 2 * 100  # Δ ≤ 3.84 ⇔ |x − 0.3| ≤ 0.196
    r = nc.profile(grid, curve)
    assert (r["best"], r["lo95"], r["hi95"]) == pytest.approx((0.3, 0.2, 0.4))
    assert not r["lo_at_edge"] and not r["hi_at_edge"]
    r = nc.profile(grid, (grid * 10) ** 2 / 10)
    assert r["best"] == 0.0 and r["lo_at_edge"]


def test_two_sided_p_capped():
    assert nc.two_sided_p(9, 9.0) == 1.0
    assert nc.two_sided_p(0, 10.0) == pytest.approx(2 * np.exp(-10.0))


def test_rows_split_by_shear_keeps_background(tmp_path):
    from astropy.table import Table

    t = Table(
        {
            "g_min": [0.1, 0.5, 0.1, 0.5],
            "strong_shear": [False, True, False, True],
            "population": ["background", "background", "unknown", "member_or_foreground"],
            "orientation_class": ["anti", "aligned", "anti", "anti"],
        }
    )
    path = tmp_path / "arcs.ecsv"
    t.write(path)
    assert list(nc.rows(path, strong=False)["g_min"]) == [0.1]
    assert list(nc.rows(path, strong=True)["g_min"]) == [0.5]
