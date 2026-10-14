"""Offline tests for scripts/m1149_sys16.py (pure helpers; no network)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import m1149_sys16 as s16  # noqa: E402

RA0, DEC0 = 177.4097, 22.3877
ARC = 1 / 3600.0


def _t(offsets):
    cd = np.cos(np.radians(DEC0))
    return Table(
        {
            "ra": [RA0 + dx * ARC / cd for dx, _ in offsets],
            "dec": [DEC0 + dy * ARC for _, dy in offsets],
            "magnification": [10.0] * len(offsets),
        }
    )


def test_assign_pair_distinct_and_rms():
    cat = _t([(0.0, 0.0), (0.6, 0.6)])
    pred = _t([(0.62, 0.6), (8.0, -6.0), (0.1, 0.0)])  # reordered, plus a far image
    idx, rms = s16.assign_pair(pred, cat)
    assert idx == [2, 0]
    assert rms == pytest.approx(np.sqrt((0.1**2 + 0.02**2) / 2), abs=1e-3)


def test_assign_pair_never_reuses_one_image():
    cat = _t([(0.0, 0.0), (0.6, 0.6)])
    pred = _t([(0.3, 0.3), (9.0, 9.0)])  # one image between both: the other must be the far one
    idx, rms = s16.assign_pair(pred, cat)
    assert len(set(idx)) == 2 and rms > 5
    assert s16.assign_pair(_t([(0.0, 0.0)]), cat) == ([], float("inf"))


def test_predicted_flux_and_colours():
    assert s16.predicted_flux(2.0, (30.0, -10.0), -4.0) == pytest.approx(0.2)
    assert np.isnan(s16.predicted_flux(2.0, (np.nan, 1.0), 4.0))
    c = s16.colours({"f150w": 1.0, "f277w": 10 ** (0.4), "f444w": 10 ** (0.4)})
    assert c == pytest.approx((1.0, 0.0))
    assert np.isnan(s16.colours({"f150w": -1.0, "f277w": 1.0, "f444w": 1.0})[0])


def test_z_verdict():
    assert s16.z_verdict(0.2, ["recovered", "undetectable"]) == "allowed"
    assert s16.z_verdict(0.2, ["absent"]).startswith("excluded: extra")
    assert s16.z_verdict(2.0, ["recovered"]).startswith("excluded: pair")
    assert s16.z_verdict(float("inf"), []).startswith("excluded: pair")
