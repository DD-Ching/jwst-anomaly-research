"""Tests for the pure helpers of scripts/epoch_compare.py."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table

_PATH = Path(__file__).resolve().parents[1] / "scripts" / "epoch_compare.py"
_spec = importlib.util.spec_from_file_location("epoch_compare", _PATH)
ec = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ec)

RA0, DEC0 = 110.0, -73.0


def _cats(n=80, shift_mas=(10.0, -5.0), target_motion_mas=(0.0, 0.0), dmag=0.1, target_dmag=0.0):
    rng = np.random.default_rng(3)
    cosd = np.cos(np.deg2rad(DEC0))
    dra = rng.uniform(-40, 40, n) / 3600 / cosd
    ddec = rng.uniform(-40, 40, n) / 3600
    dra[0] = ddec[0] = 0.0  # the target
    mag = rng.uniform(22, 25, n)
    cat1 = Table(
        {
            "ra": RA0 + dra,
            "dec": DEC0 + ddec,
            "aper50_abmag": mag,
            "aper50_abmag_err": np.full(n, 0.01),
        }
    )
    noise = rng.normal(0, 3.0, (n, 2))  # mas, galaxy centroid scatter
    move = np.array(shift_mas) + noise
    move[0] = np.array(shift_mas) + np.array(target_motion_mas)
    cat2 = Table(
        {
            "ra": cat1["ra"] + move[:, 0] / 3.6e6 / cosd,
            "dec": cat1["dec"] + move[:, 1] / 3.6e6,
            "aper50_abmag": mag + dmag + rng.normal(0, 0.02, n),
            "aper50_abmag_err": np.full(n, 0.01),
        }
    )
    cat2["aper50_abmag"][0] = mag[0] + dmag + target_dmag
    return cat1, cat2, SkyCoord(RA0, DEC0, unit="deg")


def test_static_target_has_no_residual_after_the_frame_tie():
    cat1, cat2, target = _cats()
    r = ec.compare_epochs(cat1, cat2, target)
    assert r["n_references"] == 79
    np.testing.assert_allclose(r["frame_offset_mas"], [10.0, -5.0], atol=1.5)
    assert r["target_residual_total_mas"] < 2.0
    assert abs(r["target_dmag"]) < 0.02
    assert r["label"] == "derived"


def test_moving_and_variable_target_is_detected():
    cat1, cat2, target = _cats(target_motion_mas=(30.0, 0.0), target_dmag=0.5)
    r = ec.compare_epochs(cat1, cat2, target)
    assert r["target_residual_mas"][0] == pytest.approx(30.0, abs=1.5)
    assert r["target_residual_over_tie_error"] > 20
    assert r["target_dmag"] == pytest.approx(0.5, abs=0.02)


def test_references_with_missing_magnitudes_are_skipped():
    cat1, cat2, target = _cats()
    cat1["aper50_abmag"][5] = np.nan
    cat2["aper50_abmag"][6] = np.nan
    r = ec.compare_epochs(cat1, cat2, target)
    assert r["n_references"] == 77 and np.isfinite(r["target_dmag"])


def test_missing_target_or_too_few_references_raise():
    cat1, cat2, _ = _cats()
    far = SkyCoord(RA0 + 0.1, DEC0, unit="deg")
    with pytest.raises(ValueError, match="target not within"):
        ec.compare_epochs(cat1, cat2, far)
    cat1, cat2, target = _cats(n=10)
    with pytest.raises(ValueError, match="reference sources"):
        ec.compare_epochs(cat1, cat2, target)


def test_mutual_matches_drop_contested_pairs():
    a = SkyCoord([10.0, 10.0 + 0.1 / 3600], [0.0, 0.0], unit="deg")
    b = SkyCoord([10.0 + 0.05 / 3600], [0.0], unit="deg")
    ia, ib = ec.mutual_matches(a, b, 0.3)
    assert len(ia) == 1 and len(ib) == 1
