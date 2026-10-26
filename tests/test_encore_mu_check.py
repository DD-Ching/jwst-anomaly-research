"""Rules of scripts/encore_mu_check.py (frozen plan results/lens_consistency/encore_mu/PLAN.md)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import encore_mu_check as e  # noqa: E402


def test_ln_sigma_swaps_errors_for_negative_median():
    assert math.isclose(e.ln_sigma((-20.0, 2.0, 4.0)), e.ln_sigma((20.0, 4.0, 2.0)))


def test_observed_values_give_zero_z_and_offset():
    obs = {"1a": (-e.OBS["1a"][0], 1.0, 1.0), "1b": e.OBS["1b"], "dt": e.OBS["dt"]}
    assert math.isclose(e.offset_mag(obs), 0.0, abs_tol=1e-12)
    assert all(math.isclose(e.image_z(obs, i), 0.0, abs_tol=1e-12) for i in e.IMAGES)
    assert math.isclose(e.ratio_z(obs), 0.0, abs_tol=1e-12)


def test_mass_sheet_recovers_injected_lambda():
    lam = 1.3
    model = {i: (e.OBS[i][0] * lam**2, 0.1, 0.1) for i in e.IMAGES}
    model["dt"] = (e.OBS["dt"][0] / lam, 1.0, 1.0)
    r = e.h0_mst(model)
    assert math.isclose(r["lambda_mu"], lam)
    assert math.isclose(r["lambda_t"], lam)
    assert math.isclose(r["h0_mst"], e.H0_FID)


def test_outcome_needs_every_model():
    rows = {"a": {"image_z": {"1a": 4.0, "1b": 0.0}, "ratio_z": 4.0}}
    rows["b"] = {"image_z": {"1a": 3.5, "1b": 0.0}, "ratio_z": 1.0}
    assert e.outcome(rows) == (["1a"], [], "F")
    rows["b"]["image_z"]["1a"] = 2.9
    assert e.outcome(rows) == ([], [], "N")


def test_same_sense_rule():
    assert e.same_sense([0.1] * 5 + [-0.05] * 2)
    assert not e.same_sense([0.1] * 4 + [-0.01] * 3)
    assert not e.same_sense([0.01] * 5 + [-1.0] * 2)
