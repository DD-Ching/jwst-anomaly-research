"""Tests for scripts/w1_peak_vet.py (synthetic inputs only)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w1_peak_vet as pv  # noqa: E402


def test_stripe_columns_finds_bright_band():
    rng = np.random.default_rng(0)
    img = rng.normal(0, 1, (200, 120))
    img[:, 50:55] += 3.0
    assert list(pv.stripe_columns(img)) == [50, 51, 52, 53, 54]


def test_stripe_columns_none_in_pure_noise():
    img = np.random.default_rng(1).normal(0, 1, (400, 100))
    assert len(pv.stripe_columns(img)) == 0


def test_axis_points_at_wraps_mod_180():
    assert pv.axis_points_at(np.array([5.0, 175.0, 95.0]), 0.0).tolist() == [True, True, False]
    assert pv.axis_points_at(np.array([190.0]), 10.0).tolist() == [True]


def test_stripe_columns_ignores_fill_and_zero_mad():
    img = np.zeros((50, 40))
    assert len(pv.stripe_columns(img)) == 0
    rng = np.random.default_rng(2)
    img = rng.normal(5, 1, (300, 80))
    img[:, 60:] = 0.0  # off-mosaic fill
    img[:, 20:23] += 4.0
    assert list(pv.stripe_columns(img)) == [20, 21, 22]


def test_outcome_rule():
    assert pv.outcome(0.1, 1.0, 3.6).startswith("explained")
    assert pv.outcome(3.7, 0.06, 3.6).startswith("explained")
    assert pv.outcome(3.6, 0.01, 3.6).startswith("explained")
    assert pv.outcome(3.9, 0.01, 3.6).startswith("not explained")
