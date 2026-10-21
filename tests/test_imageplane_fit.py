"""Offline tests for scripts/imageplane_fit.py (cost, starts and gate helpers; no network)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import imageplane_fit as f  # noqa: E402


def test_match_cost_pairs_one_to_one():
    cost, d, paired = f.match_cost([0, 10], [0, 0], [10.1, 0.2], [0, 0], penalty=10.0)
    assert paired.all()
    np.testing.assert_allclose(d, [0.2, 0.1])
    np.testing.assert_allclose(cost, 0.05)


def test_match_cost_penalises_missing_prediction():
    cost, d, paired = f.match_cost([0, 1], [0, 0], [0.5], [0], penalty=10.0)
    assert paired.sum() == 1
    np.testing.assert_allclose(d, [0.5, 0.5])
    np.testing.assert_allclose(cost, 0.5 + 100.0)


def test_match_cost_without_predictions():
    cost, d, paired = f.match_cost([0, 1], [0, 0], [], [], penalty=10.0)
    assert not paired.any() and np.isnan(d).all()
    assert cost > 2 * 100.0


def test_starts_merge_close_traced_sources():
    s = f.starts([0.0, 0.01, 2.0], [0.0, 0.0, 0.0], merge=0.05)
    assert s[0] == (np.mean([0.0, 0.01, 2.0]), 0.0)
    assert (2.0, 0.0) in s and (0.0, 0.0) in s and (0.01, 0.0) not in s


def test_v1_gate():
    assert f.v1_gate([0.01] * 9 + [np.nan], 0.05, 9)["passes"]
    assert not f.v1_gate([0.01] * 8 + [0.06, np.nan], 0.05, 9)["passes"]


def test_rms_ignores_nan():
    assert f.rms([3.0, 4.0, np.nan]) == np.sqrt(12.5)
    assert np.isnan(f.rms([np.nan]))
