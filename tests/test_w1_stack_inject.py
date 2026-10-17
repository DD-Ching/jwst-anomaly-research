import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import w1_stack_inject as wsi  # noqa: E402


def test_excluded_needs_level_share_above_observed():
    st = np.r_[np.full(95, 1.0), np.full(5, 0.0)]
    assert wsi.excluded(st, 0.5)
    assert not wsi.excluded(np.r_[np.full(94, 1.0), np.full(6, 0.0)], 0.5)
    assert not wsi.excluded(np.full(100, 0.5), 0.5)  # equal is not above


def test_efficiency():
    assert wsi.efficiency(np.array([0.1, 0.5, 0.9, 1.0]), 0.45) == 0.75


def test_fraction_stacks():
    s_obs = np.zeros(10)
    s_inj = np.ones(10)
    rng = np.random.default_rng(0)
    assert np.allclose(
        wsi.fraction_stacks(s_obs, s_inj, 0.25, 50, rng), 0.2
    )  # round(2.5) = 2 (half to even)
    assert np.allclose(wsi.fraction_stacks(s_obs, s_inj, 0.01, 5, rng), 0.1)  # at least one
    one = wsi.fraction_stacks(s_obs, s_inj, 1.0, 50, rng)
    assert one.shape == (1,) and one[0] == 1.0


def test_noise_stacks_adds_delta_to_each_null_draw():
    s_null = np.array([[0.0, 1.0], [0.0, -1.0]])
    rng = np.random.default_rng(0)
    out = wsi.noise_stacks(s_null, np.array([2.0, 4.0]), 1.0, rng)
    assert np.allclose(out, [3.0, 3.0])
    out = wsi.noise_stacks(s_null, np.array([np.nan, 4.0]), 1.0, rng)  # NaN centre drops out
    assert np.allclose(out, [4.0, 3.0])
    one = wsi.noise_stacks(s_null, np.array([2.0, 2.0]), 0.5, rng)
    assert np.allclose(one, [1.0, 1.0])


def test_fraction_stacks_skips_undefined_centres():
    out = wsi.fraction_stacks(
        np.zeros(2), np.array([np.nan, 2.0]), 1.0, 5, np.random.default_rng(0)
    )
    assert out[0] == 2.0
