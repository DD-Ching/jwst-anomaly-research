"""Shared-design leave-one-out stamp kernels equal per-stamp ``kernel_fit`` (synthetic, offline)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import arc_diff as ad  # noqa: E402


def _scene(seed=0, n_src=14, size=140, bright=0.0):
    rng = np.random.default_rng(seed)
    ref = rng.normal(0, 0.02, (size, size))
    for _ in range(n_src):
        y, x = rng.uniform(12, size - 12, 2)
        ad.add_star(ref, y, x, 1.0, rng.uniform(20, 80))
    if bright:
        ad.add_star(ref, 30.3, 100.6, 1.0, bright)  # one stamp dominating the Gram
    tgt = ndimage.gaussian_filter(ref, 0.7) * 1.03 + rng.normal(0, 0.02, ref.shape)
    tgt[70, 70] += 30.0  # an outlier that the clipping must remove
    var = np.full_like(ref, 0.02**2)
    stamp = ad.compact_stamps(ref, 3) & ad.compact_stamps(tgt, 3)
    lab, nlab = ndimage.label(stamp)
    return ref, tgt, var, stamp, lab, nlab


@pytest.mark.parametrize("bright", [0.0, 1e5, 1e6])
def test_shared_equals_refit_per_stamp(bright):
    ref, tgt, var, stamp, lab, nlab = _scene(bright=bright)
    assert nlab >= 6
    chi2: dict = {}
    shared = ad.kernel_fit_stamps(ref, tgt, var, 3, stamp, lab, nlab, chi2=chi2)
    n_cmp = 0
    for i in range(1, nlab + 1):
        try:
            k_old, c_old = ad.kernel_fit(ref, tgt, var, 3, mask=stamp & (lab != i))
        except ValueError:
            assert i not in shared
            continue
        k_new = shared[i]
        assert np.max(np.abs(k_new - k_old)) <= 1e-6 * np.max(np.abs(k_old))
        assert chi2[i] == pytest.approx(c_old, rel=1e-6)  # same clipping
        n_cmp += 1
    assert n_cmp >= 6


def test_half_below_one_raises():
    ref, tgt, var, stamp, lab, nlab = _scene()
    with pytest.raises(ValueError):
        ad.kernel_fit_stamps(ref, tgt, var, 0, stamp, lab, nlab)


def test_stamp_never_in_own_fit():
    # a stamp whose target is wildly wrong must not move its own kernel
    ref, tgt, var, stamp, lab, nlab = _scene(seed=1)
    base = ad.kernel_fit_stamps(ref, tgt, var, 3, stamp, lab, nlab)
    bad = tgt.copy()
    bad[lab == 1] *= 5.0
    moved = ad.kernel_fit_stamps(ref, bad, var, 3, stamp, lab, nlab)
    assert np.allclose(moved[1], base[1], rtol=0, atol=1e-9 * np.abs(base[1]).max())
