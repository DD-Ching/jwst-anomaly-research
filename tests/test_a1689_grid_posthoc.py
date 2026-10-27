"""Tests for scripts/a1689_grid_posthoc.py."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import a1689_grid_posthoc as ap  # noqa: E402
from transient_forced import compare  # noqa: E402


def test_expected_sigma_matches_compare_noise_free():
    f1 = np.array([1.0, 0.2, 5.0])
    e1, e2 = np.array([0.05, 0.03, 0.1]), np.array([0.02, 0.02, 0.05])
    inj = f1.copy()
    # compare with no zero-point refs (min_zp_refs above n) gives the same error model
    _, _, sig = compare(f1, e1, f1 + inj, e2, 0.05, len(f1) + 1)
    z = ap.expected_sigma(f1, e1, f1, e2, inj, 1.0)
    assert z == pytest.approx(sig)
    assert ap.expected_sigma(f1, e1, f1, e2, inj, 2.0) == pytest.approx(sig / 2)


def test_recovery_by_bins_counts_and_open_top():
    z = np.array([0.5, 3.5, 4.5, 5.5, 7.0, 9.0, 9.5])
    rec = np.array([False, False, True, True, True, True, False])
    b = ap.recovery_by_bins(z, rec)
    assert [x["n"] for x in b] == [1, 1, 1, 1, 1, 2]
    assert b[-1]["z_hi"] is None and b[-1]["recovered_frac"] == 0.5
    assert b[0]["recovered_frac"] == 0.0


def test_recovery_by_bins_empty_bin_is_none():
    b = ap.recovery_by_bins(np.array([7.0]), np.array([True]))
    assert b[0] == {"z_lo": 0.0, "z_hi": 3.0, "n": 0, "recovered_frac": None}


def test_expected_sigma_matches_band_pass_with_zero_point_and_scale():
    import transient_grid as tg

    f1 = np.array([1.0, 0.3, 4.0, 0.05])
    e1, e2 = np.array([0.04, 0.03, 0.1, 0.02]), np.array([0.02, 0.01, 0.05, 0.01])
    inj, zp, scale = np.array([1.0, 0.3, 4.0, 0.2]), 0.12, 1.4
    k = 10 ** (-0.4 * zp)
    f2_base = f1 * k  # no change apart from the injection, on the epoch-2 zero point
    m = {"f1": f1, "e1": e1, "f2": f2_base + inj, "e2": e2}
    sig_c = tg.band_pass(m, 0.05, 0.3, 5.0, scale, zp)[1]
    assert ap.expected_sigma(f1 * k, e1 * k, f2_base, e2, inj, scale) == pytest.approx(sig_c)
