"""Frozen stage 4 rules of scripts/m1931_muse.py (offline)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import m1931_muse as mm  # noqa: E402


def test_line_window_and_continuum():
    lam = 6900.0 + 1.25 * np.arange(128)
    lnom = 6980.0
    cm = mm.continuum_mask(lam, lnom)
    assert not cm[np.abs(lam - lnom) <= 30].any() and cm[np.abs(lam - lnom) > 30].all()
    spec = np.zeros_like(lam)
    k = int(np.argmin(np.abs(lam - (lnom + 10))))
    spec[k] = 5.0
    spec[0] = 50.0  # outside the peak search: ignored
    lo, hi, pk = mm.line_window(spec, lam, lnom)
    assert pk == k and (lo, hi) == (k - 4, k + 5)


def test_narrow_band_subtracts_continuum():
    cube = np.full((40, 20, 20), 2.0)
    cube[18:21, 10, 10] += 1.0
    lam = 7000 + 1.25 * np.arange(40) - 25
    sub = mm.subtract_continuum(cube, np.abs(lam - 7000) > 12)
    nb = mm.narrow_band(sub, 16, 25, 1.25)
    assert np.isclose(nb[10, 10], 3.75) and np.isclose(nb[0, 0], 0.0)


def test_classify_order_and_thresholds():
    assert mm.classify(False, False, 100, 10, 1) == "off_cube"
    assert mm.classify(True, False, 100, np.nan, 1) == "off_cube"
    assert mm.classify(True, True, 100, 0, 1) == "blended"
    assert mm.classify(True, False, 4.9, 0, 1) == "undetectable_local"
    assert mm.classify(True, False, 30, 10, 1) == "recovered"
    assert mm.classify(True, False, 30, 91, 1) == "confused"
    assert mm.classify(True, False, 30, 9.9, 1) == "absent"
    assert mm.classify(True, False, 9, 2, 1) == "ambiguous"


def test_gate_l():
    assert mm.gate_l(["recovered", "confused", "recovered", "absent", "blended"])["passes"]
    assert not mm.gate_l(["recovered", "confused", "recovered", "off_cube"])["passes"]  # 3 testable
    assert not mm.gate_l(["recovered", "absent", "recovered", "ambiguous"])["passes"]  # 50 %


def test_blending_and_windows():
    assert mm.is_blended([2.0, 1.5]) and not mm.is_blended([1.7]) and not mm.is_blended([])
    assert mm.windows_overlap((1, 3), (2, 4)) and not mm.windows_overlap((1, 2), (2, 3))


def test_max_in_radius_finds_offset_source():
    nb = np.zeros((60, 60))
    nb[30, 33] = 10.0  # 0.6" from the centre
    assert mm.max_in_radius(nb, 30, 30, 1.0) == 10.0
    assert mm.max_in_radius(nb, 30, 30, 0.0) == 10.0  # within the 0.8" aperture of the centre
    assert np.isnan(mm.max_in_radius(nb, 1, 1, 0.0))  # aperture leaves the image


def test_local_noise_reproducible():
    rng = np.random.default_rng(0)
    nb = rng.normal(size=(100, 100))
    a = mm.local_noise(nb, 50, 50, [50.0], [50.0], np.random.default_rng(mm.SEED))
    b = mm.local_noise(nb, 50, 50, [50.0], [50.0], np.random.default_rng(mm.SEED))
    assert a == b and 4 < a < 9  # ~sqrt(49 spaxels)
