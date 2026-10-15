"""Unit tests for scripts/arc_diff.py on synthetic images (offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

_PATH = Path(__file__).resolve().parents[1] / "scripts" / "arc_diff.py"
_spec = importlib.util.spec_from_file_location("arc_diff", _PATH)
ad = importlib.util.module_from_spec(_spec)
sys.modules["arc_diff"] = ad
_spec.loader.exec_module(ad)


def _scene(rng: np.random.Generator, n: int = 128) -> np.ndarray:
    img = np.zeros((n, n))
    for _ in range(40):
        img += ad.gaussian_star((n, n), *rng.uniform(10, n - 10, 2), rng.uniform(1.0, 3.0), 50.0)
    return img


def test_register_recovers_subpixel_shift():
    rng = np.random.default_rng(1)
    a = _scene(rng)
    b = ad.shift(a, -1.3, 0.6)  # b is a moved by (-1.3, +0.6)
    dy, dx = ad.register(a, b)
    assert abs(dy - 1.3) < 0.1 and abs(dx + 0.6) < 0.1


def test_psf_match_blurs_the_sharper_epoch():
    rng = np.random.default_rng(2)
    e1 = _scene(rng)
    e2 = ad._blur(e1, 1.0)  # epoch 2 broader
    s, k, m1, m2 = ad.psf_match(e1, 1.03 * e2, np.round(np.linspace(-1.5, 1.5, 13), 3))
    assert 0.75 <= s <= 1.25  # blur epoch 1
    assert abs(k - 1.03) < 0.01
    assert ad._robust_std(m2 - m1) < 0.2 * ad._robust_std(e2 - e1)
    assert ad._robust_std(m2 - m1) < 0.2 * ad._robust_std(1.03 * e2 - e1)


def test_point_source_detected_and_noise_calibrated():
    rng = np.random.default_rng(3)
    n = 160
    noise = 1.0
    e1 = rng.normal(0, noise, (n, n))
    e2 = rng.normal(0, noise, (n, n))
    sig = 1.5
    flux = 12 * noise * np.sqrt(4 * np.pi) * sig  # ~12 sigma matched-filter
    e2 += ad.gaussian_star((n, n), 80.3, 70.6, sig, flux)
    var = np.full((n, n), 2 * noise**2)
    s, d = ad.significance(e2 - e1, var, sig, 21)
    assert 0.8 < ad._robust_std(s) < 1.2
    pk = ad.peaks(s, 5.0, 3)
    assert any(abs(y - 80) <= 1 and abs(x - 71) <= 1 and v > 8 for y, x, v in pk)
    assert len(pk) <= 3  # pure noise elsewhere rarely reaches 5 sigma


def test_group_links_bands_within_radius():
    d = [
        {"band": "F150W", "ra": 64.0, "dec": -24.0, "sigma": 6.0},
        {"band": "F277W", "ra": 64.0 + 0.05 / 3600, "dec": -24.0, "sigma": 7.0},
        {"band": "F356W", "ra": 64.0, "dec": -24.0 + 1.0 / 3600, "sigma": -6.0},
    ]
    g = sorted(ad.group(d, 0.1), key=len)
    assert [len(x) for x in g] == [1, 2]


def test_mag_to_flux_roundtrip():
    pixar = 2.3e-14
    f = ad.mag_to_mjysr(25.0, pixar)
    mag = -2.5 * np.log10(f * 1e6 * pixar / 3631.0)
    assert abs(mag - 25.0) < 1e-9


def test_local_shift_removes_a_varying_offset():
    rng = np.random.default_rng(4)
    a = _scene(rng, 192)
    # epoch 2: left half moved by +0.6 px in x, right half by -0.4 px (a smooth-ish warp)
    left = ad.shift(a, 0.0, 0.6)
    right = ad.shift(a, 0.0, -0.4)
    w = np.clip((np.arange(192) - 64) / 64.0, 0, 1)[None, :]
    b = (1 - w) * left + w * right
    glob = ad.shift(b, *ad.register(a, b))
    loc, loc2 = ad.local_shift(a, glob, 48, glob)
    assert np.allclose(loc, loc2, equal_nan=True)
    inner = (slice(20, -20), slice(20, -20))
    assert ad._robust_std((loc - a)[inner]) < 0.6 * ad._robust_std((glob - a)[inner])
