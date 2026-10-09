"""Tests for scripts/w12_niq_hst.py (synthetic stamps, offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.modeling import models

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("w12_niq_hst", _DIR / "w12_niq_hst.py")
hst = importlib.util.module_from_spec(_spec)
sys.modules["w12_niq_hst"] = hst
_spec.loader.exec_module(hst)

PIX = 0.04


def _pair(lens_flux=0.0, seed=0):
    rng = np.random.default_rng(seed)
    yy, xx = np.indices((120, 120))
    img = rng.normal(0, 0.05, (120, 120)) + 0.05
    for x, y, a in ((35, 60, 5.0), (87, 60, 3.0)):  # 2.08" apart
        img += models.Moffat2D(a, x, y, gamma=1.8, alpha=2.0)(xx, yy)
    if lens_flux:
        g = models.Sersic2D(1, 6.0, 1, x_0=66, y_0=60)(xx, yy)
        img += lens_flux * g / g.sum()
    return img


def _measure(img):
    p = hst.Params()
    peaks = hst.find_two_peaks(img, (35, 60), p.search_arcsec / PIX, 0.5 * 2.08 / PIX)
    fit, res = hst.two_psf_fit(img, peaks)
    clean = hst.remove_halos(res, fit)
    mad = 1.4826 * np.median(np.abs(res - np.median(res)))
    return hst.residual_between(clean, fit, PIX, p, mad)


def test_two_point_sources_leave_no_lens_light():
    m = _measure(_pair())
    assert m["sep_arcsec"] == pytest.approx(2.08, abs=0.05)
    assert abs(m["resid_snr"]) < 5


def test_injected_lens_galaxy_is_recovered():
    m = _measure(_pair(lens_flux=20.0))
    assert m["resid_snr"] > 5 and m["resid_flux"] == pytest.approx(20.0, rel=0.5)


def test_remove_halos_keeps_light_between_images():
    # a symmetric halo around each image is removed; a blob between them is kept
    yy, xx = np.indices((120, 120))
    fit = models.Const2D(0) + models.Moffat2D(1, 35, 60) + models.Moffat2D(1, 87, 60)
    halo = sum(np.exp(-np.hypot(xx - x, yy - 60) / 6) for x in (35, 87))
    blob = np.exp(-(np.hypot(xx - 61, yy - 60) ** 2) / 20)
    clean = hst.remove_halos(halo + blob, fit)
    assert clean[60, 61] == pytest.approx(blob[60, 61], abs=0.1)
    assert abs(clean[60, 20]) < 0.05
