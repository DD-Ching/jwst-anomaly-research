"""Offline tests for scripts/m1149_refs.py (synthetic SIS map; no network)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import m1149_refs as mr  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

CENTRE = (177.3988, 22.3985)
THETA_E = 5.0
PIX = 0.05


def _sis_model(half: float = 20.0) -> lensmodel.MapLensModel:
    n = int(round(2 * half / PIX)) + 1
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = list(CENTRE)
    w.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    w.wcs.cdelt = [-PIX / 3600.0, PIX / 3600.0]
    g = (np.arange(n) - (n - 1) / 2) * PIX
    x, y = np.meshgrid(g, g)
    r = np.hypot(x, y)
    r[r == 0] = 1e-6
    return lensmodel.MapLensModel(
        THETA_E * x / r, THETA_E * y / r, w, 0.542, FlatLambdaCDM(H0=70, Om0=0.3), centre=CENTRE
    )


def test_second_model_images_recovers_the_sis_pair():
    model = _sis_model()
    grid = lensmodel.DeflectionGrid.compute(model, 15.0, 0.1)
    z = 2.0
    te = THETA_E * model.dls_ds(z)
    beta = np.array([1.0, 0.5])
    u = beta / np.hypot(*beta)
    a, b = beta + te * u, beta - te * u
    ra, dec = model.to_sky(np.array([a[0], b[0]]), np.array([a[1], b[1]]))
    cat = Table({"ra": ra, "dec": dec})
    imgs = mr.second_model_images(model, grid, cat, z)
    imgs = imgs[np.abs(imgs["magnification"]) >= 0.5]  # drop the map core's central image
    assert len(imgs) == 2
    for x, y in (a, b):
        assert np.min(np.hypot(imgs["x"] - x, imgs["y"] - y)) < 0.05


def test_d041_verdict():
    assert mr.d041_verdict(np.nan, False) == "model-dependent (CATS only)"
    assert mr.d041_verdict(mr.XMODEL_LOOSE_ARCSEC + 0.1, True) == "model-dependent (CATS only)"
    assert mr.d041_verdict(1.5, True) == "catalogued image misplaced by CATS"
    assert mr.d041_verdict(0.2, False) == "both models"


def test_sharon_model_rejects_unpinned_maps(tmp_path):
    p = tmp_path / "x.fits"
    p.write_bytes(b"not the pinned map")
    with pytest.raises(ValueError, match="sha256"):
        mr.sharon_model(p, p)
