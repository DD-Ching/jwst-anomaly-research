"""Offline tests for scripts/sunrise_scofield.py (synthetic maps; no network)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import sunrise_scofield as ss  # noqa: E402

STEP = 0.02  # arcsec per pixel, as the published maps
CENTRE = (24.3475, -8.4650)
THETA_E = 10.0  # arcsec


def _wcs() -> WCS:
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = list(CENTRE)
    w.wcs.crpix = [ss.NAXIS / 2, ss.NAXIS / 2]
    w.wcs.cdelt = [-STEP / 3600.0, STEP / 3600.0]
    w.wcs.set()
    return w


def _sis(wcs, name, i0, j0, half=ss.HALF):
    """SIS (alpha = theta_E along r) maps on the patch; x = West along +i, y = North along +j."""
    jj, ii = np.mgrid[j0 - half : j0 + half + 1, i0 - half : i0 + half + 1].astype(float)
    ra, dec = wcs.pixel_to_world_values(ii, jj)
    x = -(ra - CENTRE[0]) * np.cos(np.deg2rad(CENTRE[1])) * 3600.0
    y = (dec - CENTRE[1]) * 3600.0
    r = np.hypot(x, y)
    if name == "defx":
        return THETA_E * x / r
    if name == "defy":
        return THETA_E * y / r
    if name == "kappa":
        return THETA_E / (2 * r)
    return 1.0 / (1.0 - THETA_E / r)


@pytest.fixture
def sis(monkeypatch):
    wcs = _wcs()
    keys = {v[0]: k for k, v in ss.FILES.items()}
    monkeypatch.setattr(
        ss, "read_patch", lambda session, key, i0, j0, half=ss.HALF: _sis(wcs, keys[key], i0, j0)
    )
    return wcs


def test_measure_reproduces_sis_kappa_mu_and_source(sis):
    # 15" West, 5" North of the centre: kappa = theta_E / 2r, mu = r / (r - theta_E)
    ra = CENTRE[0] - 15.0 / 3600.0 / np.cos(np.deg2rad(CENTRE[1]))
    dec = CENTRE[1] + 5.0 / 3600.0
    row = ss.measure(None, sis, "x", ra, dec)
    r = np.hypot(15.0, 5.0)
    assert row["kappa_fd"] == pytest.approx(THETA_E / (2 * r), rel=2e-3)
    assert row["mu_fd"] == pytest.approx(r / (r - THETA_E), rel=5e-3)
    ax, ay = row["alpha_x"], row["alpha_y"]
    assert np.hypot(ax, ay) == pytest.approx(THETA_E, rel=1e-4)
    assert ax > 0 and ay > 0  # deflection points away from the centre (West, North)


def test_summarise_constraint_rms_and_held_out_offset():
    eye = [0.0, 0.0, 0.0]  # psi = 0: A = identity, so dtheta = dbeta
    rows = [
        {
            "image": "9.1a",
            "beta_x": 1.0,
            "beta_y": 0.0,
            "psi": eye,
            "mu_fd": 1.0,
            "constraint": True,
        },
        {
            "image": "9.1b",
            "beta_x": -1.0,
            "beta_y": 0.0,
            "psi": eye,
            "mu_fd": 1.0,
            "constraint": True,
        },
        {
            "image": "9.1c",
            "beta_x": 0.0,
            "beta_y": 3.0,
            "psi": eye,
            "mu_fd": 1.0,
            "constraint": False,
        },
    ]
    s = ss.summarise(rows)
    assert s["constraint_source_rms_arcsec"] == pytest.approx(1.0)
    held = s["9.1"]["images"]["9.1c"]
    assert held["dbeta_arcsec"] == pytest.approx(3.0)
    assert held["dtheta_linear_arcsec"] == pytest.approx(3.0)
