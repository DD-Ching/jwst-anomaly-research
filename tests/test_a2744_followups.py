from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

from jwst_anomaly import lensmodel

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
_spec = importlib.util.spec_from_file_location(
    "a2744_followups", ROOT / "scripts" / "a2744_followups.py"
)
af = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(af)


def _model():
    comps = [lensmodel.DPIE("halo", 0.0, 0.0, 0.0, 0.0, 1.0, 1000.0, 1000.0, 0.4)]
    return lensmodel.LensModel(comps, 64.0, -24.0, FlatLambdaCDM(H0=70, Om0=0.3))


def test_counter_images_round_lens():
    """A source just outside the Einstein ring of a round cored lens has a negative-parity
    counter-image on the opposite side, fainter than itself, and a demagnified central image; a
    source far outside the ring has no bright counter-image."""
    model = _model()
    grid = lensmodel.DeflectionGrid.compute(model, 60.0, 0.25)
    cat = Table({"label": [1, 2], "ra": [64.0, 64.01], "dec": [-24.0, -24.0], "mag": [20.0, 21.0]})
    cc = SkyCoord(cat["ra"], cat["dec"], unit="deg")
    # Einstein radius at z = 2, from the lens map along +x
    xs = np.linspace(1.0, 50.0, 4901)
    lm = model.lens_map(xs, np.zeros_like(xs), np.full_like(xs, 2.0))
    theta_e = float(xs[np.argmin(np.abs(lm["beta_x"]))])
    x = theta_e * 1.3
    ra, dec = model.to_sky(np.array([x]), np.array([0.0]))
    c = SkyCoord(float(ra[0]), float(dec[0]), unit="deg")
    mu, cps = af.counter_images(model, grid, x, 0.0, c, 2.0, 24.0, cat, cc)
    assert mu > 1.0
    bright = [p for p in cps if abs(p["mu"]) > 0.1]
    assert len(bright) == 1 and len(cps) <= 2
    assert bright[0]["mu"] < 0
    assert bright[0]["pred_mag_f200w"] > 24.0
    assert bright[0]["sep_from_source_arcsec"] > 1.5 * theta_e
    assert bright[0]["nearest_f200w"] == 2
    ra, dec = model.to_sky(np.array([55.0]), np.array([0.0]))
    c = SkyCoord(float(ra[0]), float(dec[0]), unit="deg")
    _, cps = af.counter_images(model, grid, 55.0, 0.0, c, 2.0, 24.0, cat, cc)
    assert all(abs(p["mu"]) < 0.1 for p in cps)


def test_knot_photometry_two_gaussians():
    """Two point-like knots of different colour: both found, background removed, fluxes ordered."""
    n = 101
    g = (np.arange(n) - n // 2) * 0.024
    xx, yy = np.meshgrid(g, g)

    def img(a1, a2):
        k1 = a1 * np.exp(-((xx - 0.3) ** 2 + yy**2) / (2 * 0.03**2))
        k2 = a2 * np.exp(-((xx + 0.3) ** 2 + yy**2) / (2 * 0.03**2))
        return 5.0 + k1 + k2

    err = np.full((n, n), 0.01)
    stamps = {"f200w": (img(1.0, 1.0), err, xx, yy), "f150w": (img(1.0, 0.1), err, xx, yy)}
    out = af.knot_photometry(stamps)
    assert len(out["knots_arcsec"]) == 2
    xs = [round(k[0], 1) for k in out["knots_arcsec"]]
    f150 = dict(zip(xs, out["flux"]["f150w"], strict=True))
    f200 = dict(zip(xs, out["flux"]["f200w"], strict=True))
    assert abs(f200[0.3][0] - f200[-0.3][0]) < 0.05 * f200[0.3][0]
    assert f150[-0.3][0] < 0.2 * f150[0.3][0]
