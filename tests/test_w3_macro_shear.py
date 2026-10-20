"""Point lens inside a macro model (exotic_sim.macro_magnification_map); coarse grids, fast."""

import math

import numpy as np
import pytest

from jwst_anomaly import exotic_sim as es


def _img(t):
    ny, nx = t.meta["shape"]
    return np.asarray(t["mu"]).reshape(ny, nx)


@pytest.mark.parametrize("sign", [1, -1])
def test_isolated_matches_closed_form(sign):
    t = es.macro_magnification_map(1.0, 1.0, sign, 3.0, 3.0, pixel=0.05, ray_step=0.01)
    assert t.meta["provenance"] == "simulated"
    r = np.hypot(t["beta_x"], t["beta_y"])
    exact = es.total_magnification(np.asarray(r), 1.0, sign)
    sel = (exact > 0.5) & (exact < 5) & (np.abs(r - 2) > 0.2)
    assert np.allclose(np.asarray(t["mu"])[sel], exact[sel], rtol=0.05)
    rad = 3.0
    got = (np.asarray(t["mu"])[r < rad] - 1).sum() * 0.05**2
    pred = math.pi * (rad * math.sqrt(rad * rad + 4 * sign) - rad * rad)
    assert got == pytest.approx(pred, rel=0.05)


def test_minimum_image_keeps_a_narrowed_umbra():
    lam1 = 0.1
    t = es.macro_magnification_map(lam1, 1.0, -1, 1.5, 3.0, pixel=0.02, ray_step=0.01)
    mu = _img(t)
    ny, nx = mu.shape
    row = mu[ny // 2] * lam1 < 0.1
    half = row.sum() * 0.02 / 2
    assert half == pytest.approx(es.macro_hole_halfwidths(lam1, 1.0)[0], abs=0.04)
    assert mu[ny // 2, nx // 2] == 0


def test_saddle_image_has_no_umbra():
    t = es.macro_magnification_map(-0.1, 1.0, -1, 1.5, 3.0, pixel=0.02, ray_step=0.01)
    assert (np.asarray(t["mu"]) * 0.1).min() > 0.5
    assert es.macro_hole_halfwidths(-0.1, 1.0) == (0.0, 2.0)


def test_rejects_macro_critical_curve():
    with pytest.raises(ValueError):
        es.macro_magnification_map(0.0, 1.0)
