"""Point lens inside a macro model (exotic_sim.macro_magnification_map); coarse grids, fast."""

import math
import sys
from pathlib import Path

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
    for lams in [(-0.1, 1.0), (0.0, 1.0), (-0.1, -1.0)]:
        with pytest.raises(ValueError):
            es.macro_hole_halfwidths(*lams)


def test_rejects_macro_critical_curve():
    with pytest.raises(ValueError):
        es.macro_magnification_map(0.0, 1.0)


def test_no_lens_control_is_the_macro_magnification():
    t = es.macro_magnification_map(0.5, 1.0, 0, 1.0, 1.0, pixel=0.05, ray_step=0.01)
    assert np.allclose(np.asarray(t["mu"]), 2.0, rtol=0.05)


def test_non_integer_bin_count_keeps_the_macro_normalisation():
    # 2 * 0.97 / 0.05 = 38.8 bins rounds to 39: mu must use the real bin width, not `pixel`
    t = es.macro_magnification_map(0.5, 1.0, 0, 0.97, 1.0, pixel=0.05, ray_step=0.01)
    assert t.meta["pixel_x"] == pytest.approx(2 * 0.97 / 39)
    assert np.mean(np.asarray(t["mu"])) == pytest.approx(2.0, rel=0.002)


def test_gate_pass_rules():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import w3_macro_shear as w

    assert w.g1_pass([0.001, 0.029]) and not w.g1_pass([0.001, 0.03])
    assert w.g2_pass(0.9, 1.0) and not w.g2_pass(0.89, 1.0) and not w.g2_pass(1.0, 0.89)
    assert w.g3_pass([0.019, -0.019]) and not w.g3_pass([0.0, -0.02])
