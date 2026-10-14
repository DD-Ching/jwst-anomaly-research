"""Member-galaxy subtraction (scripts/member_sub.py) on synthetic stamps (simulated)."""

import sys
from pathlib import Path

import numpy as np
from astropy.modeling import models

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import member_sub as ms  # noqa: E402

PIX = 0.063  # arcsec, NIRCam long-wavelength mosaics
HALF = 4.1


def _stamp(member_xy=(1.2, 0.0), member_amp=20.0, image_flux=0.0, image_xy=(0.0, 0.0), seed=1):
    n = int(round(HALF / PIX))
    jj, ii = np.mgrid[-n : n + 1, -n : n + 1].astype(float)
    xx, yy = ii * PIX, jj * PIX
    rng = np.random.default_rng(seed)
    err = np.full(xx.shape, 0.01)
    sci = rng.normal(0.0, 0.01, xx.shape)
    gal = models.Sersic2D(
        amplitude=member_amp,
        r_eff=0.25,
        n=2.0,
        x_0=member_xy[0],
        y_0=member_xy[1],
        ellip=0.2,
        theta=0.3,
    )
    sci += gal(xx, yy)
    if image_flux:
        sci += ms._gaussian(xx - image_xy[0], yy - image_xy[1], image_flux)
    return sci, err, xx, yy


NB = [(1.2, 0.0, 0.2, 0.3)]


def test_recovers_image_next_to_member():
    sci, err, xx, yy = _stamp(image_flux=2.0)
    res = ms.measure(sci, err, xx, yy, pred_flux=2.0, search=1.0, nb=NB)
    assert res["n_components"] == 1
    assert res["recovery"] > 0.9
    assert res["member_class"] == "recovered", res


def test_offset_image_within_search_is_found():
    sci, err, xx, yy = _stamp(image_flux=2.0, image_xy=(-0.5, 0.4))
    res = ms.measure(sci, err, xx, yy, pred_flux=2.0, search=1.0, nb=NB)
    assert res["member_class"] == "recovered", res
    assert abs(res["best_dx"] + 0.5) < 0.15 and abs(res["best_dy"] - 0.4) < 0.15


def test_absent_when_no_image():
    sci, err, xx, yy = _stamp()
    res = ms.measure(sci, err, xx, yy, pred_flux=2.0, search=1.0, nb=NB)
    assert res["member_class"] == "absent", res


def test_member_core_on_position_is_a_blend():
    sci, err, xx, yy = _stamp(member_xy=(0.2, 0.0))
    res = ms.measure(sci, err, xx, yy, pred_flux=2.0, search=1.0, nb=[(0.2, 0.0, 0.2, 0.3)])
    assert res["member_class"] == "blend"


def test_faint_neighbour_is_not_modelled():
    sci, err, xx, yy = _stamp(member_amp=0.001)
    res = ms.measure(sci, err, xx, yy, pred_flux=2.0, search=1.0, nb=NB)
    assert res["member_class"] == "no_member"


def test_class_thresholds():
    assert ms.member_class(0.2, 10.0) == "absent"
    assert ms.member_class(1.0, 10.0) == "recovered"
    assert ms.member_class(5.0, 10.0) == "confused"
    assert ms.member_class(1.0, 2.0) == "ambiguous"
    assert ms.member_class(float("nan"), 2.0) == "ambiguous"
