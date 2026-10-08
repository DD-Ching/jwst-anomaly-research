"""Tests for scripts/inject_shear.py (synthetic inputs only)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
sys.path.insert(0, str(Path(__file__).resolve().parent))
_spec = importlib.util.spec_from_file_location("inject_shear", _DIR / "inject_shear.py")
ish = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ish)

from test_inject_radial import SCALE, _model, _shapes  # noqa: E402


def _field(n=900, seed=0):
    rng = np.random.default_rng(seed)
    model = _model()
    x, y = rng.uniform(-40, 40, (2, n))
    shapes = _shapes(
        model,
        x,
        y,
        pa_obs=rng.uniform(0, 180, n),
        ellipticity=rng.uniform(0.1, 0.5, n),
        semimajor_px=rng.uniform(3.0, 6.0, n),
    )
    args = SimpleNamespace(
        **(
            ish.SHEAR_DEFAULTS
            | {"max_radius": 40.0, "recover_tol": 2.0, "n_random": 40, "responsivity": 1.0}
        )
    )
    return model, shapes, args


def test_shear_defaults_come_from_the_screen():
    assert ish.SHEAR_DEFAULTS == ish.es.shear_defaults()
    assert (
        ish.SHEAR_DEFAULTS["filter"] == "schirmer" and ish.SHEAR_DEFAULTS["aperture_arcsec"] == 10
    )


def test_injected_ellipticity_keeps_r_of_the_measured_change_only():
    e_src = np.array([0.1 + 0.05j, np.nan])
    raw_src = np.array([0.3 + 0.0j, np.nan])  # measured, cluster shear included
    raw_img = np.array([0.5 + 0.1j, 0.2j])
    out = ish.injected_ellipticity(e_src, raw_img, raw_src, 0.5)
    # resolved source: the measured change only; an unresolved source is never painted (NaN)
    np.testing.assert_allclose(out[0], 0.1 + 0.05j + 0.5 * (0.2 + 0.1j))
    assert np.isnan(out[1])
    # no lens change: the source's corrected shape comes back unchanged, whatever its cluster g
    same = ish.injected_ellipticity(e_src[:1], raw_src[:1], raw_src[:1], 0.45)
    np.testing.assert_allclose(same, e_src[:1])
    # |e| is capped below 1 (a large inner-image change on an already elliptical source)
    big = ish.injected_ellipticity([0.8 + 0j], [0.95 + 0j], [0.1 + 0j], 0.45)
    np.testing.assert_allclose(big, [0.99 + 0j])


def test_only_measurable_rows_are_lensed():
    model, shapes, args = _field()
    shapes["semimajor_px"][:50] = 0.5  # unresolved after deconvolving a 1 px PSF
    inj = ish.ShearInjector(model, shapes, args, psf_sigma=1.0)
    assert not inj.background[:50].any()
    assert inj.background[50:].sum() > 0.9 * (len(shapes) - 50)


def test_a_user_responsivity_outside_its_range_is_refused():
    model, shapes, _ = _field()
    for bad in (0.0, -0.3, 2.0):
        with pytest.raises(ValueError, match="outside"):
            ish.es.shear_sources(model, shapes, 1.0, 10.0, responsivity=bad)


def test_a_massive_w1_lens_is_recovered_and_no_lens_is_not():
    model, shapes, args = _field()
    inj = ish.ShearInjector(model, shapes, args, psf_sigma=1.0)
    assert inj.r == 1.0  # synthetic shapes carry no cluster shear to measure R from
    inj.draw_null(np.random.default_rng(1))
    big = np.where(inj.background, 6.0, np.nan)  # theta_E 6": dozens of radial images
    hit = inj.trial(15.0, -10.0, big, SCALE)
    assert hit["n_lensed"] > 5 and hit["n_image_used"] > 5
    assert hit["recovered"] and hit["s_best"] > np.max(inj.null)
    none = np.full(len(shapes), np.nan)  # no lensable rows: nothing painted
    miss = inj.trial(15.0, -10.0, none, SCALE)
    assert miss["n_lensed"] == 0 and miss["s_best"] < hit["s_best"]


def test_painted_images_on_a_spike_axis_are_vetoed_like_real_segments():
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy.table import Table

    model, shapes, args = _field()
    star = SkyCoord(shapes["ra"][0], shapes["dec"][0], unit="deg")
    shapes["mag"][0], shapes["is_extended"][0] = 16.0, False  # spikes reach ~19"
    for i, sep in enumerate(np.arange(2.0, 10.0), start=1):  # real segments along PA 0 (north)
        c = star.directional_offset_by(0.0 * u.deg, sep * u.arcsec)
        shapes["ra"][i], shapes["dec"][i], shapes["pa_obs"][i] = c.ra.deg, c.dec.deg, 0.0
        shapes["ellipticity"][i] = 0.6
    inj = ish.ShearInjector(model, shapes, args, psf_sigma=1.0)
    assert inj.base["spike"][1:9].all() and len(inj.spike_ref) >= 8
    # painted: on the north axis and aligned (vetoed); 30 deg off every axis, pointing at the star
    pas = [0.0, 30.0, 0.0]
    pos = [star.directional_offset_by(pa * u.deg, 6.5 * u.arcsec) for pa in (180.0, 30.0)]
    pos.append(star.directional_offset_by(0.0 * u.deg, 40.0 * u.arcsec))  # beyond the spike
    img = Table(
        {
            "ra": [c.ra.deg for c in pos],
            "dec": [c.dec.deg for c in pos],
            "pa_obs": pas,
            "label": [10**6, 10**6 + 1, 10**6 + 2],
        }
    )
    np.testing.assert_array_equal(inj.spike_veto(img), [True, False, False])
