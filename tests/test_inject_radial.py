"""Tests for scripts/inject_radial.py (synthetic catalogue and lens model only; offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from astropy.cosmology import FlatLambdaCDM, Planck18
from astropy.table import Table

from jwst_anomaly import exotic_sim
from jwst_anomaly.lensmodel import DPIE, LensModel

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("inject_radial", _DIR / "inject_radial.py")
ir = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ir)

RA0, DEC0 = 110.826989, -73.454723


def _model() -> LensModel:
    return LensModel(
        [DPIE("halo", 0.0, 0.0, 0.0, 0.0, 1e-3, 1e4, 400.0, 0.39)],
        RA0,
        DEC0,
        FlatLambdaCDM(H0=70.0, Om0=0.3),
    )


def _shapes(model, x, y, **cols) -> Table:
    ra, dec = model.to_sky(np.asarray(x, float), np.asarray(y, float))
    n = len(ra)
    base = {
        "label": np.arange(1, n + 1),
        "ra": ra,
        "dec": dec,
        "pa_obs": np.zeros(n),
        "ellipticity": np.zeros(n),
        "semimajor_px": np.full(n, 3.0),
        "area_px": np.full(n, 50.0),
        "snr": np.full(n, 200.0),
        "is_extended": np.ones(n, bool),
        "mag": np.full(n, 24.0),
    }
    base.update({k: np.asarray(v) for k, v in cols.items()})
    return Table(base)


def _args(**kw) -> SimpleNamespace:
    return SimpleNamespace(**(ir.SCREEN_DEFAULTS | {"max_radius": 60.0, "recover_tol": 2.0} | kw))


def test_mass_matches_the_d047_scaling():
    # docs/exotic_lensing.md: theta_E = 2.2" (|M| / 1e12 Msun)^(1/2), z_l = 0.4, z_s = 2, Planck18
    m = ir.theta_e_to_mass(2.2, 0.4, 2.0, Planck18)
    assert 0.95e12 < m < 1.05e12
    np.testing.assert_allclose(ir.theta_e_to_mass(4.4, 0.4) / m, 4.0, rtol=1e-9)


def test_surface_density_limit_is_three_over_effective_area():
    lim = ir.surface_density_limit([0.5, 0.25], [1e-3, 2e-3])
    np.testing.assert_allclose(lim, ir.POISSON_UL_95 / 1e-3)
    assert ir.surface_density_limit([0.0], [1.0]) == float("inf")


def test_round_source_becomes_radial_with_the_jacobian_axis_ratio():
    # image 2" east of the lens; lam_r = 0.5 and lam_t = 1.5; round source, PSF sigma 1 px
    semi, ell, pa = ir.lensed_shapes(
        np.array([3.0]),
        np.array([3.0]),
        np.array([0.0]),
        np.array([2.0]),
        np.array([0.0]),
        np.array([0.5]),
        np.array([1.5]),
        1.0,
    )
    rad, tan = 8.0 / 0.25 + 1.0, 8.0 / 2.25 + 1.0  # intrinsic sigma^2 = 9 - 1
    np.testing.assert_allclose(semi, np.sqrt(rad))
    np.testing.assert_allclose(ell, 1.0 - np.sqrt(tan / rad))
    np.testing.assert_allclose(pa, 90.0)  # along the east-west line through the lens
    img = exotic_sim.inject_images([2.5], [0.0], 1.0, n=1.0, sign=-1)
    assert np.all(img["radial"]) and np.allclose(img["pa_deg"], 90.0)


def test_paint_lens_removes_the_umbra_and_replaces_lensed_sources():
    model = _model()
    xl, yl, te = 30.0, 0.0, 2.0
    # model frame x is West: a source at x = xl - d lies d arcsec East of the lens
    x = [xl - 1.0, xl - 4.4, xl - 12.0, xl - 5.0, xl + 5.0]
    y = [0.0, 0.0, 0.0, 0.0, 0.0]
    shapes = _shapes(model, x, y, is_extended=[1, 1, 1, 0, 1], mag=[24, 24, 24, 15, 24])
    shapes["z_phot"] = [np.nan, np.nan, np.nan, np.nan, 0.1]  # the last is foreground
    shapes["z160"] = shapes["z_phot"]
    xs, ys = model.to_frame(shapes["ra"], shapes["dec"])
    bg = ir.background_mask(shapes, model.z_lens)
    assert list(bg) == [True, True, True, False, False]
    keep, img, info = ir.paint_lens(shapes, xs, ys, bg, model, xl, yl, te, 1.0, 3.0, 100)
    # umbra source (beta = 0.5) gone; beta = 2.2 source replaced by two images; the rest stay
    assert list(keep) == [False, False, True, True, True]
    assert info == {"n_umbra": 1, "n_lensed": 1, "n_images": 2}
    assert list(img["label"]) == [100, 101] and img.meta["provenance"] == "simulated"
    ix, iy = model.to_frame(img["ra"], img["dec"])
    xo = np.sort(xl - np.asarray(ix))  # east offsets of the images
    np.testing.assert_allclose(xo, te * np.array([2.2 - np.sqrt(0.84), 2.2 + np.sqrt(0.84)]) / 2)
    np.testing.assert_allclose(iy, 0.0, atol=1e-6)
    np.testing.assert_allclose(img["pa_obs"], 90.0)  # radial: along the line to the lens
    sol = exotic_sim.solve_images([2.2], 1.0, -1)
    mu = np.abs(sol["mu"][0])
    np.testing.assert_allclose(np.sort(img["mag"]), np.sort(24.0 - 2.5 * np.log10(mu)))
    np.testing.assert_allclose(np.sort(img["snr"]), np.sort(200.0 * np.sqrt(mu)))


def test_injector_recovers_a_lens_with_many_radial_images_and_not_an_empty_one():
    model = _model()
    rng = np.random.default_rng(0)
    # sparse random background arcs, plus round sources at beta = 2.1 theta_E east and west of
    # a lens site at (30, 0), where east-west is "anti" to the cluster's tangential stretch
    bx, by = rng.uniform(-55, 55, 60), rng.uniform(-55, 55, 60)
    far = np.hypot(bx - 30.0, by) > 15.0  # nothing else behind the lens site
    bx, by = bx[far][:40], by[far][:40]
    te = 3.0
    ang = np.deg2rad([0.0, 8.0, -8.0, 180.0, 172.0, 188.0])
    sx = 30.0 - 2.1 * te * np.cos(ang)
    sy = 2.1 * te * np.sin(ang)
    x, y = np.concatenate([bx, sx]), np.concatenate([by, sy])
    n = len(x)
    ell = np.concatenate([np.full(40, 0.6), np.zeros(6)])
    pa = np.concatenate([rng.uniform(0, 180, 40), np.zeros(6)])
    shapes = _shapes(model, x, y, ellipticity=ell, pa_obs=pa)
    inj = ir.RadialInjector(model, shapes, _args())
    assert inj.base_counts_info["n_elongated"] == int(np.sum(np.hypot(bx, by) <= 60.0))
    np.testing.assert_array_equal(inj.max_rand, inj.rand.reshape(len(inj.rand), -1).max(axis=1))
    hit = inj.trial(30.0, 0.0, te, np.random.default_rng(1), 1.0, 3.0)
    assert hit["n_lensed"] == 6 and hit["n_image_arcs"] >= 8
    assert hit["recovered"] and hit["n_lines"] >= 8 and hit["sep"] <= 2.0
    # a lens where nothing lies behind it changes nothing and is not recovered
    miss = inj.trial(-50.0, -50.0, 0.3, np.random.default_rng(1), 1.0, 3.0)
    assert miss["n_lensed"] == 0 and not miss["recovered"]
    assert n == len(shapes)  # the input catalogue is untouched


def test_footprint_covers_only_where_the_catalogue_has_sources():
    rng = np.random.default_rng(2)
    xs, ys = rng.uniform(0, 40, 4000), rng.uniform(-40, 40, 4000)  # half of a 40" disc
    _, _, area = ir.screened_footprint(xs, ys, 40.0)
    half_disc = 0.5 * np.pi * 40.0**2
    assert 0.95 * half_disc < area < 1.15 * half_disc  # edge grows by up to FOOTPRINT_RADIUS
