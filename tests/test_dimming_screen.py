"""Tests for scripts/dimming_screen.py on synthetic multi-epoch catalogues (offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("dimming_screen", _DIR / "dimming_screen.py")
ds = importlib.util.module_from_spec(_spec)
sys.modules["dimming_screen"] = ds  # dataclasses resolve annotations through sys.modules
_spec.loader.exec_module(ds)

RA0, DEC0 = 150.0, 2.0
BANDS = ["F200W", "F444W"]


def _field(n=400, n_epochs=3, seed=3):
    """Sources on a 60" square: fluxes (Jy) with S/N 10-300, CI 2.3 (point-like)."""
    rng = np.random.default_rng(seed)
    cosd = np.cos(np.deg2rad(DEC0))
    ra = RA0 + rng.uniform(-30, 30, n) / 3600 / cosd
    dec = DEC0 + rng.uniform(-30, 30, n) / 3600
    flux = 10 ** rng.uniform(-8.0, -6.5, n)  # 10 nJy .. 0.3 µJy
    return rng, ra, dec, flux


def _cat(rng, ra, dec, flux, err, keep, zp=1.0):
    k = np.flatnonzero(keep)
    noisy = flux[k] + rng.normal(0, err, k.size)
    t = Table(
        {
            "label": k + 1,
            "ra": ra[k],
            "dec": dec[k],
            "aper50_flux": noisy * zp,
            "aper50_flux_err": np.full(k.size, err * zp),
            "CI_70_30": np.full(k.size, 2.3),
            "nn_dist": np.full(k.size, 50.0),
        }
    )
    t.meta["pixel_scale_arcsec"] = 0.031
    return t


def _cats(n_epochs=3, vanish=None, dim=None, chroma=None, zp=(1.0, 1.1, 0.95)):
    rng, ra, dec, flux = _field(n_epochs=n_epochs)
    err = 1e-9  # 1 nJy
    cats = {}
    for k in range(n_epochs):
        for jb, b in enumerate(BANDS):
            f = flux.copy()
            keep = np.ones(len(f), bool)
            if vanish is not None and k == vanish[1]:
                keep[vanish[0]] = False
            if dim is not None and k == dim[1]:
                f[dim[0]] *= 1 - dim[2]
            if chroma is not None and k == chroma[1] and jb == 0:
                f[chroma[0]] *= 1 - chroma[2]
            cats[(k, b)] = _cat(rng, ra, dec, f, err, keep, zp[k])
    return cats, flux


def _screen(cats, n_epochs=3):
    p = ds.Params()
    lc = ds.build_light_curves(cats, n_epochs, BANDS, p)
    scale = ds.noise_scale(lc["flux"], lc["err_raw"], lc["det"])
    err = ds.scaled_errors(lc["flux"], lc["err_raw"], scale, p.sys_floor)
    res = ds.classify(ds.mask_shallow_nondetections(lc["flux"], err, lc["det"], p), err, p)
    return lc, err, res, scale


def _index(lc, ra, dec, i):
    d = np.hypot((lc["ra"] - ra[i]) * np.cos(np.deg2rad(DEC0)), lc["dec"] - dec[i]) * 3600
    return int(np.argmin(d))


def test_vanish_dim_and_chromatic_dim():
    rng, ra, dec, flux = _field()
    v = int(np.argmax(flux))  # brightest source vanishes in epoch 2
    d = int(np.argsort(flux)[-2])  # second brightest dims by 50 % in both bands in epoch 1
    c = int(np.argsort(flux)[-3])  # third dims in F200W only
    cats, _ = _cats(vanish=(v, 2), dim=(d, 1, 0.5), chroma=(c, 1, 0.5))
    lc, err, res, scale = _screen(cats)
    assert len(lc["ra"]) == len(flux)
    # zero points tied: epoch factors undo the injected (1, 1.1, 0.95)
    assert abs(lc["zero_point_factor"]["F200W_e1"] * 1.1 - 1) < 0.02
    iv, idim, ic = (_index(lc, ra, dec, i) for i in (v, d, c))
    assert res["vanish"][iv] and res["vanish_bands"][iv] == 2
    assert res["dim_achromatic"][idim] and res["dim_epoch"][idim] == 1
    assert abs(res["max_drop"][idim] - 0.5) < 0.05
    assert not res["dim_achromatic"][ic] and res["dim_bands"][ic] == 1
    others = np.ones(len(lc["ra"]), bool)
    others[[iv, idim, ic]] = False
    assert ds.flagged(res)[others].sum() <= 2  # steady sources: no flags beyond noise
    assert np.all(scale >= 1.0)


def test_shallow_nondetection_is_not_a_vanish():
    # a 12 nJy source with 1 nJy errors, missing from an epoch whose error is 3 nJy: S/N there
    # would have been 4, so the missing row says nothing
    flux = np.array([[[12.0], [0.0]]])
    err = np.array([[[1.0], [3.0]]])
    det = np.array([[[True], [False]]])
    p = ds.Params()
    res = ds.classify(ds.mask_shallow_nondetections(flux, err, det, p), err, p)
    assert not res["vanish"][0]
    err[0, 1, 0] = 1.0  # deep enough: 12 sigma expected, so the non-detection counts
    res = ds.classify(ds.mask_shallow_nondetections(flux, err, det, p), err, p)
    assert res["vanish"][0]


def test_rise_dip_rise_needs_three_epochs():
    f = np.array([[[100.0], [40.0], [100.0]], [[100.0], [100.0], [40.0]]])
    e = np.full_like(f, 2.0)
    res = ds.classify(f, e)
    assert res["rise_dip_rise"].tolist() == [True, False]
    assert res["dim_achromatic"].tolist() == [False, False]  # one band only


def test_w3_factor_umbra_capped_spike_and_noise():
    times = np.array([0.0, 1.0, 2.0])
    # row 0: lens passes behind epoch 1 with u0 = 0 and t_E = 0.3 yr: the umbra (beta < 2) lasts
    # |t - t0| < 0.6 yr, so only epoch 1 vanishes; at t = 0 and 2, beta = 3.3 and A = 1.03
    # row 1: t0 far away -> no effect
    fac = ds.w3_factor(times, np.array([1.0, 50.0]), 0.3, np.array([0.0, 0.5]), 0.01)
    assert fac[0, 1] == 0.0
    assert np.allclose(fac[0, [0, 2]], 1.03, atol=0.01)
    assert np.allclose(fac[1], 1.0, atol=1e-3)
    # near the caustic crossing (beta = 2 + rho/2) the spike reaches the D-047 peak;
    # a cap is honoured
    t_spike = 2.0 + 0.005
    spike = ds.w3_factor(np.array([t_spike]), np.array([0.0]), 1.0, np.array([0.0]), 0.01)
    assert 0.8 * ds.SPIKE_PEAK[0.01] < spike[0, 0] < 1.1 * ds.SPIKE_PEAK[0.01]
    capped = ds.w3_factor(
        np.array([t_spike]), np.array([0.0]), 1.0, np.array([0.0]), 0.01, spike_peak=2.0
    )
    assert capped[0, 0] == 2.0
    half = ds.w3_factor(times, np.array([1.0]), 0.3, np.array([0.0]), 0.01, blend=0.5)
    assert np.isclose(half[0, 1], 0.5)  # the unlensed half of the flux stays
    # applying the factor: a vanished epoch keeps only sky noise (the floor part is removed)
    flux = np.full((1, 3, 2), 100.0)
    err = np.hypot(1.0, 0.03 * flux)
    out, e2 = ds.apply_factor(flux, err, fac[:1], np.random.default_rng(0))
    assert np.all(np.abs(out[0, 1]) < 5.0) and np.allclose(e2[0, 1], 1.0, atol=0.2)
    assert np.allclose(out[0, 0], 103.0, atol=1.0)
    dimf = ds.dimming_factor(3, np.array([2]), 0.2)
    assert dimf.tolist() == [[1.0, 1.0, 0.8]]


def test_efficiency_and_limits():
    cats, _ = _cats()
    lc, err, res, _ = _screen(cats)
    n = len(lc["ra"])
    sel = np.ones(n, bool)
    veto = np.zeros(n, bool)
    times = np.array([0.0, 0.5, 1.0])
    eff = ds.efficiency_table(
        lc, err, veto, sel, times, (0.3,), (0.1,), (0.2, 1.0), seed=1, n_rep=2
    )
    assert eff.meta["provenance"] == "simulated"
    full = eff[(eff["model"] == "dimming") & (eff["depth"] == 1.0) & (eff["n_injected"] > 20)]
    assert len(full) and np.all(full["efficiency"] > 0.6)
    w3 = eff[eff["model"] == "w3"]
    assert np.isclose(w3["window_yr"][0], 1.0 + 4 * 0.3)
    n_mon = int(np.sum(w3["n_sources"]))
    lim = ds.rate_limits(eff, area_deg2=0.01, n_monitored=n_mon)
    assert lim.meta["source"]
    r = lim[0]
    assert np.isclose(r["limit_per_source_per_yr"], 3.0 / r["exposure_source_yr"])
    assert np.isclose(r["limit_tau"], r["limit_per_source_per_yr"] * np.pi * 0.3)
    assert np.isclose(
        r["limit_per_deg2_per_yr"], 3.0 / (0.01 * r["mean_efficiency"] * r["window_yr"])
    )
    assert np.isclose(r["limit_per_deg2_per_epoch"], r["limit_per_deg2_per_yr"] * np.pi * 0.3)


def test_ordinary_columns_star_self_match_and_veto():
    cats, flux = _cats()
    lc, err, res, _ = _screen(cats)
    i = 0
    gaia = Table(
        {
            "source_id": [1, 2],
            "ra": [lc["ra"][i], RA0 + 0.5],
            "dec": [lc["dec"][i], DEC0 + 0.5],
            "gmag": [20.5, 12.0],
        }
    )
    o = ds.ordinary_columns(lc, cats, BANDS, ds.Params(), gaia)
    assert o["gaia_star"][i] and not o["near_star"][i]  # a faint star is itself, not a mask
    gaia["gmag"][0] = 15.0
    o = ds.ordinary_columns(lc, cats, BANDS, ds.Params(), gaia)
    assert o["near_star"][i] and ds.catalogue_veto(o)[i]  # a saturating star is vetoed
    assert o["point_like"].mean() > 0.9


def test_area_proxy():
    ra = np.array([0.0, 0.0001, 0.01])
    dec = np.zeros(3)
    assert np.isclose(ds.area_deg2(ra, dec, np.array([2, 2, 1])), (5 / 3600) ** 2)


def test_combine_limits_adds_exposures_and_drops_uncalibrated():
    def lim(exp, area_exp, calibrated):
        t = Table(
            {
                "t_e_yr": [0.3],
                "rho": [0.1],
                "n_monitored": [10],
                "exposure_source_yr": [exp],
                "area_exposure_deg2_yr": [area_exp],
            }
        )
        t.meta.update(calibrated=calibrated)
        return t

    per = {"a": lim(10.0, 0.01, True), "b": lim(20.0, 0.02, True), "c": lim(99.0, 9.0, False)}
    r = ds.combine_limits(per)[0]
    assert np.isclose(r["limit_per_source_per_yr"], 3.0 / 30.0)
    assert np.isclose(r["limit_per_deg2_per_yr"], 3.0 / 0.03)
    assert np.isclose(r["limit_tau"], 0.1 * np.pi * 0.3)
    assert r["n_monitored"] == 20
    assert ds.combine_limits(per, calibrated_only=False)[0]["n_monitored"] == 30


def test_bright_neighbour_and_gaia_self_match():
    from astropy.coordinates import SkyCoord

    # a faint source with its own faint Gaia match, 2" from a G = 15 star: masked by the star
    pos = SkyCoord([150.0, 150.0], [2.0, 2.0 + 2.0 / 3600], unit="deg")
    gaia = Table({"ra": [150.0, 150.0], "dec": [2.0, 2.0 + 2.0 / 3600], "gmag": [20.5, 15.0]})
    near, star = ds.gaia_proximity(pos, gaia)
    assert star.all() and near[0]
    gaia["gmag"][1] = 21.0  # two faint Gaia stars (1.5" radius): only the self matches at 2"
    near, _ = ds.gaia_proximity(pos, gaia)
    assert not near.any()


def test_saturated_star_masks_neighbours_beyond_bright_neighbour_radius():
    from astropy.coordinates import SkyCoord

    # a G = 17.9 star (NaN catalogue flux, so never a bright_neighbour) and a faint source 1.8"
    # away: near_star must cover it, and check_params guards that coupling
    p = ds.Params()
    sep = p.bright_neighbour_arcsec + p.gaia_self_arcsec
    pos = SkyCoord([150.0], [2.0 + sep / 3600], unit="deg")
    gaia = Table({"ra": [150.0], "dec": [2.0], "gmag": [p.gaia_saturated_g - 0.1]})
    near, _ = ds.gaia_proximity(pos, gaia, p)
    assert near[0]
    ds.check_params(p)
    with pytest.raises(ValueError):
        ds.check_params(ds.Params(gaia_saturated_g=21.0))


def test_known_object_labels_and_failed_services():
    xm = Table(
        {
            "source_uid": ["b", "a"],
            "is_known_object": [True, False],
            "best_match_service": ["ned", ""],
            "best_match_type": ["G", ""],
            "best_match_id": ["WISEA J001412.34-302345.6 extra long identifier", ""],
            "n_simbad": [0, 0],
            "n_ned": [1, -1],
        }
    )
    labels, failed = ds.known_object_labels(xm, ["a", "b"])
    assert labels == ["", "ned:G:WISEA J001412.34-302345.6 extra long identifier"]
    assert failed == ["ned"]
    xm["n_ned"] = [1, 0]
    assert ds.known_object_labels(xm, ["a", "b"])[1] == []
