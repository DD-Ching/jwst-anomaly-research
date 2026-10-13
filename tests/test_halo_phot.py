"""Synthetic checks of scripts/halo_phot.py (no network)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import halo_phot as hp  # noqa: E402

PIX = 0.031  # arcsec, NIRCam SW i2d


def _grid(half=3.0):
    n = int(np.ceil(half / PIX))
    j, i = np.mgrid[-n : n + 1, -n : n + 1].astype(float)
    return i * PIX, j * PIX


def _halo(xx, yy, x0=1.0, y0=-0.5, amp=5.0, rs=1.5):
    return amp * np.exp(-np.hypot(xx - x0, yy - y0) / rs)


@pytest.mark.parametrize(
    ("route", "amp", "rel"),
    [
        ("medfilt", 0.5, 0.05),  # halo light in the aperture ~ the source's: 3 % bias
        ("annulus", 0.5, 0.15),  # a plane misses the cusp's curvature: +10 %
        ("medfilt", 5.0, 0.40),  # halo ~10x the source: 0.3 mag residual bias, still << raw (12x)
    ],
)
def test_compact_source_on_halo_recovered(route, amp, rel):
    xx, yy = _grid()
    flux = 30.0
    sci = _halo(xx, yy, amp=amp) + hp.gaussian(xx, yy, 0.0, 0.0, flux, 0.05)
    err = np.full(sci.shape, 0.01)
    f, e, conc = hp.measure(sci, err, xx, yy, route)
    assert f == pytest.approx(flux, rel=rel)
    assert f >= flux  # residual halo light only ever brightens the image
    assert amp > 1 or conc > hp.MIN_CONCENTRATION  # a strong residual halo also lowers conc
    assert e > 0


def test_halo_alone_adds_flux_without_subtraction():
    xx, yy = _grid()
    sci = _halo(xx, yy, amp=0.5) + hp.gaussian(xx, yy, 0.0, 0.0, 30.0, 0.05)
    raw = float(np.sum(sci[np.hypot(xx, yy) <= hp.R_AP_ARCSEC]))
    assert raw > 1.5 * 30.0  # the halo is a real bias here, so the test above is not trivial


def test_resolved_source_fails_concentration():
    xx, yy = _grid()
    sci = 2.0 * np.exp(-np.hypot(xx, yy) / 0.5)  # exponential disc, scale 0.5"
    _, _, conc = hp.measure(sci, np.full(sci.shape, 0.01), xx, yy, "annulus")
    assert conc < hp.MIN_CONCENTRATION


@pytest.mark.parametrize("fwhm", [0.05, 0.145])
def test_ring_checks(fwhm):
    xx, yy = _grid()
    sci = _halo(xx, yy, amp=0.5) + hp.gaussian(xx, yy, 0.0, 0.0, 30.0, fwhm)
    err = np.full(sci.shape, 0.01)
    checks = hp.ring_checks(sci, err, xx, yy, {"medfilt": 30.0, "annulus": 30.0}, fwhm)
    for blank, rec in checks.values():
        assert blank.shape == rec.shape == (8,)
        assert abs(np.median(rec) - 1) < 0.05  # run() uses the median; the cusp side can recentre
        assert np.median(np.abs(blank)) < 0.1 * 30.0  # a mild halo is mostly removed on the ring


def test_blank_apertures_see_a_strong_halo_residual():
    xx, yy = _grid()
    sci = _halo(xx, yy, amp=5.0, x0=0.3, y0=0.0)  # cusp next to the ring
    blank, _ = hp.ring_checks(sci, np.full(sci.shape, 0.01), xx, yy, {"annulus": 30.0}, 0.05)[
        "annulus"
    ]
    assert np.max(np.abs(blank)) > 0.5  # the check has power: it sees uncorrected halo light


def test_recentre_ignores_halo_gradient():
    xx, yy = _grid()
    sci = _halo(xx, yy, amp=1.0, x0=0.6, y0=0.0) + hp.gaussian(xx, yy, -0.05, 0.0, 30.0, 0.05)
    err = np.full(sci.shape, 0.01)
    # both routes centre on the compact image (-0.05, 0), not up the halo gradient (+x)
    for route in ("medfilt", "annulus"):
        f = hp.measure(sci, err, xx, yy, route)[0]
        f_true = hp.measure(sci, err, xx, yy, route, -0.05, 0.0, recentre=0.0)[0]
        assert f == pytest.approx(f_true, rel=0.02)


def test_nan_in_aperture_gives_nan():
    xx, yy = _grid()
    sci = hp.gaussian(xx, yy, 0.0, 0.0, 30.0, 0.05)
    sci[np.hypot(xx, yy) < 0.05] = np.nan
    f, e, c = hp.measure(sci, np.full(sci.shape, 0.01), xx, yy, "medfilt", recentre=0.0)
    assert np.isnan(f) and np.isnan(e) and np.isnan(c)


def test_unknown_route():
    xx, yy = _grid(0.5)
    with pytest.raises(ValueError):
        hp.measure(np.zeros(xx.shape), np.ones(xx.shape), xx, yy, "kron")


@pytest.fixture(autouse=True)
def _band_b(monkeypatch):
    """The synthetic band "b" gets an axis-ratio threshold, as a star-calibrated band would."""
    monkeypatch.setitem(hp.MIN_AXIS_RATIO, "b", 0.8)


def _rows(src_mags, mu=2.0, conc=0.7, err_frac=0.01, systems=None, q=0.95):
    """Synthetic measured rows: one band, both routes, given source magnitudes."""
    systems = systems or ["1"] * len(src_mags)
    rows = []
    for i, (m, s) in enumerate(zip(src_mags, systems, strict=True)):
        f = 10 ** (-0.4 * (m - 2.5 * np.log10(mu)))
        row = {"image_id": f"{s}.{i}", "system": s, "mu": mu, "b_conc": conc, "b_axis_ratio": q}
        for route in hp.ROUTES:
            row[f"b_{route}_flux"], row[f"b_{route}_err"] = f, err_frac * f
        rows.append(row)
    return rows


def test_classify_flags_only_the_worst_image():
    rows = _rows([20.0, 20.05, 21.0])  # third image 1 mag fainter than its siblings
    hp.classify(rows, ["b"])
    assert [r["b_flux_class"] for r in rows] == ["consistent", "consistent", "flux_outlier"]
    assert rows[2]["b_medfilt_resid"] == pytest.approx(0.975, abs=1e-6)


def test_classify_needs_both_routes_and_gates():
    rows = _rows([20.0, 20.05, 21.0])
    rows[2]["b_annulus_flux"] = rows[0]["b_annulus_flux"]  # second route disagrees
    hp.classify(rows, ["b"])
    assert rows[2]["b_flux_class"] == "consistent"
    rows = _rows([20.0, 20.05, 21.0])
    rows[2]["b_annulus_flux"] = float("nan")  # one route only
    hp.classify(rows, ["b"])
    assert rows[2]["b_flux_class"] == "one_route"
    rows = _rows([20.0, 20.05, 21.0], q=0.6)  # elongated (a knot on an arc): no test
    hp.classify(rows, ["b"])
    assert {r["b_flux_class"] for r in rows} == {"untested"}
    rows = _rows([20.0, 21.0], conc=0.5)  # resolved: no test
    hp.classify(rows, ["b"])
    assert {r["b_flux_class"] for r in rows} == {"untested"}
    rows = _rows([20.0, 21.0], err_frac=0.2)  # S/N 5 < 10
    hp.classify(rows, ["b"])
    assert {r["b_flux_class"] for r in rows} == {"untested"}


def test_measure_image_recovers_an_injected_dimming():
    """End to end on a synthetic stamp: corrections keep a 1 mag dimming at ~1 mag."""
    xx, yy = _grid()
    err = np.full(xx.shape, 0.01)
    outs = []
    for flux in (30.0, 30.0 * 10**-0.4):
        sci = _halo(xx, yy, amp=0.5) + hp.gaussian(xx, yy, 0.0, 0.0, flux, 0.145)
        outs.append(hp.measure_image(sci, err, xx, yy, 2.0, "f444w"))
    for route in hp.ROUTES:
        a, b = (o[f"f444w_{route}_flux"] for o in outs)
        assert -2.5 * np.log10(b / a) == pytest.approx(1.0, abs=0.1)
        assert outs[0][f"f444w_{route}_err"] > 0


@pytest.mark.parametrize("dm", [1.0, -1.0])
def test_scaled_image_changes_flux_by_dm(dm):
    """--inject-mag scales the image's own filtered light; the halo stays."""
    xx, yy = _grid()
    err = np.full(xx.shape, 0.01)
    sci = _halo(xx, yy, amp=0.5) + hp.gaussian(xx, yy, 0.0, 0.0, 30.0, 0.145)
    f0 = hp.measure(sci, err, xx, yy, "medfilt")[0]
    f1 = hp.measure(*hp._scaled_image((sci, err, xx, yy), dm)[:2], xx, yy, "medfilt")[0]
    assert -2.5 * np.log10(f1 / f0) == pytest.approx(dm, abs=0.05)


@pytest.mark.parametrize("fwhm", [0.05, 0.145])
@pytest.mark.parametrize("offset", [(0.0, 0.0), (0.5, 0.0), (0.25, 0.5)])
def test_axis_ratio_point_source_is_round(fwhm, offset):
    """Also off a pixel centre: moments about the peak pixel made such a source elongated."""
    xx, yy = _grid()
    x0, y0 = (o * PIX for o in offset)
    sci = _halo(xx, yy, amp=0.5) + hp.gaussian(xx, yy, x0, y0, 30.0, fwhm)
    f = hp.halo_filtered(sci, PIX)
    assert hp.axis_ratio(f, xx, yy) > 0.9


def test_axis_ratio_thin_arc_with_knot_is_elongated():
    """El Gordo 1c's failure mode: a knot on a thin tangential arc passes the concentration gate."""
    xx, yy = _grid()
    arc = 0.3 * np.exp(-0.5 * (yy / 0.03) ** 2) * (np.abs(xx) <= 2.5)  # 5" long, unresolved width
    sci = arc + hp.gaussian(xx, yy, 0.0, 0.0, 4.0, 0.05)
    f = hp.halo_filtered(sci, PIX)
    _, _, conc = hp.measure(sci, np.full(sci.shape, 0.01), xx, yy)
    assert conc > hp.MIN_CONCENTRATION  # the old gate lets it through
    assert hp.axis_ratio(f, xx, yy) < hp.MIN_AXIS_RATIO["f150w"]


def test_run_refuses_band_without_star_calibration(monkeypatch):
    monkeypatch.setattr(hp.lc, "has_image_list", lambda m: True)
    monkeypatch.setattr(hp.lc, "is_map_model", lambda m: False)
    with pytest.raises(SystemExit, match="axis-ratio"):
        hp.run("any", ["1"], {"f200w": "unused"})
