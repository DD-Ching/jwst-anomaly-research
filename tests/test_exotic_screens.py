"""Tests for scripts/exotic_screens.py (synthetic inputs only)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.cosmology import FlatLambdaCDM

from jwst_anomaly.lensmodel import DPIE, LensModel

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("exotic_screens", _DIR / "exotic_screens.py")
es = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(es)

RA0, DEC0 = 110.826989, -73.454723


def test_luminosity_ratio_flags_a_demagnified_sibling():
    flux = np.array([10.0, 5.0, 4.0, 0.4, 3.0, 2.0])
    err = np.full(6, 0.05)
    mu = np.array([10.0, 5.0, 4.0, 2.0, 3.0, 80.0])  # last: |mu| > 50 is unusable
    lr = es.luminosity_ratios(flux, err, mu, ["1", "1", "1", "1", "2", "2"])
    np.testing.assert_allclose(lr["lum_ratio"][:3], [1.0, 1.0, 1.0])
    assert abs(lr["lum_ratio"][3] - 0.2) < 1e-9 and lr["ratio_snr"][3] > 5
    assert np.isnan(lr["lum_ratio"][4]) and not lr["usable"][5]  # system 2 has one usable image


def test_flux_class():
    assert es.flux_class(0.2, 0.25, 8.0, 6.0) == "underluminous"
    assert es.flux_class(4.0, 3.5, 8.0, 6.0) == "overluminous"
    assert es.flux_class(0.2, 0.9, 8.0, 6.0) == "chromatic"
    assert es.flux_class(0.2, 0.25, 8.0, 3.0) == "consistent"  # not significant in band 2
    assert es.flux_class(np.nan, 1.0, 1.0, 1.0) == "untestable"


def test_radial_lines_converge_on_their_centre():
    ang = np.deg2rad([10.0, 100.0, 190.0, 280.0])
    x, y = 10 + 5 * np.cos(ang), 10 + 5 * np.sin(ang)
    # position angle E of N of the direction towards the centre (x West, y North)
    pa = np.mod(np.rad2deg(np.arctan2(-(10 - x), 10 - y)), 180.0)
    g = np.arange(0, 20.01, 0.5)
    gx, gy = np.meshgrid(g, g)
    counts = es.line_counts(x, y, pa, gx, gy, tol=0.5, max_len=10.0)
    peaks = es.convergence_peaks(counts, 3)
    assert peaks and all(counts[p] == 4 for p in peaks)
    assert all(np.hypot(gx[p] - 10, gy[p] - 10) <= 0.75 for p in peaks)
    tangential = np.mod(pa + 90.0, 180.0)  # tangential arcs do not converge
    assert not es.convergence_peaks(es.line_counts(x, y, tangential, gx, gy, 0.5, 10.0), 3)


def test_radial_magnification_of_an_isothermal_lens_is_one():
    model = LensModel(
        [DPIE("halo", 0.0, 0.0, 0.0, 0.0, 1e-3, 1e4, 900.0, 0.39)],
        RA0,
        DEC0,
        FlatLambdaCDM(H0=70.0, Om0=0.3),
    )
    ra = RA0 + np.array([10.0, 20.0]) / 3600 / np.cos(np.deg2rad(DEC0))
    mu_r = es.radial_magnification(model, ra, np.full(2, DEC0), 2.0)
    np.testing.assert_allclose(mu_r, 1.0, atol=3e-3)  # SIS: kappa = gamma, so 1 - kappa + gamma = 1


def test_spike_segments_are_dropped_but_tangential_neighbours_kept():
    from astropy.table import Table

    cosd = np.cos(np.deg2rad(DEC0))
    star = (RA0, DEC0)
    # two sources 5" N and 5" E of the star; the first along the spike (PA 0), the second across
    ras = [RA0, RA0 + 5 / 3600 / cosd]
    decs = [DEC0 + 5 / 3600, DEC0]
    src = Table({"ra": ras, "dec": decs, "pa_obs": [2.0, 0.0]})
    shapes = Table(
        {
            "ra": [star[0], *ras],
            "dec": [star[1], *decs],
            "mag": [16.0, 24.0, 24.0],
            "is_extended": [False, True, True],
        }
    )
    assert list(es.spike_segments(src, shapes)) == [True, False]
    shapes["mag"][0] = 23.0  # a faint point source has no long spikes (radius 3")
    assert not es.spike_segments(src, shapes).any()
    np.testing.assert_allclose(es.spike_radius([20.0, 17.0, 10.0]), [3.0, 3 * 10**0.6, 20.0])


def test_spike_veto_keeps_off_axis_radial_arcs():
    from astropy.table import Table

    cosd = np.cos(np.deg2rad(DEC0))
    star = (RA0, DEC0)
    # 6 spike segments on a hexagonal set at theta = 10 deg, plus one radial arc at PA 35 deg
    pas = [10.0, 70.0, 130.0, 10.0, 70.0, 130.0, 35.0]
    seps = [4.0, 4.0, 4.0, 7.0, 7.0, 7.0, 5.0]
    ras, decs = [], []
    for pa, sep in zip(pas, seps, strict=True):
        ras.append(RA0 + sep * np.sin(np.deg2rad(pa)) / 3600 / cosd)
        decs.append(DEC0 + sep * np.cos(np.deg2rad(pa)) / 3600)
    src = Table({"ra": ras, "dec": decs, "pa_obs": pas})
    shapes = Table(
        {
            "ra": [star[0], *ras],
            "dec": [star[1], *decs],
            "mag": [15.0] + [24.0] * len(ras),
            "is_extended": [False] + [True] * len(ras),
        }
    )
    veto = es.spike_segments(src, shapes)
    assert list(veto) == [True] * 6 + [False]


def test_external_stars_veto_spikes_of_stars_missing_from_the_catalogue(tmp_path):
    from astropy.table import Table

    cosd = np.cos(np.deg2rad(DEC0))
    # a saturated G = 13 star with no catalogue row; segments 30" and 33" N along its spike,
    # one 30" E across it
    ras = [RA0, RA0, RA0 + 30 / 3600 / cosd]
    decs = [DEC0 + 30 / 3600, DEC0 + 33 / 3600, DEC0]
    src = Table({"ra": ras, "dec": decs, "pa_obs": [1.0, 179.0, 0.0]})
    shapes = Table({"ra": ras, "dec": decs, "mag": [24.0] * 3, "is_extended": [True] * 3})
    assert not es.spike_segments(src, shapes).any()  # no catalogued star, nothing vetoed
    path = tmp_path / "gaia.ecsv"
    Table({"RA_ICRS": [RA0], "DE_ICRS": [DEC0], "Gmag": [13.0]}).write(path, format="ascii.ecsv")
    stars = es.read_spike_stars(path)
    assert list(stars.colnames) == ["ra", "dec", "mag"] and stars.meta["provenance"] == "observed"
    assert list(es.spike_segments(src, shapes, stars)) == [True, True, False]
    # a catalogued star of the same brightness keeps the 20" cap (D-034): the 30" segments survive;
    # passing the Gaia star as well seeds it once (deduplicated), from the catalogue
    shapes2 = Table(
        {
            "ra": [RA0, *ras],
            "dec": [DEC0, *decs],
            "mag": [13.0] + [24.0] * 3,
            "is_extended": [False] + [True] * 3,
        }
    )
    assert not es.spike_segments(src, shapes2).any()
    assert not es.spike_segments(src, shapes2, stars).any()
    np.testing.assert_allclose(es.spike_radius(12.0, cap=es.EXTERNAL_SPIKE_CAP), 60.0)
    # masked magnitudes are ignored; a missing column is a clear error
    m = Table({"ra": [RA0, RA0], "dec": [DEC0, DEC0], "Gmag": [13.0, 0.0]}, masked=True)
    m["Gmag"].mask = [False, True]
    m.write(tmp_path / "masked.ecsv", format="ascii.ecsv")
    assert np.isnan(es.read_spike_stars(tmp_path / "masked.ecsv")["mag"][1])
    Table({"ra": [RA0], "dec": [DEC0], "G": [13.0]}).write(
        tmp_path / "bad.ecsv", format="ascii.ecsv"
    )
    with pytest.raises(ValueError, match="no mag column"):
        es.read_spike_stars(tmp_path / "bad.ecsv")


def test_gaia_star_table_conversion():
    from astropy.table import MaskedColumn, Table

    spec = importlib.util.spec_from_file_location("gaia_stars", _DIR / "gaia_stars.py")
    gs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gs)
    v = Table({"RA_ICRS": [1.0, 2.0], "DE_ICRS": [3.0, 4.0]})
    v["Gmag"] = MaskedColumn([12.5, 0.0], mask=[False, True])
    out = gs.to_star_table(v)
    assert (
        out.colnames == ["ra", "dec", "mag"] and out["mag"][0] == 12.5 and np.isnan(out["mag"][1])
    )
    assert len(gs.to_star_table(None)) == 0


@pytest.mark.network
def test_gaia_stars_live_query(tmp_path):
    spec = importlib.util.spec_from_file_location("gaia_stars", _DIR / "gaia_stars.py")
    gs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gs)
    out = tmp_path / "stars.ecsv"
    gs.main(
        ["39.97134", "-1.58226", "--radius-arcmin", "4", "--out", str(out)]
    )  # Abell 370 (D-043)
    stars = es.read_spike_stars(out)
    assert len(stars) >= 5 and np.nanmin(stars["mag"]) < 14


# ------------------------------------------------------------------------------- shear (D-050)


def test_aperture_filters_vanish_outside_the_annulus():
    x = np.array([0.05, 0.2, 0.5, 1.0, 1.2])
    for kind in ("pointmass", "schirmer", "tophat"):
        q = es.aperture_filter(x, kind, 0.1)
        assert q[0] == 0 and q[-1] == 0 and np.all(q[1:4] >= 0)
    np.testing.assert_allclose(es.aperture_filter(x, "pointmass", 0.1)[1:4], 1 / x[1:4] ** 2)
    with pytest.raises(ValueError):
        es.aperture_filter(x, "gauss", 0.1)


def test_intrinsic_ellipticity_phase_and_psf():
    # q = 0.5 ellipses (sigma 4 x 2 px), no PSF: |eps| = (1 - q)/(1 + q) = 1/3, phase 2 PA
    eps = es.intrinsic_ellipticity([4.0, 4.0], [2.0, 2.0], [0.0, 45.0], 0.0)
    np.testing.assert_allclose(eps, [1 / 3, 1j / 3], atol=1e-12)
    # a PSF-sized round source has no measurable shape; an unresolved one is NaN
    assert abs(es.intrinsic_ellipticity([2.0], [2.0], [30.0], 1.0)[0]) < 1e-12
    assert np.isnan(es.intrinsic_ellipticity([1.0], [0.8], [0.0], 1.0)[0])


def test_cluster_shear_removal_inverts_the_lens_mapping():
    rng = np.random.default_rng(0)
    eps_s = 0.3 * rng.uniform(0, 1, 50) * np.exp(1j * rng.uniform(0, 2 * np.pi, 50))
    g = 0.3 * np.exp(1j * 0.7)
    eps_obs = (eps_s + g) / (1 + np.conj(g) * eps_s)  # Seitz & Schneider 1997, |g| < 1
    np.testing.assert_allclose(es.remove_cluster_shear(eps_obs, g), eps_s, atol=1e-12)


def _ring(n, r, radial: bool):
    """Sources on a ring about the origin (model frame, x West), radial or tangential."""
    phi = np.linspace(0, 2 * np.pi, n, endpoint=False)  # PA east of north
    x, y = -r * np.sin(phi), r * np.cos(phi)
    pa = phi if radial else phi + np.pi / 2
    return x, y, 0.2 * np.exp(2j * pa)


def test_aperture_mass_sign_radial_is_positive_tangential_negative():
    x, y, e_rad = _ring(12, 3.0, radial=True)
    _, _, e_tan = _ring(12, 3.0, radial=False)
    ap = es.ApertureMass(x, y, [0.0, 50.0], [0.0, 0.0], 6.0, "pointmass", 1.0, min_n=5)
    s_rad, sx_rad = ap.snr(e_rad)
    s_tan, _ = ap.snr(e_tan)
    # perfectly aligned: S = sqrt(2 n) = sqrt(24)
    np.testing.assert_allclose([s_rad[0], s_tan[0]], [np.sqrt(24), -np.sqrt(24)])
    assert abs(sx_rad[0]) < 1e-9
    assert np.isnan(s_rad[1])  # no sources in the second aperture


def test_aperture_mass_null_rarely_reaches_a_strong_radial_ring():
    rng = np.random.default_rng(1)
    bx, by = rng.uniform(-30, 30, (2, 600))
    be = 0.25 * np.exp(1j * rng.uniform(0, 2 * np.pi, 600))
    rx, ry, re = _ring(16, 2.5, radial=True)
    x, y, e = np.r_[bx, rx + 10], np.r_[by, ry], np.r_[be, re]
    gx, gy = es.radial_grid(30.0, 1.0)
    ap = es.ApertureMass(x, y, gx, gy, 6.0, "pointmass", 1.0)
    s, _ = ap.snr(e)
    s = s.reshape(gx.shape)
    iy, ix = np.unravel_index(np.nanargmax(s), s.shape)
    assert np.hypot(gx[iy, ix] - 10, gy[iy, ix]) <= 1.5
    null = ap.null_max(e, 100, np.random.default_rng(2))
    assert np.mean(null >= s[iy, ix]) < 0.05


def test_shear_responsivity_recovers_a_diluted_shear():
    rng = np.random.default_rng(3)
    n = 4000
    g = rng.uniform(0.05, 0.4, n) * np.exp(1j * rng.uniform(0, 2 * np.pi, n))
    eps = 0.25 * np.exp(1j * rng.uniform(0, 2 * np.pi, n)) + 0.5 * g
    r, err = es.shear_responsivity(eps, g)
    assert abs(r - 0.5) < 3 * err and err < 0.05
    with pytest.raises(ValueError):  # too few rows: refuse rather than assume R = 1
        es.shear_responsivity(eps[:5], g[:5])
