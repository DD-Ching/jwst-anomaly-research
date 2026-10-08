"""Tests for scripts/exotic_screens.py (synthetic inputs only)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
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
    # a catalogued star of the same brightness is capped at 20", so the 30" segments survive
    shapes2 = Table(
        {
            "ra": [RA0, *ras],
            "dec": [DEC0, *decs],
            "mag": [13.0] + [24.0] * 3,
            "is_extended": [False] + [True] * 3,
        }
    )
    assert not es.spike_segments(src, shapes2).any()
    np.testing.assert_allclose(es.spike_radius(12.0, cap=es.EXTERNAL_SPIKE_CAP), 60.0)
