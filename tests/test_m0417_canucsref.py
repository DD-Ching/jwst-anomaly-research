"""Column mapping of scripts/m0417_canucsref.py (pixel scale, sky PA, depth, catalogue form)."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

_SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_SCRIPTS))
_spec = importlib.util.spec_from_file_location("m0417_canucsref", _SCRIPTS / "m0417_canucsref.py")
mc = importlib.util.module_from_spec(_spec)
sys.modules["m0417_canucsref"] = mc
_spec.loader.exec_module(mc)


def _grid(scale=0.04):
    x, y = np.meshgrid(np.arange(0, 1000, 100.0), np.arange(0, 1000, 100.0))
    x, y = x.ravel(), y.ravel()
    dec0 = -24.0
    ra = 64.0 - (x - 500) * scale / 3600.0 / np.cos(np.radians(dec0))  # east left
    dec = dec0 + (y - 500) * scale / 3600.0
    return x, y, ra, dec


def test_pixel_scale_north_up():
    assert abs(mc.pixel_scale(*_grid(0.04)) - 0.04) < 1e-6


def test_sky_pa_convention():
    # major axis along +y (North) is PA_image 90 -> sky 0; along +x (West) -> sky 90
    assert mc.sky_pa([90.0, 0.0, -45.0]).tolist() == [0.0, 90.0, 45.0]


def test_depth_mag():
    f = np.array([3631e-9 * 10 ** (-0.4 * 3), 1.0, 2.0])  # first row: AB 25.5
    s = np.array([5.0, 100.0, 3.0])
    assert abs(mc.depth_mag(f, s) - (-2.5 * np.log10(f[0] / 3631.0))) < 1e-9
    assert np.isnan(mc.depth_mag(f[1:], s[1:]))


def test_to_cat_mapping():
    x, y, ra, dec = _grid()
    n = len(x)
    flux = np.full(n, 100.0)
    flux[0] = np.nan
    t = Table(
        {
            "SOURCE": np.arange(n),
            "X": x,
            "Y": y,
            "RA": ra,
            "DEC": dec,
            "A": np.full(n, 4.0),
            "B": np.full(n, 1.0),
            "PA": np.full(n, 90.0),
            "FLUX_KRON_TOTAL_F200W": flux,
            "FLUXERR_KRON_TOTAL_F200W": np.full(n, 10.0),
            "FLAG_POINTSRC": np.zeros(n, bool),
        }
    )
    c = mc.to_cat(t)
    assert len(c["flux"]) == n - 1  # non-finite flux dropped
    assert np.allclose(c["flux"], 100e-9) and np.allclose(c["snr"], 10.0)
    assert np.allclose(c["a"], 0.16) and np.allclose(c["b"], 0.04)
    assert np.allclose(c["pa"], 0.0) and np.allclose(c["ell"], 0.75)
    assert c["ext"].all()


def test_run_requires_all_overrides():
    import m0417_arcref as ar
    import pytest

    with pytest.raises(ValueError):
        ar.run(Path("."), cat={}, depth=None, cat_inputs=None)


def test_posthoc_aperture_edges_and_nan():
    spec = importlib.util.spec_from_file_location(
        "m0417_canucsref_posthoc", _SCRIPTS / "m0417_canucsref_posthoc.py"
    )
    ph = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ph)
    d = np.ones((50, 50))
    assert ph.aperture_sum(d, 25, 25, 2.0) == 13.0  # pixels within r = 2
    assert np.isnan(ph.aperture_sum(d, 1, 25, 2.0))  # leaves the image
    d[25, 25] = np.nan
    assert np.isnan(ph.aperture_sum(d, 25, 25, 2.0))


def test_to_cat_pa_conventions():
    x, y, ra, dec = _grid()
    n = len(x)
    t = Table(
        {
            "SOURCE": np.arange(n),
            "X": x,
            "Y": y,
            "RA": ra,
            "DEC": dec,
            "A": np.full(n, 4.0),
            "B": np.full(n, 1.0),
            "PA": np.full(n, 30.0),
            "FLUX_KRON_TOTAL_F200W": np.full(n, 100.0),
            "FLUXERR_KRON_TOTAL_F200W": np.full(n, 10.0),
            "FLAG_POINTSRC": np.zeros(n, bool),
        }
    )
    assert np.allclose(mc.to_cat(t)["pa"], 120.0)
    assert np.allclose(mc.to_cat(t, 0.0, -1.0)["pa"], 150.0)
