import importlib.util
from pathlib import Path

import numpy as np
from astropy.table import Table

_P = Path(__file__).resolve().parents[1] / "scripts" / "w1_psf_floor.py"
_spec = importlib.util.spec_from_file_location("w1_psf_floor", _P)
w1p = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(w1p)
es = w1p.es


def test_rules():
    assert w1p.r1_map_feasible(np.full(50, 9.0), 10)
    assert not w1p.r1_map_feasible(np.full(50, 9.0), 9)
    assert not w1p.r1_map_feasible(np.full(50, 11.0), 30)
    assert w1p.r2_adopt(3.79) and not w1p.r2_adopt(3.80) and not w1p.r2_adopt(float("nan"))


def test_anisotropic_psf_removes_its_own_distortion():
    # a round source of sigma 3 px convolved with an elongated PSF: isotropic deconvolution leaves
    # the PSF's direction in, the anisotropic one with the PSF's chi returns a round source
    psf_s, chi = 1.2, 0.3 * np.exp(2j * np.deg2rad(130.0))
    pee, pnn, pen = es.psf_moments(psf_s, chi)
    # the PSF's minor-axis variance is psf_sigma^2 and its distortion is chi
    lam = np.linalg.eigvalsh(np.array([[pee, pen], [pen, pnn]]))
    assert np.isclose(lam[0], psf_s**2) and np.isclose((lam[1] - lam[0]) / lam.sum(), abs(chi))
    qee, qnn, qen = 9 + pee, 9 + pnn, pen
    tr, det = qnn + qee, qnn * qee - qen**2
    lam1 = tr / 2 + np.sqrt(tr**2 / 4 - det)
    lam2 = tr / 2 - np.sqrt(tr**2 / 4 - det)
    pa = np.degrees(0.5 * np.arctan2(2 * qen, qnn - qee))
    a, b = np.sqrt(lam1), np.sqrt(lam2)
    iso = es.intrinsic_ellipticity([a], [b], [pa], psf_s)[0]
    ani = es.intrinsic_ellipticity([a], [b], [pa], psf_s, chi)[0]
    assert abs(iso) > 0.01
    assert np.isclose(np.mod(0.5 * np.degrees(np.angle(iso)), 180), 130.0)
    assert abs(ani) < 1e-9
    # moments_chi of the PSF-only object recovers chi
    got = w1p.moments_chi([np.sqrt(lam[1])], [np.sqrt(lam[0])], [130.0])[0]
    assert np.isclose(got, chi)


def test_lensed_cov_uses_the_same_psf():
    # an unlensed (identity Jacobian) image of a PSF-convolved source keeps its moments for any chi
    ir = w1p.ish.ir
    chi = 0.2 * np.exp(2j * np.deg2rad(40.0))
    a, b, pa = np.array([3.0]), np.array([2.0]), np.array([70.0])
    for c in (0j, chi):
        cov = ir.lensed_cov(a, b, pa, [1.0], [0.0], [1.0], [1.0], 1.2, c)
        sm, ell, p2 = ir.cov_shape(cov)
        assert np.allclose([sm[0], ell[0]], [3.0, 1 / 3]) and np.isclose(np.mod(p2[0], 180), 70.0)


def test_default_is_unchanged():
    a, b, pa = np.array([3.0, 2.0]), np.array([1.5, 1.9]), np.array([10.0, 100.0])
    assert np.allclose(
        es.intrinsic_ellipticity(a, b, pa, 1.2), es.intrinsic_ellipticity(a, b, pa, 1.2, 0j)
    )


def test_idw_map_and_star_selection():
    xs, ys = np.array([0.0, 10.0]), np.array([0.0, 0.0])
    chi = np.array([0.1 + 0j, 0.3 + 0j])
    m = w1p.idw_map(xs, ys, chi, np.array([0.0, 5.0]), np.array([0.0, 0.0]), k=2)
    assert np.isclose(m[0].real, 0.1, atol=1e-3) and np.isclose(m[1].real, 0.2)
    assert np.isclose(w1p.idw_map(xs, ys, chi, np.array([9.0]), np.array([0.0]), k=1)[0], 0.3)
    shapes = Table(
        {
            "ra": [10.0, 10.001, 10.002],
            "dec": [0.0, 0.0, 0.0],
            "snr": [100.0, 100.0, 20.0],
            "mag": [20.0, 18.0, 21.0],
        }
    )
    gaia = Table({"ra": [10.0, 10.001, 10.002], "dec": [0.0, 0.0, 0.0], "mag": [20.5, 19.0, 20.0]})
    assert list(w1p.select_stars(shapes, gaia)) == [0]  # row 1 too bright, row 2 low S/N


def test_r3_rule():
    old = {m: {"upper_limit_deg2": 10.0} for m in w1p.R3_MASSES}
    lower = {m: {"upper_limit_deg2": 9.0} for m in w1p.R3_MASSES}
    assert w1p.r3_adopt(lower, old)
    mixed = dict(lower, **{"2e+13": {"upper_limit_deg2": 10.0}})
    assert not w1p.r3_adopt(mixed, old)
    assert not w1p.r3_adopt(dict(lower, **{"8e+12": {"upper_limit_deg2": None}}), old)
