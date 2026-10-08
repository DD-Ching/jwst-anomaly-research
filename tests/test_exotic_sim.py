"""Tests for jwst_anomaly.exotic_sim against statements in the verified papers (D-047). Offline."""

from __future__ import annotations

import numpy as np
import pytest

from jwst_anomaly import exotic_sim as es
from jwst_anomaly import schema


def _paczynski(u):
    """Standard point-lens total magnification (u² + 2) / (u sqrt(u² + 4))."""
    return (u**2 + 2) / (u * np.sqrt(u**2 + 4))


def test_n1_attractive_is_the_standard_point_lens():
    u = np.array([0.05, 0.3, 1.0, 2.0, 5.0])
    np.testing.assert_allclose(es.total_magnification(u, 1, 1), _paczynski(u), rtol=1e-10)
    # Kitamura, Nakajima & Asada 2013 (arXiv:1211.0379): n = 1 never demagnifies.
    assert np.all(es.total_magnification(np.linspace(0.01, 20, 500), 1, 1) > 1)


def test_numeric_roots_agree_with_closed_form_at_n1():
    b = np.array([0.1, 1.0, 2.5, 3.0, 7.0])
    for sign in (1, -1):
        closed = es._roots_n1(b, sign)
        numeric = es._roots_numeric(b, 1.0, sign)
        np.testing.assert_allclose(numeric, closed, rtol=1e-9, equal_nan=True)


def test_negative_point_mass_umbra_caustic_and_magnification():
    """Safonova, Torres & Romero 2002 (gr-qc/0105070): no images for u < 2, caustic at u = 2 and
    total magnification (u² - 2) / (u sqrt(u² - 4)) outside it."""
    assert es.caustic_beta(1) == pytest.approx(2.0)
    assert np.all(es.total_magnification([0.0, 0.5, 1.99], 1, -1) == 0)
    u = np.array([2.01, 2.5, 3.0, 6.0])
    np.testing.assert_allclose(
        es.total_magnification(u, 1, -1), (u**2 - 2) / (u * np.sqrt(u**2 - 4)), rtol=1e-10
    )
    # magnification diverges at the caustic
    assert es.total_magnification(2.0 + 1e-8, 1, -1)[0] > 1e3
    # both images on the source's side and radially elongated (Izumi et al. 2013)
    sol = es.solve_images([2.2, 3.0], 1, -1)
    assert np.all(sol["x"] > 0)
    assert np.all(np.abs(sol["lam_r"]) < np.abs(sol["lam_t"]))


def test_inverted_microlensing_light_curve_dip_and_spikes():
    t = np.linspace(-3, 3, 1201)
    a = es.light_curve(t, 0.0, 1.0, 1.0, n=1, sign=-1)
    centre = np.abs(t) < 1.0  # beta < 2 for |t| < sqrt(3)
    assert np.all(a[centre] == 0)
    i_peak = np.argmax(a)
    assert (
        a[i_peak] > 5 and abs(abs(t[i_peak]) - np.sqrt(3)) < 0.02
    )  # spike at the caustic crossing
    assert abs(a[0] - 1) < 0.1  # baseline far outside
    # a finite source caps the spike
    af = es.light_curve(t, 0.0, 1.0, 1.0, n=1, sign=-1, rho=0.1)
    assert 1.5 < af.max() < a.max()


def test_ellis_wormhole_statements_of_abe_2010():
    """Abe 2010 (arXiv:1009.6084): inner image at -0.618 and -0.532 for beta = 2 and 3, carrying
    A2/A = 0.034 and 0.013; gutters about 4 % deep."""
    sol = es.solve_images([2.0, 3.0], 2, 1)
    np.testing.assert_allclose(sol["x"][:, 1], [-0.618, -0.532], atol=1e-3)
    frac = np.abs(sol["mu"][:, 1]) / np.nansum(np.abs(sol["mu"]), axis=1)
    np.testing.assert_allclose(frac, [0.034, 0.013], atol=1e-3)
    beta = np.linspace(0.5, 4, 3501)
    depth = 1 - es.total_magnification(beta, 2, 1).min()
    assert 0.03 < depth < 0.05


def test_demagnification_onsets_of_kitamura_2013():
    """KNA13: n = 10 total magnification < 1 for beta > 0.187 (numerical); n = 10 depletes up to
    ~60 % near beta ~ 0.7; n = 3 depletes ~10 % near beta ~ 1.1 (we find 14 %)."""
    beta = np.linspace(0.05, 3, 2951)
    a10 = es.total_magnification(beta, 10, 1)
    assert beta[np.argmax(a10 < 1)] == pytest.approx(0.187, abs=0.003)
    assert 1 - a10.min() == pytest.approx(0.6, abs=0.03) and beta[a10.argmin()] == pytest.approx(
        0.7, abs=0.05
    )
    a3 = es.total_magnification(beta, 3, 1)
    assert 0.05 < 1 - a3.min() < 0.2 and beta[a3.argmin()] == pytest.approx(1.1, abs=0.1)
    assert es.demagnification_threshold(10) == pytest.approx(2 / 11)


def test_convergence_and_shear_signs_of_izumi_2013():
    x = np.array([0.8, 1.5, 3.0])
    k1, _ = es.convergence_shear(x, 1, 1)
    np.testing.assert_allclose(k1, 0, atol=1e-12)
    k2, g2 = es.convergence_shear(x, 2, 1)  # Ellis: negative convergence, tangential shear
    np.testing.assert_allclose(k2, -0.5 / x**3)
    np.testing.assert_allclose(g2, -1.5 / x**3)
    k2r, g2r = es.convergence_shear(x, 2, -1)  # repulsive n = 2: positive kappa, radial shear
    assert np.all(k2r > 0) and np.all(g2r > 0)


def test_ellis_einstein_radius_matches_abe():
    """Abe 2010: R_E = (pi/4 D_L D_LS / D_S a²)^(1/3), theta_E = R_E / D_L."""
    a, dl, ds = 1e8, 1.2e20, 2.4e20  # m (a = 1e5 km; bulge-like distances)
    dls = ds - dl
    re = (np.pi / 4 * dl * dls / ds * a**2) ** (1 / 3)
    th = es.einstein_radius(es.eps_bar_ellis(a), 2, dl, ds, dls)
    assert th == pytest.approx(re / dl, rel=1e-12)
    # n = 1 reduces to the usual sqrt(4GM/c² D_LS/(D_L D_S))
    th1 = es.einstein_radius(es.eps_bar_point_mass(1.0), 1, dl, ds, dls)
    assert th1 == pytest.approx(np.sqrt(4 * es.G_SI * es.MSUN_KG / es.C_SI**2 * dls / (dl * ds)))
    assert es.deflection_integral(1) == pytest.approx(1.0)
    assert es.deflection_integral(2) == pytest.approx(np.pi / 4)


def test_finite_source_matches_point_source_and_disk_formula():
    # uniform disk centred on a point lens: A = sqrt(rho² + 4) / rho (image of the disk edge
    # spans radii (sqrt(rho²+4) ± rho)/2, so the image area over the disk area is that ratio)
    rho = 0.5
    assert es.finite_source_magnification(0.0, rho, 1, 1)[0] == pytest.approx(
        np.sqrt(rho**2 + 4) / rho, rel=0.01
    )
    b = np.array([1.5, 3.0])
    np.testing.assert_allclose(
        es.finite_source_magnification(b, 1e-4, 2, 1), es.total_magnification(b, 2, 1), rtol=1e-3
    )


def test_count_ratio_deficit_inside_einstein_radius():
    counts = lambda s: s**-1.5  # noqa: E731  (steep counts, ASSUMPTION for the test)
    r = es.count_ratio(np.array([0.3, 3.0]), counts, 1.0, n=1, sign=-1)
    assert r[0] < 0.2 and r[1] > 1  # deficit of demagnified inner images, mild excess outside


@pytest.mark.parametrize(("n", "sign"), [(1, 1), (2, 1), (1, -1), (3, -1)])
def test_inject_images_round_trip(n, sign):
    """Mapping each injected image back through the lens equation recovers its source."""
    rng = np.random.default_rng(1)
    theta_e = 1.7
    sx, sy = rng.uniform(-6, 6, 200), rng.uniform(-6, 6, 200)
    tab = es.inject_images(sx, sy, theta_e, n=n, sign=sign, flux=2.0)
    assert tab.meta["provenance"] == schema.Provenance.SIMULATED.value
    assert "source" in tab.meta
    x, y = np.asarray(tab["dx"]) / theta_e, np.asarray(tab["dy"]) / theta_e
    r = np.hypot(x, y)
    defl = sign / r ** (n + 1)  # beta_vec = theta_vec (1 - sign / |theta|^(n+1))
    bx, by = x * (1 - defl) * theta_e, y * (1 - defl) * theta_e
    sid = np.asarray(tab["source_id"])
    np.testing.assert_allclose(bx, sx[sid], atol=1e-8)
    np.testing.assert_allclose(by, sy[sid], atol=1e-8)
    np.testing.assert_allclose(tab["flux"], 2.0 * np.abs(tab["mu"]))
    beta = np.hypot(sx, sy) / theta_e
    if sign == 1:
        assert len(tab) == 2 * sx.size and not np.any(tab["radial"])
    else:
        assert len(tab) == 2 * np.sum(beta >= es.caustic_beta(n))
        # outer images (all images for n = 1) are radial: long axis along the lens→image direction.
        # For n > 1 the inner image at x < 1 has |lam_r| > lam_t, i.e. it is tangential in shape.
        rad = (tab["image"] == 0) | (n == 1)
        assert np.all(tab["radial"][rad])
        pa = np.degrees(np.arctan2(tab["dx"], tab["dy"])) % 180
        np.testing.assert_allclose(tab["pa_deg"][rad], pa[rad])


def test_inject_light_curve_round_trip_and_blend():
    t = np.linspace(-50, 50, 41)
    lc = es.inject_light_curve(t, 10.0, t0=0.0, t_e=20.0, u0=0.5, n=1, sign=-1, blend=0.6)
    assert lc.meta["provenance"] == schema.Provenance.SIMULATED.value
    a = es.light_curve(t, 0.0, 20.0, 0.5, 1, -1)
    np.testing.assert_allclose(lc["flux"], 10.0 * (0.6 * a + 0.4))
    np.testing.assert_allclose(
        np.asarray(lc["flux"])[np.abs(t) < 20], 4.0
    )  # umbra leaves the blend
    np.testing.assert_allclose(lc["beta"], np.sqrt(0.25 + (t / 20.0) ** 2))


def test_bad_parameters_raise():
    with pytest.raises(ValueError):
        es.solve_images(1.0, n=0, sign=1)
    with pytest.raises(ValueError):
        es.solve_images(1.0, n=1, sign=2)
    with pytest.raises(ValueError):
        es.solve_images(-1.0)
