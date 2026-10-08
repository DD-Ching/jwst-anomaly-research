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


def _ray_shoot(beta, rho, n, sign, xr, yr, pix):
    """Independent inverse ray shooting: fraction of image-plane rays landing in the source disk."""
    x = np.arange(xr[0], xr[1], pix) + pix / 2
    y = np.arange(yr[0], yr[1], pix) + pix / 2
    xx, yy = np.meshgrid(x, y)
    f = 1 - sign / np.hypot(xx, yy) ** (n + 1)  # vector lens equation (Izumi et al. 2013)
    inside = (xx * f - beta) ** 2 + (yy * f) ** 2 < rho**2
    return inside.sum() * pix**2 / (np.pi * rho**2)


@pytest.mark.parametrize(
    ("beta", "rho", "n", "sign", "xr", "yr", "pix"),
    [
        (2.0, 0.1, 1, -1, (0.3, 1.8), (-0.1, 0.1), 0.0015),  # disk straddles the caustic
        (2.07, 0.1, 1, -1, (0.3, 1.8), (-0.1, 0.1), 0.0015),  # near the spike peak
        (2.2, 0.1, 1, -1, (0.3, 1.8), (-0.1, 0.1), 0.0015),
        (1.95, 0.01, 1, -1, (0.8, 1.25), (-0.01, 0.01), 0.0002),  # tiny source, spike
        (0.2, 0.3, 2, 1, (-1.6, 1.9), (-1.5, 1.5), 0.006),  # disk covers the Ellis lens
    ],
)
def test_finite_source_agrees_with_inverse_ray_shooting(beta, rho, n, sign, xr, yr, pix):
    """Tolerance 1 %: the ray-shooting pixel noise at these grids is ~0.1-0.5 %."""
    expect = _ray_shoot(beta, rho, n, sign, xr, yr, pix)
    assert es.finite_source_magnification(beta, rho, n, sign)[0] == pytest.approx(expect, rel=0.01)


def test_converged_caustic_spike_heights():
    """Peak disk-averaged magnification of a negative point mass (simulated; docs W3, D-047)."""
    for rho, peak in ((0.01, 7.013), (0.1, 2.353), (0.3, 1.529)):
        b = np.linspace(1.95, 2.0 + 2 * rho, 4001)
        a = es.finite_source_magnification(b, rho, 1, -1)
        assert a.max() == pytest.approx(peak, abs=0.002)
        # converged: a 4x finer quadrature changes the peak by < 0.1 %
        i = a.argmax()
        fine = es.finite_source_magnification(b[i], rho, 1, -1, n_nodes=192)[0]
        assert fine == pytest.approx(a[i], rel=1e-3)


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
    """KNA13 (arXiv:1211.0379): n = 10 onset 0.187 numerically against 2/(n+1) = 0.182 at leading
    order; n = 10 depletes ~60 % at beta ~ 0.7. For n = 3 the text says ~10 % at beta ~ 1.1, but
    their Fig. 2c (n = 3, beta_0 = 0.1) bottoms out at A ≈ 0.865 ± 0.01 (read from the figure
    pixels), i.e. 13-14 %: the text rounds. We test the code's own 14.3 %."""
    assert es.demagnification_onset(10) == pytest.approx(0.187, abs=0.001)
    assert es.demagnification_onset_approx(10) == pytest.approx(0.182, abs=0.001)
    assert es.demagnification_onset(2) == pytest.approx(1.1111, abs=1e-3)
    assert es.demagnification_onset(3) == pytest.approx(0.6429, abs=1e-3)
    assert np.isinf(es.demagnification_onset(1)) and np.isinf(es.demagnification_onset(2, -1))
    a_on = es.total_magnification([es.demagnification_onset(3)], 3, 1)[0]
    assert a_on == pytest.approx(1.0, abs=1e-9)
    beta = np.linspace(0.05, 3, 5901)
    a10 = es.total_magnification(beta, 10, 1)
    assert 1 - a10.min() == pytest.approx(0.587, abs=0.002)
    assert beta[a10.argmin()] == pytest.approx(0.70, abs=0.01)
    a3 = es.total_magnification(beta, 3, 1)
    assert 1 - a3.min() == pytest.approx(0.143, abs=0.002)
    assert beta[a3.argmin()] == pytest.approx(1.12, abs=0.01)


def test_convergence_and_shear_signs_of_izumi_2013():
    x = np.array([0.8, 1.5, 3.0])
    k1, _ = es.convergence_shear(x, 1, 1)
    np.testing.assert_allclose(k1, 0, atol=1e-12)
    k2, g2 = es.convergence_shear(x, 2, 1)  # Ellis: negative convergence, tangential shear
    np.testing.assert_allclose(k2, -0.5 / x**3)
    np.testing.assert_allclose(g2, -1.5 / x**3)
    k2r, g2r = es.convergence_shear(x, 2, -1)  # repulsive n = 2: positive kappa, radial shear
    assert np.all(k2r > 0) and np.all(g2r > 0)


def test_ellis_einstein_radius_matches_abe_tables():
    """Abe 2010 Tables 1-2 (bulge: D_S = 8 kpc, D_L = 4 kpc; bound v_T = 220 km/s):
    a = 1e3 km -> R_E 3.64e7 km, 0.061 mas; a = 1e5 km -> R_E 7.85e8 km, 1.31 mas, t_E 41.3 d."""
    kpc = 1e3 * es.PC_M
    dl, ds = 4 * kpc, 8 * kpc
    for a_km, re_km, th_mas, te_day in ((1e3, 3.64e7, 0.061, 1.92), (1e5, 7.85e8, 1.31, 41.3)):
        th = es.einstein_radius(es.eps_bar_ellis(a_km * 1e3), 2, dl, ds, ds - dl)
        assert th * dl / 1e3 == pytest.approx(re_km, rel=0.005)
        assert np.degrees(th) * 3.6e6 == pytest.approx(th_mas, rel=0.01)
        assert th * dl / 220e3 / 86400 == pytest.approx(te_day, rel=0.005)
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


def test_physical_scales_quoted_in_docs():
    """docs/exotic_lensing.md scales; ASSUMPTIONs: Planck18, z_l = 0.4, z_s = 2, v = 1000 km/s."""
    from astropy.cosmology import Planck18

    dl = Planck18.angular_diameter_distance(0.4).to_value("m")
    ds = Planck18.angular_diameter_distance(2.0).to_value("m")
    dls = Planck18.angular_diameter_distance(0.4, 2.0).to_value("m")
    arcsec = np.pi / 180 / 3600
    th12 = es.einstein_radius(es.eps_bar_point_mass(1e12), 1, dl, ds, dls) / arcsec
    assert th12 == pytest.approx(2.2, abs=0.05)
    for a_pc, th in ((0.1, 0.033), (1, 0.15), (10, 0.72), (100, 3.3)):
        th_e = es.einstein_radius(es.eps_bar_ellis(a_pc * es.PC_M), 2, dl, ds, dls) / arcsec
        assert th_e == pytest.approx(th, rel=0.05)
    th1 = es.einstein_radius(es.eps_bar_point_mass(1.0), 1, dl, ds, dls)
    t_e_yr = th1 * dl / 1e6 / (365.25 * 86400)
    assert t_e_yr == pytest.approx(12, rel=0.05)


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
        # For n > 1 the inner image is tangential in shape where x^(n+1) < (n - 1)/2.
        x_in = np.hypot(tab["dx"], tab["dy"]) / theta_e
        rad = (tab["image"] == 0) | (x_in ** (n + 1) >= (n - 1) / 2)
        assert np.all(~np.asarray(tab["radial"])[~rad])
        assert np.all(tab["radial"][rad])
        pa = np.degrees(np.arctan2(tab["dx"], tab["dy"])) % 180
        np.testing.assert_allclose(tab["pa_deg"][rad], pa[rad])


@pytest.mark.parametrize("n", [2, 3])
def test_repulsive_inner_image_shape_switch(n):
    """Tangential (|lam_r| > lam_t) iff x^(n+1) < (n - 1)/2: x < 0.794 (n = 2), x < 1 (n = 3)."""
    x_sw = ((n - 1) / 2) ** (1 / (n + 1))
    for x, radial in ((0.97 * x_sw, False), (1.03 * x_sw, True)):
        beta = x + x**-n  # source radius whose inner image sits at x
        tab = es.inject_images([beta], [0.0], 1.0, n=n, sign=-1)
        inner = tab[tab["image"] == 1]
        assert inner["dx"][0] == pytest.approx(x, rel=1e-9)
        assert bool(inner["radial"][0]) is radial


def test_n1_inner_image_precise_at_large_beta():
    x = es.solve_images([1e8], 1, 1)["x"][0]
    assert x[1] == pytest.approx(-1e-8, rel=1e-9)  # (beta - sqrt(beta² + 4))/2 would cancel to 0


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
    with pytest.raises(ValueError):
        es.finite_source_magnification(2.0, -0.1)
    with pytest.raises(ValueError):
        es.light_curve([0.0], 0.0, 1.0, 0.5, rho=-0.1)
    for t_e in (0.0, -1.0):
        with pytest.raises(ValueError):
            es.impact_track([0.0], 0.0, t_e, 0.5)
        with pytest.raises(ValueError):
            es.inject_light_curve([0.0], 1.0, 0.0, t_e, 0.5)
    with pytest.raises(ValueError):
        es.inject_images([1.0, 2.0], [0.0, 0.0], 1.0, source_id=[7])
    with pytest.raises(ValueError):
        es.inject_light_curve([0.0], 1.0, 0.0, 1.0, 0.5, blend=1.5)
