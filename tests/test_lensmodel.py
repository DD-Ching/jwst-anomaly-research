"""Offline tests for the Lenstool dPIE lens-model stage (D-024)."""

from __future__ import annotations

import numpy as np
import pytest
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

from jwst_anomaly import lensmodel
from jwst_anomaly.lensmodel import DPIE, LensModel, UnsupportedModelError

RA0, DEC0 = 110.826989, -73.454723
COSMO = FlatLambdaCDM(H0=70.0, Om0=0.3)

PAR = """\
runmode
    reference 3 {ra} {dec}
    end
cosmologie
    H0 70.0
    omegaM 0.3
    omegaX 0.7
    end
image
    z_m_limit 1 4.0 0 2.17 0.0 0.0
    z_m_limit 2 6.0 1 1.0 3.0 0.1
    sigposArcsec 0.44
    end
potentiel O1
    profil 81
    x_centre 0.5
    y_centre -1.0
    ellipticite 0.3
    angle_pos 30.0
    core_radius_kpc 10.0
    cut_radius 900.0
    v_disp 800.0
    z_lens 0.39
    end
fini
"""


def _write(tmp_path, text, name="best.par"):
    path = tmp_path / name
    path.write_text(text)
    return path


def _dpie(**kw):
    base = dict(
        name="p",
        x=0.3,
        y=-0.7,
        ellipticity=0.35,
        angle_pos=30.0,
        r_core=0.8,
        r_cut=150.0,
        v_disp=700.0,
        z_lens=0.39,
    )
    base.update(kw)
    return DPIE(**base)


def test_parse_lenstool_par_french_keywords(tmp_path):
    parsed = lensmodel.parse_lenstool_par(_write(tmp_path, PAR.format(ra=RA0, dec=DEC0)))
    assert parsed["reference"] == {"mode": 3, "ra": RA0, "dec": DEC0}
    assert parsed["cosmology"] == {"H0": 70.0, "Om0": 0.3}
    assert parsed["z_m_limit"] == {"4": 2.17}  # system 6 is free (flag 1), so not fixed
    assert parsed["sigpos_arcsec"] == 0.44
    (pot,) = parsed["potentials"]
    assert pot["name"] == "O1" and pot["profile"] == 81
    assert pot["ellipticity"] == 0.3 and pot["core_radius_kpc"] == 10.0
    assert np.isnan(pot["core_radius"])
    assert len(parsed["sha256"]) == 64


def test_kpc_radius_uses_model_cosmology(tmp_path):
    model = LensModel.from_par(_write(tmp_path, PAR.format(ra=RA0, dec=DEC0)))
    (comp,) = model.components
    kpc_per_arcsec = COSMO.kpc_proper_per_arcmin(0.39).value / 60.0
    assert comp.r_core == pytest.approx(10.0 / kpc_per_arcsec)
    assert comp.r_cut == 900.0


@pytest.mark.parametrize(
    "old, new, err",
    [
        ("profil 81", "profil 12", UnsupportedModelError),
        ("reference 3", "reference 2", UnsupportedModelError),
        ("omegaX 0.7", "omegaX 0.6", UnsupportedModelError),
        ("    end\nfini", "fini", ValueError),
    ],
)
def test_parse_rejects_unsupported(tmp_path, old, new, err):
    text = PAR.format(ra=RA0, dec=DEC0).replace(old, new)
    with pytest.raises(err):
        lensmodel.parse_lenstool_par(_write(tmp_path, text))


def test_potfile_is_not_expanded(tmp_path):
    text = PAR.format(ra=RA0, dec=DEC0).replace(
        "fini", "potfile\n    filein 3 gal.cat\n    end\nfini"
    )
    with pytest.raises(UnsupportedModelError, match="potfile"):
        LensModel.from_par(_write(tmp_path, text))


def test_frame_roundtrip_and_orientation():
    model = LensModel([_dpie()], RA0, DEC0, COSMO)
    x, y = model.to_frame(RA0 - 1.0 / 3600, DEC0 + 2.0 / 3600)  # West and North of the reference
    assert x == pytest.approx(np.cos(np.deg2rad(DEC0)))
    assert y == pytest.approx(2.0)
    ra, dec = model.to_sky(x, y)
    assert ra == pytest.approx(RA0 - 1.0 / 3600, abs=1e-12)
    assert dec == pytest.approx(DEC0 + 2.0 / 3600, abs=1e-12)
    # +y (North) -> PA 0; +x (West) and -x (East) -> PA 90 (axes are mod 180).
    np.testing.assert_allclose(LensModel.frame_angle_to_pa([90.0, 0.0, 180.0]), [0.0, 90.0, 90.0])


def test_axis_offset_deg():
    np.testing.assert_allclose(
        lensmodel.axis_offset_deg([10, 170, 0, 45], [170, 10, 90, 45]), [20, 20, 90, 0]
    )


def test_hessian_trace_equals_analytic_kappa():
    comp = _dpie()
    rng = np.random.default_rng(1)
    x, y = rng.uniform(-30, 30, 500), rng.uniform(-30, 30, 500)
    hxx, _, hyy = comp.hessian(x, y)
    np.testing.assert_allclose(0.5 * (hxx + hyy), comp.kappa(x, y), rtol=1e-7)


def test_deflection_gradient_matches_hessian():
    comp = _dpie()
    rng = np.random.default_rng(2)
    x, y = rng.uniform(-20, 20, 200), rng.uniform(-20, 20, 200)
    h = 1e-4
    axp, ayp = comp.deflection(x + h, y)
    axm, aym = comp.deflection(x - h, y)
    bxp, byp = comp.deflection(x, y + h)
    bxm, bym = comp.deflection(x, y - h)
    hxx, hxy, hyy = comp.hessian(x, y)
    np.testing.assert_allclose((axp - axm) / (2 * h), hxx, rtol=1e-5, atol=1e-8)
    np.testing.assert_allclose((byp - bym) / (2 * h), hyy, rtol=1e-5, atol=1e-8)
    np.testing.assert_allclose((bxp - bxm) / (2 * h), hxy, rtol=1e-5, atol=1e-8)
    np.testing.assert_allclose((ayp - aym) / (2 * h), hxy, rtol=1e-5, atol=1e-8)  # curl-free


def test_isothermal_limit():
    # A round dPIE with a tiny core and a huge cut is an SIS: |alpha| = b0 and kappa = b0 / (2 r).
    comp = _dpie(x=0.0, y=0.0, ellipticity=0.0, r_core=1e-4, r_cut=1e5)
    ax, ay = comp.deflection(np.array([5.0, 0.0]), np.array([0.0, -8.0]))
    np.testing.assert_allclose(np.hypot(ax, ay), comp.b0, rtol=2e-3)
    np.testing.assert_allclose(
        [ax[0], ay[1]], [comp.b0, -comp.b0], rtol=2e-3
    )  # along theta: beta = theta - alpha
    np.testing.assert_allclose(
        comp.kappa(np.array([5.0]), np.array([0.0])), comp.b0 / 10.0, rtol=2e-3
    )


def test_dls_ds_and_evaluate():
    model = LensModel(
        [_dpie(x=0.0, y=0.0, ellipticity=0.0, r_core=1e-3, r_cut=1e4)], RA0, DEC0, COSMO
    )
    d_a = COSMO.angular_diameter_distance(2.0).value
    d_ls = (
        COSMO.comoving_transverse_distance(2.0) - COSMO.comoving_transverse_distance(0.39)
    ).value / 3.0
    expected = d_ls / d_a
    np.testing.assert_allclose(model.dls_ds([0.2, 2.0]), [0.0, expected])
    theta_e = model.components[0].b0 * expected
    # North of the centre, outside the Einstein radius: tangential arcs run East-West (PA 90).
    t = model.evaluate(RA0, DEC0 + 2 * theta_e / 3600, 2.0)
    assert t.meta["provenance"] == "model_prediction"
    assert t["tangential_pa"][0] == pytest.approx(90.0, abs=1e-6)
    assert t["magnification"][0] > 1
    # Inside theta_E / 2, kappa > 1 and the stretch is radial (North-South, PA 0) with odd parity.
    t = model.evaluate(RA0, DEC0 + 0.3 * theta_e / 3600, 2.0)
    assert t["kappa"][0] > 1
    assert min(t["tangential_pa"][0], 180 - t["tangential_pa"][0]) == pytest.approx(0.0, abs=1e-6)
    # In front of the lens nothing is lensed.
    t = model.evaluate(RA0, DEC0 + 1 / 3600, 0.2)
    assert t["kappa"][0] == 0 and t["magnification"][0] == pytest.approx(1.0)


def test_backtrace_consistent_system(tmp_path):
    # SIS-like lens; a source at beta = (0.5, 0) arcsec has images at x = beta +- theta_E.
    model = LensModel(
        [_dpie(x=0.0, y=0.0, ellipticity=0.0, r_core=1e-4, r_cut=1e5)], RA0, DEC0, COSMO
    )
    theta_e = model.components[0].b0 * model.dls_ds(2.0)
    xs = np.array([0.5 + theta_e, 0.5 - theta_e])
    ra, dec = model.to_sky(xs, np.zeros(2))
    arcs = tmp_path / "arcs.dat"
    arcs.write_text(
        "#REFERENCE 0\n"
        + "".join(
            f"1.{k + 1} {r:.9f} {d:.9f} 0.1 0.1 0 2.0 0\n"
            for k, (r, d) in enumerate(zip(ra, dec, strict=True))
        )
        + f"2.1 {RA0:.9f} {DEC0 + 3 / 3600:.9f} 0.1 0.1 0 0 0\n"  # free redshift, not in z_m_limit
    )
    images = lensmodel.load_lenstool_images(arcs)
    assert list(images["system"]) == ["1", "1", "2"]
    assert images.meta["provenance"] == "observed"
    bt = lensmodel.backtrace_images(model, images, z_m_limit={})
    assert np.all(bt["dbeta_arcsec"][:2] < 2e-3 * theta_e)
    assert np.all(bt["dtheta_arcsec"][:2] < 1e-2)
    assert np.isnan(bt["z_used"][2]) and np.isnan(bt["dtheta_arcsec"][2])
    assert bt["magnification"][0] > 0 > bt["magnification"][1]  # opposite parities


def test_load_images_rejects_relative_coordinates(tmp_path):
    arcs = tmp_path / "arcs.dat"
    arcs.write_text("#REFERENCE 3 110.8 -73.4\n1.1 1.0 2.0 0.1 0.1 0 2.0 0\n")
    with pytest.raises(UnsupportedModelError):
        lensmodel.load_lenstool_images(arcs)


def test_multi_plane_rejected():
    with pytest.raises(UnsupportedModelError):
        LensModel([_dpie(), _dpie(z_lens=0.5)], RA0, DEC0, COSMO)


def test_subset_tables_keep_meta():
    model = LensModel([_dpie()], RA0, DEC0, COSMO)
    for method in (
        model.kappa,
        model.gamma,
        model.magnification,
        model.deflection,
        model.tangential_pa,
    ):
        t = method([RA0, RA0 + 1e-3], [DEC0, DEC0], 2.0)
        assert isinstance(t, Table) and len(t) == 2
        assert t.meta["provenance"] == "model_prediction"


def test_system_key_keeps_non_integer_ids():
    assert [lensmodel.system_key(v) for v in ("4", "4.0", 4.0, "4.10", "c2", "1a")] == [
        "4",
        "4",
        "4",
        "4.10",
        "c2",
        "1a",
    ]


def test_image_redshifts_prefer_z_m_limit(tmp_path):
    arcs = tmp_path / "arcs.dat"
    arcs.write_text(
        "1.1 110.8 -73.4 0.1 0.1 0 1.5 0\n"  # catalogued z, overridden by the model's fixed value
        "2.1 110.8 -73.4 0.1 0.1 0 2.5 0\n"  # catalogued z, not in z_m_limit
        "3.1 110.8 -73.4 0.1 0.1 0 0 0\n"  # neither
        "4.10.1 110.8 -73.4 0.1 0.1 0 0 0\n"  # system "4.10", not "4.1"
    )
    images = lensmodel.load_lenstool_images(arcs)
    assert list(images["system"]) == ["1", "2", "3", "4.10"]
    z = lensmodel.image_redshifts(images, {"1": 1.7, "4.10": 3.0, "4.1": 9.0})
    np.testing.assert_allclose(z, [1.7, 2.5, np.nan, 3.0])


def test_unknown_redshift_stays_unknown():
    model = LensModel([_dpie()], RA0, DEC0, COSMO)
    assert np.isnan(model.dls_ds(np.nan))
    t = model.evaluate([RA0, RA0], [DEC0 + 5 / 3600] * 2, [np.nan, 2.0])
    assert np.isnan(t["kappa"][0]) and np.isnan(t["magnification"][0])
    assert np.isfinite(t["kappa"][1])


def test_unknown_potential_keyword_is_refused(tmp_path):
    text = PAR.format(ra=RA0, dec=DEC0).replace(
        "    v_disp 800.0", "    v_disp 800.0\n    ellip_pot 0.1"
    )
    with pytest.raises(UnsupportedModelError, match="ellip_pot"):
        lensmodel.parse_lenstool_par(_write(tmp_path, text))


def test_inconsistent_core_radius_is_refused(tmp_path):
    text = PAR.format(ra=RA0, dec=DEC0).replace(
        "core_radius_kpc 10.0", "core_radius_kpc 10.0\n    core_radius 9.0"
    )
    with pytest.raises(UnsupportedModelError, match="disagrees"):
        LensModel.from_par(_write(tmp_path, text))


def test_latin1_comments_are_accepted(tmp_path):
    text = "# mod\u00e8le de r\u00e9f\u00e9rence\n" + PAR.format(ra=RA0, dec=DEC0)
    path = tmp_path / "best.par"
    path.write_bytes(text.encode("latin-1"))
    assert len(lensmodel.parse_lenstool_par(path)["potentials"]) == 1


def test_outputs_meet_schema_contracts(tmp_path):
    from jwst_anomaly import schema

    model = LensModel([_dpie()], RA0, DEC0, COSMO)
    schema.validate(model.evaluate(RA0, DEC0, 2.0), schema.LENS_PREDICTION_COLUMNS)
    arcs = tmp_path / "arcs.dat"
    arcs.write_text("1.1 110.8 -73.4 0.1 0.1 0 2.0 0\n1.2 110.81 -73.41 0.1 0.1 0 2.0 0\n")
    bt = lensmodel.backtrace_images(model, lensmodel.load_lenstool_images(arcs), {})
    schema.validate(bt, schema.BACKTRACE_COLUMNS)


def test_find_images_sis_pair():
    # SIS-like lens: a source at beta = (1, 0) has images at x = beta +- theta_E.
    model = LensModel(
        [_dpie(x=0.0, y=0.0, ellipticity=0.0, r_core=1e-4, r_cut=1e5)], RA0, DEC0, COSMO
    )
    theta_e = model.components[0].b0 * float(model.dls_ds(2.0))
    grid = lensmodel.DeflectionGrid.compute(model, half_width=2 * theta_e, step=0.2)
    imgs = lensmodel.find_images(model, grid, 1.0, 0.0, 2.0)
    bright = imgs[np.abs(imgs["magnification"]) > 0.1]
    assert len(bright) == 2
    np.testing.assert_allclose(sorted(bright["x"]), [1.0 - theta_e, 1.0 + theta_e], atol=2e-3)
    np.testing.assert_allclose(bright["y"], 0.0, atol=2e-3)
    # SIS magnifications: mu = 1 / (1 - theta_E / |x|), positive outside, negative inside the ring
    xs = np.asarray(bright["x"])
    np.testing.assert_allclose(bright["magnification"], 1 / (1 - theta_e / np.abs(xs)), rtol=5e-3)
    assert imgs.meta["provenance"] == "model_prediction"


def test_find_images_elliptical_quad():
    # An elliptical lens with the source near the centre gives four bright images (plus a faint
    # central one for a cored profile).
    model = LensModel(
        [_dpie(x=0.0, y=0.0, ellipticity=0.4, angle_pos=20.0, r_core=0.5, r_cut=300.0)],
        RA0,
        DEC0,
        COSMO,
    )
    grid = lensmodel.DeflectionGrid.compute(model, half_width=60.0, step=0.25)
    imgs = lensmodel.find_images(model, grid, 0.3, -0.2, 2.0)
    bright = imgs[np.abs(imgs["magnification"]) > 0.5]
    assert len(bright) == 4
    assert (bright["magnification"] > 0).sum() == 2  # two of each parity
    assert np.all(imgs["residual_arcsec"] < 1e-5)


def test_deflection_grid_cache(tmp_path):
    model = LensModel([_dpie()], RA0, DEC0, COSMO, sha256="abc")
    path = tmp_path / "grid.npz"
    g1 = lensmodel.DeflectionGrid.cached(model, path, half_width=5.0, step=0.5)
    g2 = lensmodel.DeflectionGrid.cached(model, path, half_width=5.0, step=0.5)
    np.testing.assert_array_equal(g1.alpha_x, g2.alpha_x)
    other = LensModel([_dpie(v_disp=500.0)], RA0, DEC0, COSMO, sha256="def")
    g3 = lensmodel.DeflectionGrid.cached(other, path, half_width=5.0, step=0.5)
    assert not np.allclose(g3.alpha_x, g1.alpha_x)  # a different model recomputes


def test_image_family_lenstool_ids():
    ids = ("23a", "1.1a", "A200.1a", "4.1", "4.10", "4.0", "4", "7", "c2")
    assert [lensmodel.image_family(i) for i in ids] == [
        "23",
        "1.1",
        "A200.1",
        "4",
        "4",
        "4",
        "4",
        "7",
        "c2",
    ]


def test_z_m_limit_letter_ids_and_shared_redshift(tmp_path):
    # Bergamini+2023b: image ids with letters, and one line fixing three families' redshift.
    text = PAR.format(ra=RA0, dec=DEC0).replace(
        "    z_m_limit 1 4.0 0 2.17 0.0 0.0\n",
        "    z_m_limit 1 4.0 0 2.17 0.0 0.0\n"
        "    z_m_limit 1 301.1a 0 5.29 0.0 0.0\n"
        "    z_m_limit 1 A200.1a B200.2a 0 7.39 0.0 0.0\n"
        "    z_m_limit 1 9a 4 1.0 3.0 0.1\n"  # boundary flag 4 and parabolic -n: free, not fixed
        "    z_m_limit 1 10a -2 1.0 3.0 0.1\n",
    )
    zml = lensmodel.parse_lenstool_par(_write(tmp_path, text))["z_m_limit"]
    assert zml == {"4": 2.17, "301.1": 5.29, "A200.1": 7.39, "B200.2": 7.39}


def test_letter_suffixed_images_and_error_column(tmp_path):
    arcs = tmp_path / "arcs.dat"
    arcs.write_text(
        "#REFERENCE 0\n"
        "23a 110.8 -73.4 0.62 0.62 0 2.19 25\n"
        "23b 110.8 -73.4 1.24 1.24 0 2.19 25\n"
        "1.1a 110.8 -73.4 0.38 0.38 0 1.69 25\n"
    )
    images = lensmodel.load_lenstool_images(arcs)
    assert list(images["system"]) == ["23", "23", "1.1"]
    np.testing.assert_allclose(images["a"], [0.62, 1.24, 0.38])


def test_kpc_radius_tolerates_six_decimal_rounding(tmp_path):
    # Bergamini+2023b writes core_radius 0.000021 next to core_radius_kpc 0.000097 (2% apart).
    per_arcsec = COSMO.kpc_proper_per_arcmin(0.39).value / 60.0
    core = round(9.7e-5 / per_arcsec, 6)
    text = PAR.format(ra=RA0, dec=DEC0).replace(
        "core_radius_kpc 10.0", f"core_radius_kpc 0.000097\n    core_radius {core:.6f}"
    )
    model = LensModel.from_par(_write(tmp_path, text))
    assert model.components[0].r_core == pytest.approx(core)


def test_imageplane_residuals_sis_pair_and_offset(tmp_path):
    # Exact images of an SIS give zero residual; moving one catalogued image by 0.3" gives ~0.15"
    # for both (the mean source shifts), never more than the offset.
    model = LensModel(
        [_dpie(x=0.0, y=0.0, ellipticity=0.0, r_core=1e-4, r_cut=1e5)], RA0, DEC0, COSMO
    )
    theta_e = model.components[0].b0 * float(model.dls_ds(2.0))
    grid = lensmodel.DeflectionGrid.compute(model, half_width=2 * theta_e, step=0.2)
    for shift, expect in ((0.0, 0.0), (0.3, 0.15)):
        xs = np.array([1.0 + theta_e + shift, 1.0 - theta_e])
        ra, dec = model.to_sky(xs, np.zeros(2))
        arcs = _write(
            tmp_path,
            "".join(
                f"1.{k + 1} {r:.9f} {d:.9f} 0.1 0.1 0 2.0 0\n"
                for k, (r, d) in enumerate(zip(ra, dec, strict=True))
            ),
            "arcs.dat",
        )
        bt = lensmodel.backtrace_images(model, lensmodel.load_lenstool_images(arcs), {})
        ip = lensmodel.imageplane_residuals(model, grid, bt)
        np.testing.assert_allclose(ip["dtheta_arcsec"], expect, atol=5e-3)
        assert np.all(ip["n_predicted"] >= 2) and ip.meta["provenance"] == "model_prediction"


def test_z_m_limit_malformed_flag_is_refused(tmp_path):
    text = PAR.format(ra=RA0, dec=DEC0).replace(
        "    z_m_limit 1 4.0 0 2.17 0.0 0.0\n", "    z_m_limit 1 A200.1a B200.2a 0 7.39\n"
    )
    with pytest.raises(ValueError, match="z_m_limit"):
        lensmodel.parse_lenstool_par(_write(tmp_path, text))


def _find_images_reference(model, grid, beta_x, beta_y, z_s, newton_steps=12, tol=1e-5):
    """The original one-seed-at-a-time solver over every grid triangle (no pre-filter)."""
    s = float(model.dls_ds(z_s))
    g = grid.x
    bx = g[None, :] - s * grid.alpha_x
    by = g[:, None] - s * grid.alpha_y
    n = len(g) - 1
    seeds = []
    for (a0, a1), (b0, b1), (c0, c1) in (((0, 0), (0, 1), (1, 0)), ((1, 1), (1, 0), (0, 1))):
        ax_, ay_ = bx[a0 : a0 + n, a1 : a1 + n], by[a0 : a0 + n, a1 : a1 + n]
        bx_, by_ = bx[b0 : b0 + n, b1 : b1 + n], by[b0 : b0 + n, b1 : b1 + n]
        cx_, cy_ = bx[c0 : c0 + n, c1 : c1 + n], by[c0 : c0 + n, c1 : c1 + n]
        det = (by_ - cy_) * (ax_ - cx_) + (cx_ - bx_) * (ay_ - cy_)
        with np.errstate(divide="ignore", invalid="ignore"):
            l1 = ((by_ - cy_) * (beta_x - cx_) + (cx_ - bx_) * (beta_y - cy_)) / det
            l2 = ((cy_ - ay_) * (beta_x - cx_) + (ax_ - cx_) * (beta_y - cy_)) / det
        for iy, ix in zip(*np.nonzero((l1 >= 0) & (l2 >= 0) & (1 - l1 - l2 >= 0)), strict=True):
            w1, w2 = l1[iy, ix], l2[iy, ix]
            w3 = 1 - w1 - w2
            seeds.append(
                (
                    w1 * g[ix + a1] + w2 * g[ix + b1] + w3 * g[ix + c1],
                    w1 * g[iy + a0] + w2 * g[iy + b0] + w3 * g[iy + c0],
                )
            )
    found = []
    for x0, y0 in seeds:
        x, y = np.array([x0]), np.array([y0])
        for _ in range(newton_steps):
            f = model.fields_xy(x, y)
            rx, ry = beta_x - (x - s * f["alpha_x"]), beta_y - (y - s * f["alpha_y"])
            a11, a12, a22 = 1 - s * f["psi_xx"], -s * f["psi_xy"], 1 - s * f["psi_yy"]
            det = a11 * a22 - a12 * a12
            x, y = x + (a22 * rx - a12 * ry) / det, y + (-a12 * rx + a11 * ry) / det
        f = model.fields_xy(x, y)
        if np.hypot(beta_x - (x - s * f["alpha_x"]), beta_y - (y - s * f["alpha_y"]))[0] < tol:
            if all(np.hypot(x[0] - u, y[0] - v) > 0.05 for u, v in found):
                found.append((float(x[0]), float(y[0])))
    return found


@pytest.mark.parametrize("beta", [(0.3, -0.2), (1.0, 0.0), (0.0, 0.0), (6.0, 2.5)])
def test_find_images_matches_the_unfiltered_solver(beta):
    # includes sources on the symmetry axes, where mapped grid values can equal beta exactly
    model = LensModel(
        [_dpie(x=0.0, y=0.0, ellipticity=0.4, angle_pos=20.0, r_core=0.5, r_cut=300.0)],
        RA0,
        DEC0,
        COSMO,
    )
    grid = lensmodel.DeflectionGrid.compute(model, half_width=40.0, step=0.5)
    imgs = lensmodel.find_images(model, grid, *beta, 2.0)
    ref = _find_images_reference(model, grid, *beta, 2.0)
    assert len(imgs) == len(ref)
    got = sorted(zip(imgs["x"], imgs["y"], strict=True))
    np.testing.assert_allclose(got, sorted(ref), atol=1e-9)


def test_map_lens_model_reproduces_the_analytic_model():
    from astropy.wcs import WCS

    model = LensModel(
        [_dpie(x=0.0, y=0.0, ellipticity=0.3, angle_pos=20.0, r_core=1.0, r_cut=300.0)],
        RA0,
        DEC0,
        COSMO,
    )
    n, pix = 801, 0.1  # 80" maps, north-up, east-left
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [RA0, DEC0]
    w.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    w.wcs.cdelt = [-pix / 3600, pix / 3600]
    jj, ii = np.mgrid[0:n, 0:n]
    x, y = model.to_frame(*w.pixel_to_world_values(ii, jj))  # each pixel's true sky position
    ax, ay = model.deflection_xy(x, y)  # +x = West along +i, as the class documents
    mm = lensmodel.MapLensModel(ax, ay, w, model.z_lens, COSMO, source="synthetic maps")
    assert abs(mm.ra0 - RA0) < 1e-9 and abs(mm.dec0 - DEC0) < 1e-9
    px, py = np.array([5.3, -12.1, 20.7]), np.array([-7.4, 3.3, 15.2])
    fa, fm = model.fields_xy(px, py), mm.fields_xy(px, py)
    for k in ("alpha_x", "alpha_y"):
        np.testing.assert_allclose(fm[k], fa[k], atol=2e-3)
    for k in ("psi_xx", "psi_xy", "psi_yy"):
        np.testing.assert_allclose(fm[k], fa[k], atol=5e-3)
    assert np.isnan(mm.fields_xy(100.0, 0.0)["alpha_x"])  # off the maps
    grid = lensmodel.DeflectionGrid.compute(mm, half_width=35.0, step=0.25)
    got = lensmodel.find_images(mm, grid, 0.6, -0.3, 2.0)
    ref = lensmodel.find_images(
        model, lensmodel.DeflectionGrid.compute(model, 35.0, 0.25), 0.6, -0.3, 2.0
    )
    bright_g = got[np.abs(got["magnification"]) > 0.5]
    bright_r = ref[np.abs(ref["magnification"]) > 0.5]
    assert len(bright_g) == len(bright_r) >= 2
    np.testing.assert_allclose(
        sorted(zip(bright_g["x"], bright_g["y"], strict=True)),
        sorted(zip(bright_r["x"], bright_r["y"], strict=True)),
        atol=0.01,
    )


def test_read_z_m_limit_without_parsing_the_model(tmp_path):
    path = _write(
        tmp_path,
        "potfile\n  filein 3 members.cat\n  end\nimage\n"
        "  z_m_limit 1 5.0 0 3.948 0.000 0.0000\n  z_m_limit 2 7a 0 1.954 0.0 0.0\n"
        "  z_m_limit 3 9.0 1 1.0 3.0 0.1\n  end\n",
        "params.txt",
    )
    assert lensmodel.read_z_m_limit(path) == {"5": 3.948, "7": 1.954}


def test_load_images_accepts_parenthesised_redshifts(tmp_path):
    path = _write(tmp_path, f"#REFERENCE 0\n1.1 {RA0} {DEC0} 0.1 0.1 0.0 (2.16) 0.0\n", "arcs.dat")
    assert lensmodel.load_lenstool_images(path)["z"][0] == 2.16


def test_find_images_refines_cells_on_a_fold():
    # D-040: a merging pair straddling a critical curve inside one coarse grid cell (here about
    # 0.1" apart in a 0.5" cell) is missed by the plain triangle scan; refinement recovers it
    model = LensModel(
        [_dpie(x=0.0, y=0.0, ellipticity=0.4, angle_pos=20.0, r_core=0.5, r_cut=300.0)],
        RA0,
        DEC0,
        COSMO,
    )
    s = float(model.dls_ds(2.0))

    def det(x):
        f = model.fields_xy(np.atleast_1d(x), np.zeros(1))
        a11, a12, a22 = 1 - s * f["psi_xx"], -s * f["psi_xy"], 1 - s * f["psi_yy"]
        return float((a11 * a22 - a12 * a12)[0])

    xs = np.linspace(1.0, 30.0, 2901)  # the outer (tangential) critical curve on the +x axis
    signs = np.sign([det(x) for x in xs])
    i = int(np.flatnonzero(signs[:-1] != signs[1:])[-1])
    x_img = xs[i + 1] + 0.02  # just outside the curve; its partner is about 0.1" away
    f = model.fields_xy(np.array([x_img]), np.zeros(1))
    beta = (float(x_img - s * f["alpha_x"][0]), float(-s * f["alpha_y"][0]))
    coarse = lensmodel.DeflectionGrid.compute(model, half_width=40.0, step=0.5)
    plain = lensmodel.find_images(model, coarse, *beta, 2.0, refine_arcsec=0)
    refined = lensmodel.find_images(model, coarse, *beta, 2.0)

    def near(t):
        return t[np.hypot(t["x"] - x_img, t["y"]) < 0.5]

    assert len(near(plain)) < 2
    pair = near(refined)
    assert len(pair) == 2  # one image on each side of the curve, opposite parity
    assert np.sign(pair["magnification"][0]) != np.sign(pair["magnification"][1])
    assert np.min(np.hypot(pair["x"] - x_img, pair["y"])) < 1e-3


def _cluster_and_galaxy():
    cluster = _dpie(
        name="cl", x=0.0, y=0.0, ellipticity=0.4, angle_pos=20.0, r_core=0.5, r_cut=300.0
    )
    galaxy = _dpie(name="gal", x=9.0, y=4.0, ellipticity=0.1, r_core=0.05, r_cut=20.0, v_disp=180.0)
    return LensModel([cluster, galaxy], RA0, DEC0, COSMO)


def test_split_planes_without_moves_matches_the_single_plane_model():
    model = _cluster_and_galaxy()
    multi = model.split_planes({})
    assert multi.z_planes == (0.39,)
    grid = lensmodel.DeflectionGrid.compute(model, half_width=40.0, step=0.25)
    mgrid = lensmodel.DeflectionGrid.compute(multi, half_width=40.0, step=0.25)
    a = lensmodel.find_images(model, grid, 0.3, -0.2, 2.0)
    b = lensmodel.find_images(multi, mgrid, 0.3, -0.2, 2.0)
    assert len(a) == len(b) >= 4
    np.testing.assert_allclose(b["x"], a["x"], atol=1e-6)
    np.testing.assert_allclose(b["magnification"], a["magnification"], rtol=1e-6)
    ra, dec = model.to_sky(a["x"][:3], a["y"][:3])
    images = Table(
        {
            "image_id": ["1.1", "1.2", "1.3"],
            "system": ["1"] * 3,
            "ra": ra,
            "dec": dec,
            "z": [2.0] * 3,
        }
    )
    bt_a = lensmodel.backtrace_images(model, images, {"1": 2.0})
    bt_b = lensmodel.backtrace_images(multi, images, {"1": 2.0})
    for col in ("beta_x", "beta_y", "magnification", "dtheta_arcsec"):
        np.testing.assert_allclose(bt_b[col], bt_a[col], rtol=1e-6, atol=1e-8)


def test_multiplane_jacobian_matches_finite_differences():
    multi = _cluster_and_galaxy().split_planes({"gal": 0.2}, v_disp={"gal": 250.0})
    assert multi.z_planes == (0.2, 0.39)
    x, y = np.array([8.0, 10.5, -3.0]), np.array([3.0, 5.5, 12.0])
    m = multi.lens_map(x, y, 3.0)
    h = 1e-5
    bxp, byp = multi.source_points(x + h, y, 3.0)
    bxm, bym = multi.source_points(x - h, y, 3.0)
    np.testing.assert_allclose(m["a11"], (bxp - bxm) / (2 * h), atol=1e-6)
    np.testing.assert_allclose(m["a21"], (byp - bym) / (2 * h), atol=1e-6)
    bxp, byp = multi.source_points(x, y + h, 3.0)
    bxm, bym = multi.source_points(x, y - h, 3.0)
    np.testing.assert_allclose(m["a12"], (bxp - bxm) / (2 * h), atol=1e-6)
    np.testing.assert_allclose(m["a22"], (byp - bym) / (2 * h), atol=1e-6)
    assert np.max(np.abs(m["a12"] - m["a21"])) > 1e-3  # two planes: A is not symmetric
    np.testing.assert_allclose(m["beta_x"], multi.source_points(x, y, 3.0)[0])


def test_multiplane_plane_behind_the_source_does_not_lens():
    model = _cluster_and_galaxy()
    multi = model.split_planes({"gal": 2.5})
    alone = LensModel([model.components[0]], RA0, DEC0, COSMO)
    x, y = np.linspace(-20, 20, 7), np.linspace(15, -15, 7)
    m, ref = multi.lens_map(x, y, 2.0), alone.lens_map(x, y, 2.0)
    for key in ("beta_x", "beta_y", "a11", "a12", "a21", "a22"):
        np.testing.assert_allclose(m[key], ref[key], atol=1e-12)
    grid = lensmodel.DeflectionGrid.compute(multi, half_width=30.0, step=0.25)
    imgs = lensmodel.find_images(multi, grid, 0.3, -0.2, 2.0)
    assert np.all(imgs["residual_arcsec"] < 1e-5)
    assert imgs.meta["lens_planes"][1] == {"z_lens": 2.5, "n_potentials": 1, "potentials": ["gal"]}


def test_split_planes_rejects_unknown_potentials():
    model = _cluster_and_galaxy()
    with pytest.raises(ValueError, match="unknown or unmoved"):
        model.split_planes({"nope": 0.2})
    with pytest.raises(ValueError, match="unknown or unmoved"):
        model.split_planes({"gal": 0.2}, v_disp={"cl": 100.0})


def test_multiplane_grid_serves_every_source_redshift(tmp_path):
    model = LensModel(list(_cluster_and_galaxy().components), RA0, DEC0, COSMO, sha256="abc")
    multi = model.split_planes({"gal": 0.2}, v_disp={"gal": 250.0})
    grid = lensmodel.DeflectionGrid.cached(multi, tmp_path / "g.npz", half_width=10.0, step=0.5)
    assert grid.alpha_x.shape == (2, 41, 41)
    again = lensmodel.DeflectionGrid.cached(multi, tmp_path / "g.npz", half_width=10.0, step=0.5)
    np.testing.assert_array_equal(again.alpha_x, grid.alpha_x)
    gx, gy = np.meshgrid(grid.x, grid.x)
    for z_s in (0.3, 1.0, 4.0):  # one grid, any source redshift (the rays do not depend on z_s)
        bx, by = multi.source_grid(grid, z_s)
        ex, ey = multi.source_points(gx, gy, z_s)
        np.testing.assert_allclose(bx, ex, atol=1e-10)
        np.testing.assert_allclose(by, ey, atol=1e-10)
    single = lensmodel.DeflectionGrid.compute(model, half_width=10.0, step=0.5)
    with pytest.raises(ValueError, match="another model"):
        lensmodel.find_images(multi, single, 0.3, -0.2, 2.0)
    assert multi.planes[0].sha256 not in ("", model.sha256)  # a plane is not the whole file


def test_multiplane_shift_frame_leaves_the_caller_planes_alone():
    model = _cluster_and_galaxy()
    other = LensModel([_dpie(name="g2", z_lens=0.6, v_disp=150.0)], RA0, DEC0, COSMO)
    multi = lensmodel.MultiPlaneLensModel([model, other])
    multi.shift_frame(0.5, -0.2)
    assert (model.ra0, model.dec0) == (RA0, DEC0)
    assert multi.planes[0].ra0 != RA0 and multi.ra0 == multi.planes[1].ra0


def test_multiplane_identity_tracks_the_cosmology_and_grids_check_their_model():
    comps = list(_cluster_and_galaxy().components)
    a = LensModel(comps, RA0, DEC0, COSMO, sha256="f").split_planes({"gal": 0.2})
    b = LensModel(comps, RA0, DEC0, FlatLambdaCDM(H0=50.0, Om0=0.5), sha256="f")
    b = b.split_planes({"gal": 0.2})
    assert a.sha256 and a.sha256 != b.sha256
    assert LensModel(comps, RA0, DEC0, COSMO).split_planes({"gal": 0.2}).sha256 == ""
    before = a.sha256
    a.shift_frame(0.1, 0.0)
    assert a.sha256 != before
    named = LensModel(comps[1:], RA0, DEC0, FlatLambdaCDM(H0=70.0, Om0=0.3, name="x"))
    lensmodel.MultiPlaneLensModel(
        [LensModel(comps[:1], RA0, DEC0, COSMO), named.split_planes({"gal": 0.2}).planes[0]]
    )
    single = LensModel(comps, RA0, DEC0, COSMO)  # in memory: no sha256, so only the shape check
    grid = lensmodel.DeflectionGrid.compute(a, half_width=5.0, step=0.5)
    with pytest.raises(ValueError, match="single-plane model"):
        lensmodel.find_images(single, grid, 0.3, -0.2, 2.0)
    with pytest.raises(AttributeError, match="z_planes"):
        _ = a.z_lens
