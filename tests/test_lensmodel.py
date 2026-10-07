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
