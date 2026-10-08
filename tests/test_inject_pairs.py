"""Tests for scripts/inject_pairs.py (synthetic standard-layout catalogue only; offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from jwst_anomaly import exotic_sim

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("inject_pairs", _DIR / "inject_pairs.py")
ip = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ip)
op = ip.op

BANDS = ["F435W", "F606W", "F814W", "F090W", "F115W", "F150W", "F200W", "F277W", "F356W", "F444W"]
WAVE = np.array([0.435, 0.606, 0.814, 0.90, 1.15, 1.50, 2.00, 2.77, 3.56, 4.44])
RA0, DEC0 = 64.0, -24.0


def standard_catalogue(seed: int = 5) -> Table:
    """A sparse grid of faint sources plus one bright source (row 0) at the origin."""
    rng = np.random.default_rng(seed)
    gx, gy = np.meshgrid(np.arange(-40.0, 41.0, 8.0), np.arange(-40.0, 41.0, 8.0))
    x = np.r_[0.0, gx.ravel() + 4.0]
    y = np.r_[0.0, gy.ravel() + 4.0]
    n = len(x)
    slope = rng.uniform(-1.5, 1.5, n)
    flux = (WAVE[None, :] ** slope[:, None]) * np.r_[20000.0, np.full(n - 1, 200.0)][:, None]
    err = np.full_like(flux, 5.0)
    t = Table(
        {
            "src_id": np.arange(1, n + 1),
            "ra": RA0 + x / 3600.0 / np.cos(np.deg2rad(DEC0)),
            "dec": DEC0 + y / 3600.0,
        }
    )
    for c in ("x_pix", "y_pix", "xmin", "xmax", "ymin", "ymax", "pa_deg"):
        t[c] = np.zeros(n)
    t["x_pix"], t["y_pix"] = -x / 0.04, y / 0.04
    t["xmin"], t["xmax"] = t["x_pix"] - 3, t["x_pix"] + 3
    t["ymin"], t["ymax"] = t["y_pix"] - 3, t["y_pix"] + 3
    t["ap_radius"] = np.full(n, 0.15)
    for c in ("deblend", "bcg", "pointsrc"):
        t[c] = np.zeros(n, bool)
    t["use"] = np.ones(n, bool)
    t["z_low"], t["z16"], t["z84"], t["z_best"] = 1.6, 1.8, 2.2, 2.0
    t["z_spec"], t["mu"] = -1.0, 1.0
    t["flux"], t["err"] = flux, err
    t.meta.update(bands=list(BANDS), format="test", pixscale=0.04, source="synthetic")
    return t


def test_mass_and_throat_conversions_invert_einstein_radius():
    z_l = 0.4
    m = ip.theta_e_to_mass(1.0, z_l)
    assert np.isclose(ip.theta_e_to_mass(2.0, z_l) / m, 4.0)
    d_l, d_s, d_ls = ip._distances(z_l)
    th = exotic_sim.einstein_radius(exotic_sim.eps_bar_point_mass(m), 1.0, d_l, d_s, d_ls)
    assert np.isclose(np.rad2deg(th) * 3600.0, 1.0)
    a_m = ip.theta_e_to_throat_pc(0.3, z_l) * exotic_sim.PC_M
    th2 = exotic_sim.einstein_radius(exotic_sim.eps_bar_ellis(a_m), 2.0, d_l, d_s, d_ls)
    assert np.isclose(np.rad2deg(th2) * 3600.0, 0.3)


def test_poisson_limits():
    assert np.isclose(ip.poisson_signal_ul(0, 0.0), ip.POISSON_UL_95, rtol=1e-3)
    assert ip.poisson_signal_ul(10, 10.0) > ip.poisson_signal_ul(10, 12.0)
    assert np.isinf(ip.surface_density_limit(0.0))


def test_point_lens_on_a_bright_isolated_source_is_recovered_as_an_orphan():
    field = ip.Field("synthetic", std=standard_catalogue(), z_lens=0.4)
    assert field.lensable[0] and field.selected[0]
    rng = np.random.default_rng(1)
    theta_e = 0.7
    lx, ly = field.x[0] - 0.3 * theta_e, field.y[0]  # beta 0.3 theta_E: images ~1.4" apart
    removed, images, src = ip.paint(field, lx, ly, theta_e, "point", rng)
    assert list(src) == [0] and 0 in removed and len(images) == 2
    assert np.all(np.asarray(images["src_id"]) < 0)
    mu = exotic_sim.solve_images(0.3, 1.0, 1)["mu"][0]
    ratio = np.nansum(images["flux"][0]) / np.nansum(images["flux"][1])
    assert np.isclose(ratio, abs(mu[0] / mu[1]), rtol=0.05)
    res = ip.evaluate(field, lx, ly, 3 * theta_e, removed, images)[0]
    assert res["n_selected"] == 2 and res["matched"] and res["cls"] == "orphan"
    assert res["recovered"] and 1.0 < res["sep"] < 2.0


def test_close_images_merge_into_one_row_and_are_not_recovered():
    field = ip.Field("synthetic", std=standard_catalogue(), z_lens=0.4)
    field.d_blend = 1.0  # force the merge rule
    removed, images, _ = ip.paint(
        field, field.x[0] - 0.1, field.y[0], 0.3, "point", np.random.default_rng(2), rows=[0]
    )
    assert len(images) == 1 and images["inj_img"][0] == 2
    mu = np.abs(exotic_sim.solve_images(0.1 / 0.3, 1.0, 1)["mu"][0]).sum()
    assert np.isclose(np.nansum(images["flux"][0]) / np.nansum(field.std["flux"][0]), mu, rtol=0.05)
    res = ip.evaluate(field, field.x[0] - 0.1, field.y[0], 0.9, removed, images)[0]
    assert not res.get("recovered", False)


def test_w1_removes_the_umbra_and_puts_both_images_on_the_source_side():
    field = ip.Field("synthetic", std=standard_catalogue(), z_lens=0.4)
    theta_e = 1.0
    # source at beta = 1.5 theta_E: inside the umbra, removed without images
    removed, images, src = ip.paint(
        field, field.x[0] - 1.5, field.y[0], theta_e, "W1", np.random.default_rng(3)
    )
    assert 0 in removed and images is None and len(src) == 0
    # beta = 2.5 theta_E: two images, both east of the lens (the source's side)
    removed, images, src = ip.paint(
        field, field.x[0] - 2.5, field.y[0], theta_e, "W1", np.random.default_rng(3)
    )
    assert list(src) == [0] and len(images) == 2
    x, _ = op.tangent_xy(images["ra"], images["dec"], field.ra0, field.dec0)
    assert np.all(x > field.x[0] - 2.5)


def test_per_lens_efficiency_is_bounded_and_counts_lensed_sources():
    field = ip.Field("synthetic", std=standard_catalogue(), z_lens=0.4)
    res = ip.per_lens(field, "point", 0.7, 30, np.random.default_rng(4))
    assert 0.0 <= res["efficiency"] <= res["frac_with_images"] <= 1.0
    src = ip.per_source(field, "point", 0.7, 5, np.random.default_rng(4))
    assert len(src) == 5 and set(src.colnames) >= {"mag", "recovered", "cls"}


def test_cli_refuses_to_combine_without_every_orphan_summary(tmp_path, capsys):
    (tmp_path / "orph" / "goodsn-dja").mkdir(parents=True)
    (tmp_path / "orph" / "goodsn-dja" / "summary.json").write_text("{}")
    argv = ["--combine-only", "--orphans", str(tmp_path / "orph"), "--out", str(tmp_path)]
    with pytest.raises(SystemExit) as e:
        ip.main([*argv, "--fields", "goodsn-dja", "macs0416-ncf"])
    assert e.value.code == 2
    assert "macs0416-ncf" in capsys.readouterr().err
    with pytest.raises(SystemExit):  # cluster fields have their own lens redshift
        ip.main([*argv, "--fields", "macs1149"])


def test_combine_rejects_fields_at_another_lens_redshift():
    with pytest.raises(ValueError, match="macs1149"):
        ip.combine([{"field": "macs1149", "z_lens": 0.543, "runs": {}}], {})


def test_mass_conversion_matches_inject_radial():
    spec = importlib.util.spec_from_file_location("inject_radial", _DIR / "inject_radial.py")
    ir = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ir)
    te = np.array([0.3, 1.5, 6.0])
    assert np.allclose(ip.theta_e_to_mass(te, 0.4), ir.theta_e_to_mass(te, 0.4), rtol=1e-10)
