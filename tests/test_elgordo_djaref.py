"""Unit tests of the El Gordo reference-catalogue rules (PLAN_1.md rule 2, PLAN_2.md rule 2)."""

import io
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from astropy.table import Table  # noqa: E402
from elgordo_djaref import (  # noqa: E402
    SHAPE_COLS,
    flux_unit_jy,
    grid_rotation_deg,
    header_cards,
    reference_inputs,
    to_cat_dja,
)


def _blob(names) -> bytes:
    cols = [fits.Column(name=n, format="E", array=np.zeros(3)) for n in names]
    buf = io.BytesIO()
    hdus = [fits.PrimaryHDU(), fits.BinTableHDU.from_columns(cols)]
    hdus[1].header["APER_0"] = 9.0
    fits.HDUList(hdus).writeto(buf)
    return buf.getvalue()


def test_header_cards_reads_primary_and_table():
    heads = header_cards(_blob(["id", "ra"]))
    assert [h.get("XTENSION") for h in heads] == [None, "BINTABLE"]
    assert heads[1]["NAXIS2"] == 3


def test_passes_with_shape_and_total_flux_columns():
    r = reference_inputs(_blob([*SHAPE_COLS, "f200w_flux_aper_0", "f200w_tot_1"]), 100)
    assert r["passes"] and r["rows"] == 3 and r["total_flux_columns"] == ["f200w_tot_1"]


def test_fails_without_total_flux_column_and_lists_tot_corr():
    r = reference_inputs(_blob([*SHAPE_COLS, "f200w_flux_aper_0", "f200w_tot_corr"]), 100)
    assert not r["passes"] and r["shape_columns_present"]
    assert r["post_hoc_tot_corr_columns"] == ["f200w_tot_corr"]
    assert r["post_hoc_detection_aperture_pix"] == {"APER_0": 9.0}


def test_fails_without_shape_columns():
    r = reference_inputs(_blob(["x_image", "y_image", "f200w_tot_1"]), 100)
    assert not r["passes"] and not r["shape_columns_present"]


def test_no_table_header():
    buf = io.BytesIO()
    fits.PrimaryHDU().writeto(buf)
    assert reference_inputs(buf.getvalue(), 1) == {
        "reachable": True,
        "bytes": 1,
        "table_header": False,
        "passes": False,
    }


def _unit_table(zp):
    f = np.array([1.0, 10.0, 100.0, -1.0])
    return Table({"flux_auto": f, "mag_auto": zp - 2.5 * np.log10(np.abs(f))})


def test_flux_unit_from_mag_auto():
    assert flux_unit_jy(_unit_table(23.9))["unit_jy"] == 1e-6
    assert flux_unit_jy(_unit_table(28.9))["unit_jy"] == 1e-8
    r = flux_unit_jy(_unit_table(25.0))
    assert not r["passes"] and r["unit_jy"] is None
    assert flux_unit_jy(_unit_table(23.9))["n_rows"] == 3  # the negative flux is left out


def _grid_table(rot_deg=0.0, parity=1.0, pix=0.04):
    rng = np.random.default_rng(1)
    x, y = rng.uniform(0, 4000, 200), rng.uniform(0, 4000, 200)
    c, s = np.cos(np.radians(rot_deg)), np.sin(np.radians(rot_deg))
    # east-left north-up for parity 1: xi = -x, eta = y, then rotated by rot_deg toward east
    xi = (-parity * x * c + y * s) * pix
    eta = (parity * x * s + y * c) * pix
    dec0 = -49.25
    return Table(
        {
            "x_image": x,
            "y_image": y,
            "ra": 15.72 + xi / 3600.0 / np.cos(np.radians(dec0)),
            "dec": dec0 + eta / 3600.0,
        }
    )


def test_grid_north_up_east_left_passes():
    r = grid_rotation_deg(_grid_table())
    assert r["passes"]
    assert abs(r["rotation_deg"]) < 0.05
    assert abs(r["pixel_scale_arcsec"] - 0.04) < 1e-3


def test_grid_rotated_or_flipped_fails():
    r = grid_rotation_deg(_grid_table(rot_deg=5.0))
    assert not r["passes"] and abs(r["rotation_deg"] - 5.0) < 0.1
    assert not grid_rotation_deg(_grid_table(parity=-1.0))["passes"]


def test_to_cat_dja_total_flux_and_drops():
    n = 4
    t = Table(
        {
            "id": np.arange(1, n + 1),
            "ra": np.full(n, 15.7),
            "dec": np.full(n, -49.25),
            "a_image": np.full(n, 4.0),
            "b_image": np.full(n, 2.0),
            "theta_image": np.zeros(n),
            "flux_auto": np.array([3.0, 3.0, 3.0, 3.0]),
            "flux_aper_1": np.array([1.5, 1.5, 0.0, 1.5]),
            "f200w_flux_aper_1": np.array([2.0, -1.0, 2.0, 2.0]),
            "f200w_fluxerr_aper_1": np.full(n, 0.5),
            "f200w_tot_corr": np.full(n, 1.1),
            "f200w_mask_aper_1": np.array([0.0, 0.0, 0.0, 100.0]),
        },
        meta={"APER_1": 12.5},
    )
    rc = to_cat_dja(t, 1e-6, 0.04)
    # row 2: negative flux; row 3: detection aperture flux 0; row 4: masked > 0.2 x aperture area
    assert list(rc["label"]) == [1]
    np.testing.assert_allclose(rc["flux"], [2.0 * 3.0 / 1.5 * 1.1 * 1e-6])
    np.testing.assert_allclose(rc["snr"], [4.0])
    np.testing.assert_allclose(rc["a"], [0.16])
    np.testing.assert_allclose(rc["ell"], [0.5])
    assert rc["ext"].all()
