"""Unit tests of the El Gordo stage-1 reference-input check (PLAN_1.md rule 2)."""

import io
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from elgordo_djaref import SHAPE_COLS, header_cards, reference_inputs  # noqa: E402


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
