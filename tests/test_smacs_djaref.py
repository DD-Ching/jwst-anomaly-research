"""Unit tests of the SMACS 0723 reference-catalogue rules (PLAN_1.md 0, 2; PLAN_2.md D2, D3)."""

import io
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from smacs_djaref import STAGE2_COLS, fixed_z_systems, split_z_edge, stage2_inputs  # noqa: E402

PAR = """runmode
\tz_m_limit  1 4.0 1 1.0 7.0 0.1000
\tz_m_limit  2 6 1 1.0 7.0 0.1000
\tz_m_limit  3 12.0 1 1.0 7.0 0.1000
\tend
"""


def _blob(names, cards=("APER_1",)) -> bytes:
    cols = [fits.Column(name=n, format="E", array=np.zeros(3)) for n in names]
    hdus = [fits.PrimaryHDU(), fits.BinTableHDU.from_columns(cols)]
    for k in cards:
        hdus[1].header[k] = 7.0
    buf = io.BytesIO()
    fits.HDUList(hdus).writeto(buf)
    return buf.getvalue()


def test_fixed_z_systems_drops_every_z_m_limit_system():
    assert fixed_z_systems(PAR, ["1", "4", "6", "7", "12", "2", "1"]) == ["1", "2", "7"]


def test_fixed_z_systems_ignores_other_keywords():
    assert fixed_z_systems("z_m 1 4.0\n", [4, 5]) == ["4", "5"]


def test_stage2_inputs_pass_with_every_column_and_card():
    r = stage2_inputs(_blob(STAGE2_COLS), 10)
    assert r["passes"] and r["missing"] == [] and r["rows"] == 3 and r["aper_1_pix"] == 7.0


def test_stage2_inputs_fail_on_a_missing_column_or_card():
    r = stage2_inputs(_blob([c for c in STAGE2_COLS if c != "f200w_tot_corr"]), 10)
    assert not r["passes"] and r["missing"] == ["f200w_tot_corr"]
    r = stage2_inputs(_blob(STAGE2_COLS, cards=()), 10)
    assert not r["passes"] and r["missing"] == ["APER_1"]


def test_stage2_inputs_without_a_table_fail():
    buf = io.BytesIO()
    fits.PrimaryHDU().writeto(buf)
    assert not stage2_inputs(buf.getvalue(), 10)["passes"]


def test_split_z_edge_reports_edge_families_and_labels_the_rest():
    flags = ["8@110.8,-73.4", "3@110.9,-73.5", "17@110.7,-73.45"]
    r = split_z_edge(flags, {"8": 10.0, "3": 1.99, "17": 2.05}, ["1", "3"])
    assert r["flags_z_edge"] == ["8@110.8,-73.4"]
    assert r["flags"] == [
        {"flag": "3@110.9,-73.5", "z": "fixed_z"},
        {"flag": "17@110.7,-73.45", "z": "fitted_z"},
    ]


def test_split_z_edge_keeps_a_family_just_below_the_edge():
    r = split_z_edge(["9@1,2"], {"9": 9.99}, [])
    assert r == {"flags": [{"flag": "9@1,2", "z": "fitted_z"}], "flags_z_edge": []}
