import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from m0308_grid_posthoc import classify, n_bits, n_exposures  # noqa: E402

COV_OK = {"f300m_e1_n_exp": 3, "f300m_e2_n_exp": 4, "f410m_e1_n_exp": 3, "f410m_e2_n_exp": 4}


def _star(sep, pa):
    return {"groups": {"0": {"sep_arcsec": sep, "pa_deg": pa}}}


def test_n_bits_and_all_con_planes():
    assert n_bits(0) == 0
    assert n_bits(0b1011) == 3
    con = np.zeros((2, 1, 1), dtype=np.int32)
    con[0, 0, 0], con[1, 0, 0] = 0b1, 0b11
    assert n_exposures(con, 0, 0) == 3


def test_n_bits_counts_bit_31_of_int32():
    assert n_bits(np.int32(-1)) == 32  # all 32 exposures of a plane
    assert n_bits(np.int32(-1431658496)) == 10  # an Abell 1689 F356W CON plane


def test_star_wings():
    assert classify({"0": COV_OK}, _star(1.5, 120.1)) == {"0": ["star_wings"]}


def test_fixed_offset_needs_separation_and_angle():
    assert classify({"0": COV_OK}, _star(8.30, 120.5)) == {"0": ["fixed_offset_star_image"]}
    assert classify({"0": COV_OK}, _star(8.30, 122.0)) == {"0": ["unexplained"]}
    assert classify({"0": COV_OK}, _star(8.60, 120.1)) == {"0": ["unexplained"]}


def test_every_matching_label_is_kept():
    cov = dict(COV_OK, f410m_e1_n_exp=1)
    assert classify({"0": cov}, _star(8.33, 120.1)) == {
        "0": ["fixed_offset_star_image", "single_exposure"]
    }


def test_single_exposure_and_off_footprint_are_distinct():
    assert classify({"0": dict(COV_OK, f410m_e2_n_exp=1)}, _star(30.0, 10.0)) == {
        "0": ["single_exposure"]
    }
    assert classify({"0": dict(COV_OK, f300m_e1_n_exp=0)}, _star(30.0, 10.0)) == {
        "0": ["off_footprint"]
    }
    assert classify({"0": COV_OK}, {"groups": {}}) == {"0": ["unexplained"]}
