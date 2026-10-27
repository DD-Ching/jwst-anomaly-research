import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from m0308_grid_posthoc import classify, n_bits  # noqa: E402

COV_OK = {"f300m_e1_n_exp": 3, "f300m_e2_n_exp": 4, "f410m_e1_n_exp": 3, "f410m_e2_n_exp": 4}


def _star(sep, pa):
    return {"groups": {"0": {"sep_arcsec": sep, "pa_deg": pa}}}


def test_n_bits():
    assert n_bits(0) == 0
    assert n_bits(0b1011) == 3


def test_star_wings_first():
    assert classify({"0": COV_OK}, _star(1.5, 120.1)) == {"0": "star_wings"}


def test_fixed_offset_needs_separation_and_angle():
    assert classify({"0": COV_OK}, _star(8.30, 120.5)) == {"0": "fixed_offset_star_image"}
    assert classify({"0": COV_OK}, _star(8.30, 122.0)) == {"0": "unexplained"}
    assert classify({"0": COV_OK}, _star(8.60, 120.1)) == {"0": "unexplained"}


def test_single_exposure_any_band_or_epoch():
    cov = dict(COV_OK, f410m_e2_n_exp=1)
    assert classify({"0": cov}, _star(30.0, 10.0)) == {"0": "single_exposure"}
    assert classify({"0": COV_OK}, _star(30.0, 10.0)) == {"0": "unexplained"}
