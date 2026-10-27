import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from m2129_grid_posthoc import classify  # noqa: E402

COV_OK = {"f277w_e1_n_exp": 3, "f277w_e2_n_exp": 4, "f444w_e1_n_exp": 3, "f444w_e2_n_exp": 4}


def _star(sep, pa):
    return {"groups": {"0": {"sep_arcsec": sep, "pa_deg": pa}}}


def test_star_wings_and_single_exposure():
    assert classify({"0": COV_OK}, _star(2.7, 17.0), 77.16) == {"0": ["star_wings"]}
    cov = dict(COV_OK, f444w_e1_n_exp=1)
    assert classify({"0": cov}, _star(30.0, 10.0), 77.16) == {"0": ["single_exposure"]}
    assert classify({"0": COV_OK}, {"groups": {}}, 77.16) == {"0": ["unexplained"]}


def test_fixed_offset_vector_rotates_with_roll():
    # 120.1 deg at PA_APER 75.12 becomes 122.14 deg at PA_APER 77.16
    assert classify({"0": COV_OK}, _star(8.33, 122.1), 77.16) == {"0": ["fixed_offset_star_image"]}
    assert classify({"0": COV_OK}, _star(8.33, 120.1), 77.16) == {"0": ["unexplained"]}
    # wraps through 360 deg
    assert classify({"0": COV_OK}, _star(8.33, 0.2), 315.0) == {"0": ["fixed_offset_star_image"]}


def test_off_footprint():
    cov = dict(COV_OK, f277w_e1_n_exp=0)
    assert classify({"0": cov}, _star(30.0, 10.0), 77.16) == {"0": ["off_footprint"]}
