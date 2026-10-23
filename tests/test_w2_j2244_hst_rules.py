import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w2_j2244_hst_rules as r  # noqa: E402


def test_threshold_from_development_negatives():
    assert r.threshold([3.0, 10.0, 7.5]) == 25.0
    assert r.threshold([3.0, 40.0]) == 80.0
    assert r.threshold([]) == 25.0
    assert r.threshold([None, float("nan"), 12.0]) == 25.0


def test_status_order():
    assert r.status(None, 25, 0.1, 1.0, True) == "fit unreliable"
    assert r.status(100, 25, 0.1, 1.0, False) == "fit unreliable"
    assert r.status(100, 25, 0.1, 10.5, True) == "fit unreliable"
    assert (
        r.status(100, 25, 0.1, 1.0, True, sersic_r=0.05, sersic_flux=10.0, band="F814W")
        == "unresolved"
    )
    assert r.status(100, 25, 0.1, 1.0, True, sersic_r=0.3, sersic_flux=0.0, band="F814W") == "none"
    assert r.status(20, 25, 0.1, 1.0, True) == "none"
    assert r.status(30, 25, 0.7, 1.0, True) == "extended light offset"
    assert r.status(30, 25, 0.3, 1.0, True) == "extended light"


def test_gate():
    neg = ["none"] * 5 + ["fit unreliable"]
    pos = ["extended light"] * 3 + ["none"]
    g = r.gate(neg, pos, True)
    assert g["n_neg"] == 5 and g["n_pos"] == 4 and g["pos_detected"] == 3
    assert g["passed"] is False  # 3 / 4 < 0.8
    assert r.gate(neg, ["extended light"] * 4 + ["none"], True)["passed"] is True  # 4 / 5 = 0.8
    assert r.gate(neg[:-2] + ["extended light"], ["extended light"] * 5, True)["passed"] is False
    assert r.gate(neg, ["extended light"] * 5, False)["passed"] is False
    assert r.gate(["none"] * 4, ["extended light"] * 5, True)["enough_controls"] is False


def test_contiguous_limit():
    assert r.contiguous_limit([19, 20, 21, 22], [True, True, False, True]) == 20
    assert r.contiguous_limit([19, 20], [False, True]) is None


def test_geometry_rule():
    assert r.geometry_ok(1.02) and r.geometry_ok(0.7) and r.geometry_ok(1.8)
    assert not r.geometry_ok(2.27) and not r.geometry_ok(0.6)
    assert not r.geometry_ok(None) and not r.geometry_ok(float("nan"))
