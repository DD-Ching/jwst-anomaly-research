import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w2_j2244_stpsf_rules as r  # noqa: E402


def test_hookup():
    assert r.hookup(0.05) == "pass"
    assert r.hookup(0.2) == "pass (flagged)"
    assert r.hookup(0.31) == "stop"
    assert r.hookup(math.nan) == "stop"


def test_star_ok_needs_every_exposure_and_annulus():
    assert r.star_ok([1.0, 1.9], [[0.5, -2.9], [1.0, 2.0]])
    assert not r.star_ok([1.0, 2.1], [[0.5], [1.0]])
    assert not r.star_ok([1.0, 1.0], [[0.5], [-3.1]])


def test_pair_and_injection():
    assert r.pair_ok(2.9) and not r.pair_ok(-3.0) and not r.pair_ok(math.nan)
    assert r.injection_ok(5.0, 0.5) and not r.injection_ok(4.9, 0.9)
    assert not r.injection_ok(9.0, 0.49)


def test_worth_validation_and_combine():
    assert r.worth_validation("pass (flagged)", True, True, True, 5.0)
    assert not r.worth_validation("stop", True, True, True, 9.0)
    assert not r.worth_validation("pass", True, True, True, 4.0)
    assert math.isclose(r.combine([1.0, 1.0], [1.0, 1.0]), math.sqrt(2.0))
