import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import a2744_s700_nirspec as m  # noqa: E402


def test_decide_gate_first():
    assert m.decide(False, True, True, ["8204"]) == "not validated"


def test_decide_both_ends():
    assert m.decide(True, True, True, []) == "both ends measured"


def test_decide_new_input():
    assert m.decide(True, True, False, ["8204"]) == "new input, not in DJA"


def test_decide_none():
    assert m.decide(True, True, False, []) == "no new input"
    assert m.decide(True, False, False, []) == "no new input"


def test_ends_disjoint():
    assert not set(m.END_A) & set(m.END_B)


def test_release_mjd_masked_and_empty():
    import numpy as np

    assert m.release_mjd(np.ma.masked) is None
    assert m.release_mjd(None) is None
    assert m.release_mjd("") is None
    assert m.release_mjd(float("nan")) is None
    assert m.release_mjd(61000.5) == 61000.5
