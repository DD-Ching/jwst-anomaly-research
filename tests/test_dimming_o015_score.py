import importlib.util
import sys
from pathlib import Path


def _load():
    p = Path("results/dimming/o015_k1/score.py")
    sys.path.insert(0, "results/dimming/nan_k1")
    spec = importlib.util.spec_from_file_location("o015_score", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _c(k, n, kp=0, np_=100):
    m = _load()
    return {
        "le4": {"on_sky": k, "n": n, "wilson": m.wilson(k, n)},
        "persistence": {"wilson": m.wilson(kp, np_)},
    }


def test_decide_rules():
    m = _load()
    assert m.decide(True, _c(5, 9)) == "R3"
    assert m.decide(True, _c(14, 20)) == "R1"
    assert m.decide(True, _c(4, 20)) == "R2"  # interval below 0.584
    assert m.decide(False, _c(14, 20)) == "R2"  # G0 fails
    assert m.decide(True, _c(14, 20, kp=10)) == "R2"  # persistence too common
