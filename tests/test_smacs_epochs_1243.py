import importlib.util
import sys
from pathlib import Path

_path = Path(__file__).resolve().parents[1] / "scripts" / "smacs_epochs_1243.py"
sys.path.insert(0, str(_path.parent))
_spec = importlib.util.spec_from_file_location("smacs_epochs_1243", _path)
mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mod)


def _st(dm=0.0, sig=0.05, r=5.0, p95=10.0, p99=15.0, n=20):
    return {
        "n_peers": n,
        "target_dmag": dm,
        "peer_dmag_sigma": sig,
        "target_residual_mas": r,
        "peer_residual_p95_mas": p95,
        "peer_residual_p99_mas": p99,
    }


def test_not_measured():
    assert mod.decide(None) == {"variability": "not measured", "motion": "not measured"}


def test_too_few_peers():
    assert mod.decide(_st(n=9)) == {"variability": "undetermined", "motion": "undetermined"}


def test_constant_no_motion():
    assert mod.decide(_st(dm=0.1, r=10.0)) == {
        "variability": "constant",
        "motion": "no motion detected",
    }


def test_variable_moves():
    assert mod.decide(_st(dm=-0.25, r=15.1)) == {"variability": "variable", "motion": "moves"}


def test_between_thresholds_undetermined():
    assert mod.decide(_st(dm=0.15, r=12.0)) == {
        "variability": "undetermined",
        "motion": "undetermined",
    }


def test_exact_thresholds():
    # dm = 2 / 5 sigma -> constant / variable; r = p95 / p99 -> no motion / undetermined
    assert mod.decide(_st(dm=0.1, r=10.0))["variability"] == "constant"
    assert mod.decide(_st(dm=0.25))["variability"] == "variable"
    assert mod.decide(_st(r=10.0))["motion"] == "no motion detected"
    assert mod.decide(_st(r=15.0))["motion"] == "undetermined"
    assert mod.decide(_st(n=10))["variability"] == "constant"


def test_missing_or_zero_sigma_and_missing_dm():
    assert mod.decide(_st(sig=0.0, dm=0.0))["variability"] == "undetermined"
    assert mod.decide(_st(sig=None))["variability"] == "undetermined"
    assert mod.decide(_st(dm=None))["variability"] == "not measured"
    assert mod.decide(_st(p95=None, p99=None))["motion"] == "undetermined"
