import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import arcs_pa_convention as apc  # noqa: E402


def test_medians_identify_each_convention():
    rng = np.random.default_rng(1)
    pa = rng.uniform(0, 180, 50)
    noisy = pa + rng.normal(0, 5, 50)
    assert apc.outcome(50, apc.medians(noisy, pa)) == "convention agrees"
    assert apc.outcome(50, apc.medians(180 - noisy, pa)) == "convention differs: mirrored"
    assert apc.outcome(50, apc.medians(noisy + 90, pa)) == "convention differs: rotated90"


def test_outcome_undetermined_cases():
    good = {"identity": 5.0, "mirrored": 45.0, "rotated90": 80.0}
    assert apc.outcome(7, good) == "undetermined"
    assert apc.outcome(8, good) == "convention agrees"
    assert (
        apc.outcome(20, {"identity": 25.0, "mirrored": 45.0, "rotated90": 60.0}) == "undetermined"
    )
    assert (
        apc.outcome(20, {"identity": 15.0, "mirrored": 25.0, "rotated90": 60.0}) == "undetermined"
    )
    rng = np.random.default_rng(2)
    assert apc.outcome(60, apc.medians(rng.uniform(0, 180, 60), rng.uniform(0, 180, 60))) == (
        "undetermined"
    )
