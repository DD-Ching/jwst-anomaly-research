import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import smacs_zphot_stratify as s  # noqa: E402


def test_stratum_ranges_follow_plan():
    z = np.array([0.25, 0.39, 0.55, 0.6, 0.8, 3.0, -1.0, np.nan, 0.1])
    assert list(s.stratum_of(z)) == [
        "member",
        "member",
        "member",
        "other",
        "background",
        "background",
        "other",
        "other",
        "other",
    ]


def test_decision_needs_both_enrichment_and_count():
    assert s.decide(4, 0.2, 0.1).startswith("population effect")
    assert s.decide(3, 0.15, 0.05) == "no rule change"  # fewer than 4 members
    assert s.decide(5, 0.25, 0.21) == "no rule change"  # not 2x the stratum share
    assert s.decide(2, 0.1, 0.209) == "no rule change"  # the observed case
