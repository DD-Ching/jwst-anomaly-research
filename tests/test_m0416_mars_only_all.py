import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from m0416_mars_only_all import gate_r, sheet_rows, split_families  # noqa: E402


def test_split_families_partitions_and_drops_nonfinite():
    t = Table(
        {
            "system": ["A", "A", "B", "C", "C", "D"],
            "beta_x": [0.0, 0.1, 1.0, np.nan, 2.0, 3.0],
            "beta_y": [0.0, 0.0, 1.0, 2.0, 2.0, 3.0],
        }
    )
    pred = Table({"system": ["A", "C"]})
    assert split_families(t, pred) == (["A"], ["B", "D"])


def test_gate_r():
    known = ("D115a", 64.044742, -24.060901)
    hit = {
        "system": "D115a",
        "ra": 64.044742 + 0.05 / 3600 / np.cos(np.radians(-24.06)),
        "dec": -24.060901,
    }
    assert gate_r([hit], known)
    assert not gate_r([{**hit, "dec": hit["dec"] + 0.2 / 3600}], known)
    assert not gate_r([{**hit, "system": "D115b"}], known)
    assert not gate_r([hit, hit], known)
    assert not gate_r([], known)


def test_sheet_rows_orders_and_caps():
    rows = [
        {"reading": "no_reading", "pred_over_sigma": 99.0},
        {"reading": "recovered_x", "pred_over_sigma": 2.0},
        {"reading": "lensing_residual_input", "pred_over_sigma": 5.0},
        {"reading": "lensing_residual_input", "pred_over_sigma": float("nan")},
    ]
    out = sheet_rows(rows, n=2)
    assert [r["pred_over_sigma"] for r in out] == [5.0, 2.0]
