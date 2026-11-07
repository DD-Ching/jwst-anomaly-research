"""Frozen outcome rules of scripts/smacs_epoch3_384.py."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import smacs_epoch3_384 as m  # noqa: E402

PREV = {"f150w": [-23.0, 9.0], "f444w": [-22.0, 4.0]}


def _band(res, med=5.0, p95=20.0):
    return {
        "target_residual_mas": list(res),
        "target_residual_total_mas": float(np.hypot(*res)),
        "peer_residual_median_mas": med,
        "peer_residual_p95_mas": p95,
    }


def test_reproduced_needs_both_bands_above_p95_and_aligned():
    both = {"f150w": _band((-24, 10)), "f444w": _band((-21, 6))}
    assert m.decide(both, PREV) == "shift reproduced"
    weak = {"f150w": _band((-24, 10)), "f444w": _band((-15, 5))}
    assert m.decide(weak, PREV) == "undetermined"
    turned = {"f150w": _band((-24, 10)), "f444w": _band((21, 6))}
    assert m.decide(turned, PREV) == "undetermined"


def test_not_reproduced_and_missing_band():
    none = {"f150w": _band((1, 2)), "f444w": _band((-3, 1))}
    assert m.decide(none, PREV) == "shift not reproduced"
    assert m.decide({"f150w": _band((1, 2)), "f444w": _band((-6, 1))}, PREV) == "undetermined"
    assert m.decide({"f150w": _band((1, 2)), "f444w": None}, PREV) == "undetermined"
