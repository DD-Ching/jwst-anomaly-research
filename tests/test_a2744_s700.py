"""Offline tests for scripts/a2744_s700.py (pure helpers; no network)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import a2744_s700 as s7  # noqa: E402

RA0, DEC0 = 3.5797, -30.3577
ARC = 1 / 3600.0


def _t(offsets):
    cd = np.cos(np.radians(DEC0))
    return Table(
        {
            "ra": [RA0 + dx * ARC / cd for dx, _ in offsets],
            "dec": [DEC0 + dy * ARC for _, dy in offsets],
            "magnification": [-3.0] * len(offsets),
        }
    )


def test_assign_is_one_to_one():
    cat = _t([(0, 0), (1.7, 0), (15, 5)])
    # one predicted image between a and b: nearest is shared, the assignment is not
    pred = _t([(0.85, 0), (15.2, 5), (40, 0)])
    idx, sep = s7.assign(pred, cat)
    assert len(set(idx)) == 3
    assert sep[0] < 1.0 and sep[2] < 1.0 and sep[1] > 20
    near = s7.nearest(pred, cat)
    assert np.all(near < 1.0)


def test_assign_too_few_images():
    idx, sep = s7.assign(_t([(0, 0)]), _t([(0, 0), (2, 0)]))
    assert idx == [] and np.all(np.isnan(sep))


def test_verdict_outcomes():
    rows = [
        {"z": 1.0, "rms_arcsec": 3.0, "reproduced": False},
        {"z": 1.1, "rms_arcsec": None, "reproduced": False},
    ]
    v = s7.verdict(rows)
    assert v["outcome"].startswith("B") and v["best_assignment"]["z"] == 1.0
    rows.append({"z": 1.2, "rms_arcsec": 0.4, "reproduced": True})
    v = s7.verdict(rows)
    assert v["outcome"].startswith("A") and v["z_reproduced"] == [1.2, 1.2]
