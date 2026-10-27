"""Tests for scripts/grid_g1_coverage_posthoc.py."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import grid_g1_coverage_posthoc as gp  # noqa: E402


def test_split_recovery_counts_each_cell():
    inj = Table(
        {
            "frac": [1.0, 1.0, 1.0, 1.0, 0.5],
            "nominal": [True, True, False, False, True],
            "e1_n_exp_min": [1, 3, 2, 1, 1],
            "recovered": [True, False, True, False, True],
        }
    )
    s = gp.split_recovery(inj, 1.0)
    assert s["nominal_e1n_1"] == {"n": 1, "recovered_frac": 1.0}
    assert s["nominal_e1n_>=2"] == {"n": 1, "recovered_frac": 0.0}
    assert s["floor_e1n_>=2"] == {"n": 1, "recovered_frac": 1.0}
    assert s["floor_e1n_1"] == {"n": 1, "recovered_frac": 0.0}


def test_floor_significance_at_equal_errors():
    e = np.full(10, 0.02)
    assert gp.floor_significance(e, e, 1.0, 0.2) == pytest.approx(10 / np.sqrt(2))
    assert gp.floor_significance(e, e, 2.0, 0.2) == pytest.approx(5 / np.sqrt(2))
