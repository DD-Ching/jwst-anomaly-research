"""Frozen stage 3c chance-match control of scripts/m1931_stage3c.py (offline)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import m1931_stage3c as s3c  # noqa: E402


def _grid():
    ra = 292.95 + np.arange(5) * 20 / 3600.0
    dec = np.full(5, -26.57)
    return ra, dec


def test_shifted_counts_zero_for_exact_catalogue():
    ra, dec = _grid()
    assert s3c.shifted_counts(ra, dec, ra, dec) == [0] * s3c.CHANCE_N_DIRECTIONS


def test_shifted_counts_finds_rows_at_the_shift():
    ra, dec = _grid()
    north = dec + s3c.CHANCE_SHIFT_ARCSEC / 3600.0  # a row 5" north of every image
    counts = s3c.shifted_counts(ra, dec, ra, north)
    assert counts[0] == 5 and sum(counts[1:]) == 0


def test_chance_passes_threshold():
    assert s3c.chance_passes(20, [10] * 8)
    assert not s3c.chance_passes(19, [10] * 8)
    assert s3c.chance_passes(0, [0] * 8)


def test_add_chance_needs_both_rules():
    ra, dec = _grid()
    images = {"ra": ra, "dec": dec}
    cat = {"ra": ra, "dec": dec}
    f = {"passes": True, "n_matched": 5}
    s3c.add_chance(f, images, cat)
    assert f["passes"] and f["chance_passes"] and f["passes_stage3_rule"]
    f = {"passes": False, "n_matched": 5}
    s3c.add_chance(f, images, cat)
    assert not f["passes"] and f["chance_passes"]
