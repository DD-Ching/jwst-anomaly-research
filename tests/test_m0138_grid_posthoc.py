"""Tests for scripts/m0138_grid_posthoc.py."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import a1689_grid_posthoc as ap  # noqa: E402
import m0138_grid_posthoc as mp  # noqa: E402


def test_configure_points_reused_functions_at_this_pair():
    saved = (ap.E1, ap.E2, ap.BANDS, ap.STAR_BAND)
    try:
        mp.configure()
        assert (ap.E1, ap.E2) == (mp.E1, mp.E2)
        assert ap.BANDS == ("F150W", "F444W") and ap.STAR_BAND == "F444W"
        assert ap.STAR_BAND in ap.BANDS
    finally:
        ap.E1, ap.E2, ap.BANDS, ap.STAR_BAND = saved
