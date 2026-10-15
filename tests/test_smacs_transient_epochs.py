"""Tests for scripts/smacs_transient_epochs.py (offline: the aperture photometry only)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location(
    "smacs_transient_epochs", _DIR / "smacs_transient_epochs.py"
)
ste = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ste)


def test_aperture_flux_recovers_point_source_on_flat_background() -> None:
    pix = 0.03
    rng = np.random.default_rng(1)
    sci = 0.5 + rng.normal(0.0, 0.01, (81, 81))
    sci[40, 40] += 10.0  # all source flux in one pixel at the aperture centre
    f, e = ste.aperture_flux(sci, 40.0, 40.0, pix)
    k = pix**2 * ste.MJYSR_ARCSEC2_TO_UJY
    assert abs(f - 10.0 * k) < 5 * e
    assert 0 < e < 0.1 * f


def test_aperture_flux_is_zero_without_a_source() -> None:
    rng = np.random.default_rng(2)
    sci = 3.0 + rng.normal(0.0, 0.01, (81, 81))
    f, e = ste.aperture_flux(sci, 40.0, 40.0, 0.03)
    assert abs(f) < 4 * e
