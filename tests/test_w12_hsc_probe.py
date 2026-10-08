"""Tests for scripts/w12_hsc_probe.py ``classify`` (tiny synthetic tables, offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("w12_hsc_probe", _DIR / "w12_hsc_probe.py")
probe = importlib.util.module_from_spec(_spec)
sys.modules["w12_hsc_probe"] = probe  # dataclasses look their module up
_spec.loader.exec_module(probe)

D = 1 / 3600
P = probe.Params()


def _src(offsets, ci):
    """Sources at (dx, dy) arcsec offsets from (10, 0) with concentration indices ci."""
    return Table(
        {
            "MatchRA": [10 + dx * D for dx, _ in offsets],
            "MatchDec": [dy * D for _, dy in offsets],
            "CI": ci,
        }
    )


def test_pair_with_lens_between_is_deflector():
    s = _src([(-0.6, 0), (0.6, 0), (0.05, 0)], [1.0, 1.05, 2.0])
    status, n, sep = probe.classify(10.0, 0.0, s, P)
    assert (status, n) == ("deflector", 2)
    assert abs(sep - 1.2) < 1e-3


def test_pair_without_extended_source_is_none():
    s = _src([(-0.6, 0), (0.6, 0), (2.5, 2.5)], [1.0, 1.05, 2.0])  # extended source far away
    assert probe.classify(10.0, 0.0, s, P)[0] == "none"


def test_extended_source_on_an_image_is_not_the_lens():
    s = _src([(-0.6, 0), (0.6, 0), (0.65, 0)], [1.0, 1.05, 2.0])
    assert probe.classify(10.0, 0.0, s, P)[0] == "none"


def test_fewer_than_two_point_images_is_undecided():
    s = _src([(0, 0), (0.3, 0)], [1.0, 2.0])
    assert probe.classify(10.0, 0.0, s, P)[:2] == ("undecided", 1)
    assert probe.classify(10.0, 0.0, Table(), P)[0] == "undecided"
