"""Tests for scripts/elgordo_spike_check.py (synthetic; no network)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location(
    "elgordo_spike_check", _DIR / "elgordo_spike_check.py"
)
esc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(esc)


def test_classify_on_spike_rule():
    assert esc.classify(179.0, 3.5, [0.5, 180.5]) == (True, 1.5)
    assert esc.classify(175.0, 3.5, [0.5, 180.5])[0] is False  # 5.5° away
    assert esc.classify(180.0, 7.0, [180.5])[0] is False  # too far
    assert esc.classify(359.5, 2.0, [0.5])[0] is True  # wraps around 0°


def test_ray_angles_finds_spikes():
    yy, xx = np.indices((101, 101), dtype=float)
    img = np.zeros((101, 101))
    for ang in (0, 90, 180, 270):
        t = np.radians(ang)
        for r in range(5, 50):
            img[int(round(50 + r * np.sin(t))), int(round(50 + r * np.cos(t)))] = 10.0
    a = esc.ray_angles(img, 50.0, 50.0, (10, 40), n=4)
    for want in (0.5, 90.5, 180.5, 270.5):
        assert min(min(abs(x - want), 360 - abs(x - want)) for x in a) <= 1.0


def test_across_samples_both_sides():
    S = np.zeros((21, 21))
    S[9, 10], S[11, 10] = -3.0, 4.0  # horizontal spike through (10, 10): sides are ±y
    v = esc.across(S, 10.0, 10.0, 0.0, off_px=1.0)
    assert v[0] * v[1] < 0
