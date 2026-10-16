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


def test_across_line_samples_both_sides_of_the_line():
    S = np.zeros((21, 21))
    S[9, 13], S[11, 13] = -3.0, 4.0  # horizontal spike y = 10 through source (5, 10)
    v = esc.across_line(S, 5.0, 10.0, 0.0, 13.0, 12.0, off_px=1.0)  # entry 2 px off the line
    assert v == [-3.0, 4.0]
    assert esc.across_line(S, 5.0, 10.0, 0.0, 30.0, 10.0) is None  # off the map


def test_symmetric_spikes_drop_unpaired_rays():
    rays = {"A": [0.5, 30.5, 89.5, 210.5], "B": [180.5, 121.5, 270.5]}
    keep = esc.symmetric_spikes(rays)
    assert 121.5 not in keep and 30.5 in keep and 270.5 in keep and 0.5 in keep


def test_onspike_class_and_grid():
    assert esc.onspike_class(-3.0, [-2.5, -2.1, -3.0, -2.0, 1.0, -1.9, -0.5, 2.5]) == (
        4,
        "explained by the spike shift",
    )
    assert esc.onspike_class(-3.0, [-2.5, 2.5, 3.0, -1.0] + [0.0] * 4)[1] == "undetermined"
    pts = esc.onspike_grid(15.72722, -49.25394, 15.72734, -49.25296)
    assert len(pts) == 8
