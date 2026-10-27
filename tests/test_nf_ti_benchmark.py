"""Offline tests for the time-integrated SkyLLH benchmark (Neutrino Frontier benchmark 2)."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import nf_ti_benchmark as ti  # noqa: E402


def test_agrees_uses_relative_ns_tolerance_for_large_ns():
    target = {"ns": 80.1, "gamma": 3.2}
    assert ti.agrees({"ns": 84.0, "gamma": 3.25}, target)  # |dn_s| 3.9 <= 4.005
    assert not ti.agrees({"ns": 84.2, "gamma": 3.2}, target)
    assert not ti.agrees({"ns": 80.1, "gamma": 3.31}, target)


def test_agrees_floor_of_one_event_for_small_ns():
    target = {"ns": 8.8, "gamma": 2.0}
    assert ti.agrees({"ns": 9.7, "gamma": 2.0}, target)
    assert not ti.agrees({"ns": 9.9, "gamma": 2.0}, target)


def test_targets_are_the_papers_four_sources():
    assert set(ti.SOURCES) == {"NGC 1068", "PKS 1424+240", "TXS 0506+056", "GB6 J1542+6129"}
    assert set(ti.DV_VERSION) == set(ti.SAMPLES)
