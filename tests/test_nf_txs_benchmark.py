"""Offline tests for the TXS 0506+056 SkyLLH benchmark helpers (Neutrino Frontier benchmark 1)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import nf_txs_benchmark as tb  # noqa: E402


def test_rounded_window_is_centred():
    assert tb.rounded_window(57020.0, 185.0) == (56927.5, 57112.5)


def test_angsep():
    assert tb.angsep_deg(10.0, 0.0, 11.0, 0.0) == pytest.approx(1.0)
    assert tb.angsep_deg(0.0, 89.0, 180.0, 89.0) == pytest.approx(2.0)


def test_snap_edges_moves_to_nearest_on_source_event_only():
    mjd = np.array([56927.0, 56927.86, 56928.0, 57112.653, 57112.4, 57000.0])
    sep = np.array([0.5, 0.77, 2.0, 0.38, 1.5, 0.1])
    # 56927.86 is nearer than 56927.0; 56928.0 and 57112.4 are off source
    assert tb.snap_edges(mjd, sep, (56927.5, 57112.5)) == (56927.86, 57112.653)


def test_snap_edges_keeps_edge_without_event_in_reach():
    mjd = np.array([56925.0, 57115.0])
    sep = np.array([0.1, 0.1])
    assert tb.snap_edges(mjd, sep, (56927.5, 57112.5)) == (56927.5, 57112.5)


def test_snap_edges_rejects_empty_window():
    mjd = np.array([100.4])
    sep = np.array([0.1])
    with pytest.raises(ValueError):
        tb.snap_edges(mjd, sep, (100.0, 100.5))


def test_agrees_tolerance():
    assert tb.agrees({"ns": 12.72, "gamma": 2.26}, tb.PAPER["skyllh"])
    assert not tb.agrees({"ns": 11.1, "gamma": 2.22}, tb.PAPER["skyllh"])


def test_manifest_covers_both_irf_versions():
    man = Table.read(tb.MANIFEST, format="ascii.ecsv")
    for v in tb.SAMPLES:
        files = set(man["file"][man["irf_version"] == v])
        assert files == {
            "events/IC86_IV_exp.csv",
            "uptime/IC86_IV_exp.csv",
            "irfs/IC86_effectiveArea.csv",
            "irfs/IC86_smearing.csv",
        }
    assert all(len(m) == 32 for m in man["md5"])
