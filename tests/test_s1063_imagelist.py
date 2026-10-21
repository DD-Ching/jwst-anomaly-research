"""Offline tests for scripts/s1063_imagelist.py (no network)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import MaskedColumn, Table

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import s1063_imagelist as si  # noqa: E402

OUT = ROOT / "results" / "lens_consistency" / "abells1063_images"


def _rows():
    return Table(
        {
            "ID": ["1.1", "1.2", "2.1", "2.2", "3a.1"],
            "RAJ2000": ["22 48 46.668"] * 5,
            "DEJ2000": ["-44 31 37.21", "-44 31 38", "-44 31 39", "-44 31 40", "-44 31 37.21"],
            "zspec": MaskedColumn(
                [1.229, 0.0, 0.0, 0.0, 3.0], mask=[False, True, True, True, False]
            ),
            "Geoz": [1.45, 1.45, 2.5, 2.4, 3.1],
            "CL": [3, 3, 1, 1, 2],
        }
    )


def test_system_of():
    assert si.system_of("12.3") == "12"
    assert si.system_of("1a.1") == "1a"


def test_build_propagates_system_zspec_and_falls_back_to_geoz():
    out, summary = si.build(_rows(), "2026-10-21")
    assert list(out["z_spec"]) == [True, True, False, False, True]
    np.testing.assert_allclose(out["z"], [1.229, 1.229, 2.5, 2.4, 3.0])
    np.testing.assert_allclose(out["ra"][0], 342.19445, atol=1e-6)
    assert summary["n_systems"] == 3
    assert summary["n_systems_zspec"] == 2
    assert summary["n_images_system_zspec"] == 3
    assert summary["cl_counts"] == {"1": 2, "2": 1, "3": 2}
    assert out.meta["provenance"] == "observed"
    assert list(out["dup_position"]) == [True, False, False, False, True]
    assert summary["dup_position_ids"] == ["1.1", "3a.1"]


def test_committed_list_matches_summary():
    tab = Table.read(OUT / "chan20_images.ecsv")
    summary = json.loads((OUT / "summary.json").read_text())
    assert len(tab) == summary["n_images"]
    assert int(np.sum(tab["z_spec"])) == summary["n_images_system_zspec"]
    assert len({si.system_of(i) for i in tab["image_id"]}) == summary["n_systems"]
    assert list(tab["image_id"][tab["dup_position"]]) == summary["dup_position_ids"]
