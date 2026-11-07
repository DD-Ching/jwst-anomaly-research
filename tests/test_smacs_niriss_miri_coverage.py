"""Footprint and redshift-catalogue matching of scripts/smacs_niriss_miri_coverage.py."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pytest
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import smacs_niriss_miri_coverage as m  # noqa: E402

CAT = """#Spectroscopic redshift catalogue (synthetic)
ID RA DEC Z_GRISM ZERR_GRISM Q_FLAG Z_SPEC REF_SPEC
1 10.0 -70.0 1.5 0.001 2 "" ""
2 10.01 -70.0 "" "" -99 0.4 nirspec
"""


def test_footprints_keep_level3_only():
    t = Table(
        {
            "obs_id": ["a", "a", "b"],
            "calib_level": [3, 2, 3],
            "s_region": ["POLYGON 0 0 1 0 1 1", "POLYGON 5 5 6 5 6 6", "POLYGON 0 0 2 0 2 2"],
        }
    )
    fp = m.footprints(t)
    assert sorted(fp) == ["a", "b"] and len(fp["a"]) == 1


def test_covered_by():
    cov = {"x": {"inside": True}, "y": {"inside": False}}
    assert m.covered_by(cov) == ["x"]


def test_noirot_matches(tmp_path, monkeypatch):
    p = tmp_path / "cat.txt"
    p.write_text(CAT)
    monkeypatch.setattr(m, "NOIROT_SHA256", hashlib.sha256(p.read_bytes()).hexdigest())
    recs = [{"ra": 10.0, "dec": -70.0 + 0.1 / 3600}, {"ra": 10.0, "dec": -69.99}]
    out = m.noirot_matches(p, recs)
    assert out[0]["matched"] and out[0]["nearest_id"] == 1 and out[0]["q_flag"] == 2
    assert out[0]["z_spec"] is None
    assert not out[1]["matched"]


def test_noirot_sha_mismatch(tmp_path):
    p = tmp_path / "cat.txt"
    p.write_text(CAT)
    with pytest.raises(ValueError, match="sha256"):
        m.noirot_matches(p, [{"ra": 10.0, "dec": -70.0}])
