"""Offline tests for scripts/w3_short_te.py (frozen dev / validation split)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("w3_short_te", _DIR / "w3_short_te.py")
st = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(st)


def test_split_depends_on_the_injection_not_on_row_order():
    rows = [
        {
            "event_id": f"gb12-R-1-1-{i}",
            "tE": 3.0,
            "rho": 0.01,
            "u0": i / 997,
            "t0": 5000.0 + i,
            "Is": 19.0,
        }  # fmt: skip
        for i in range(2000)
    ]
    a = [st.split_of(r, 3003) for r in rows]
    b = [st.split_of(r, 3003) for r in reversed(rows)][::-1]
    assert a == b  # order-free and deterministic
    assert abs(a.count("dev") - 1000) < 100
    assert a != [st.split_of(r, 3004) for r in rows]  # the seed enters


def test_tracked_sample_uses_the_job_key():
    from astropy.table import Table

    p = st.OUT / "injections_gb12_te3_seed3003.ecsv.gz"
    tab = Table.read(p, format="ascii.ecsv")
    assert list(tab["split"]) == [st.split_of(r, 3003) for r in tab]
