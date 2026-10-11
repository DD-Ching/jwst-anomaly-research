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


def test_split_is_frozen_and_balanced():
    a = [st.split_of(i, 3003) for i in range(2000)]
    assert a == [st.split_of(i, 3003) for i in range(2000)]  # deterministic
    assert set(a) == {"dev", "validation"}
    assert abs(a.count("dev") - 1000) < 100
    # the tracked sample's split (998 dev / 1002 validation) is the one this function gives
    assert a.count("dev") == 998
