"""Tests for scripts/w1_row_veto_eval.py (frozen adoption rule)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w1_row_veto_eval as ev  # noqa: E402


def test_adopt_needs_k1r_and_every_field_cheap():
    ok = {f: True for f in ev.NCF}
    assert ev.adopt(True, ok)
    assert not ev.adopt(False, ok)
    assert not ev.adopt(True, {**ok, "abell370-ncf": False})
    assert not ev.adopt(True, {f: True for f in ev.NCF[:4]})  # a missing field is not a pass
