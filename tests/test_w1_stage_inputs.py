"""Tests for scripts/w1_stage_inputs.py (offline; synthetic tables only)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import astrodeep_zout  # noqa: E402
import inject_shear as ish  # noqa: E402
import w1_stage_inputs as wsi  # noqa: E402


def _zout(n: int, ra0: float) -> Table:
    t = Table({c: np.arange(n, dtype=float) for c in astrodeep_zout.ZOUT_COLUMNS})
    t["id"] = np.arange(1, n + 1)
    t["ra"] = ra0 + t["ra"] * 1e-3
    t["extra"] = np.zeros(n)  # other zout columns are dropped
    return t


def test_stack_zouts_keeps_order_and_only_the_match_columns():
    a, b = _zout(3, 34.0), _zout(2, 35.0)
    s = wsi.stack_zouts([a, b])
    assert s.colnames == list(astrodeep_zout.ZOUT_COLUMNS)
    assert len(s) == 5 and list(s["ra"][:3]) == list(a["ra"]) and list(s["ra"][3:]) == list(b["ra"])
    sha = astrodeep_zout.content_sha256(s)
    assert sha == astrodeep_zout.content_sha256(wsi.stack_zouts([a, b]))
    assert sha != astrodeep_zout.content_sha256(wsi.stack_zouts([b, a]))


def test_primer_photoz_sources_are_pinned_stacks_and_each_file_is_staged_once():
    for survey in ("uds", "cosmos"):
        src = wsi.PHOTOZ[Path(ish.PRIMER_ZOUT.format(survey)).name]
        assert src["sha256"] and len(src["stack"]) == 2
        for _, member in src["stack"]:
            assert member["sha256"] and member["member_sha256"]
    jobs = wsi.plan(sorted(ish.PRIMER_FIELDS))
    targets = [j[1] for j in jobs]
    assert len(targets) == len(set(targets)) == 4 + 2  # four mosaic catalogues, two stacks
