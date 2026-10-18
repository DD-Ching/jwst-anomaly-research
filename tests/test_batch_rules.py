"""The batch pass rule and class rule (scripts/batch_run.py, scripts/batch_review.py), offline."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import batch_pairs as bp  # noqa: E402
import batch_review as br  # noqa: E402
import batch_run as bru  # noqa: E402


def test_box_passes_boundaries():
    assert bru.box_passes(7, 10, 2, 0.70, 2)  # exactly at both limits
    assert not bru.box_passes(6, 10, 0, 0.70, 2)  # vanish fraction below
    assert not bru.box_passes(10, 10, 3, 0.70, 2)  # too many false hits
    assert not bru.box_passes(0, 0, 0, 0.70, 2)  # no vanish sources: undeterminable


def test_classify():
    assert br.classify(True, {"F150W": -30, "F444W": -10}, None) == "artifact (edge)"
    assert br.classify(False, {"F150W": -2.9, "F444W": 2.0}, None) == "no significant change"
    assert br.classify(False, {"F150W": -27, "F444W": -13}, None) == "real change"
    assert br.classify(False, {"F150W": -5, "F444W": 4}, None) == "undetermined"
    assert br.classify(False, {"F150W": -5}, None).startswith("undetermined (photometry")
    assert br.classify(False, {"F150W": -27, "F444W": -13}, "veto: spike").startswith("artifact")
    assert (
        br.classify(False, {"F150W": 1, "F444W": 1}, "looked at: fine") == "no significant change"
    )


def test_wilson():
    assert bru.wilson(0, 0) is None
    lo, hi = bru.wilson(62, 71)
    assert lo < 62 / 71 < hi


def test_touching_boxes_do_not_overlap():
    import math

    # batch-1 A370_1 (stored rounded) and its east grid neighbour: 23.999999999999" apart, which the
    # old strict "< BOX" test called overlapping
    a370 = (39.9754, -1.576)
    east = (39.98206918946676, -1.576)
    c = math.cos(math.radians(east[1]))
    assert abs((east[0] - a370[0]) * c * 3600) < bp.BOX  # the old rule's false overlap
    assert not bp.overlaps(*east, [a370])
    step = bp.BOX / 3600 / c
    assert bp.overlaps(a370[0] + step / 2, a370[1], [a370])  # half a box: real overlap
    assert bp.overlaps(*a370, [a370])


def test_run_timeout_kills_and_flags(tmp_path):
    import sys as _sys

    slow = [_sys.executable, "-c", "import time; time.sleep(30)"]
    r = bru._run(slow, tmp_path / "slow.log", timeout=1.0)
    assert r["rc"] != 0 and r.get("timed_out") is True and r["wall_s"] < 10
    fast = [_sys.executable, "-c", "pass"]
    r = bru._run(fast, tmp_path / "fast.log", timeout=30.0)
    assert r["rc"] == 0 and "timed_out" not in r


def test_too_late():
    assert not bru.too_late(100.0, None)
    assert not bru.too_late(100.0, 200.0)
    assert bru.too_late(300.0, 200.0)


def test_rank_key_prefers_extra_epochs_within_depth_tier():
    a = bp.rank_key(1245.0, 0, 0, 0.0)
    b = bp.rank_key(1300.0, 3, 5, 1.0)  # same factor-2 depth tier, more extra epochs
    c = bp.rank_key(2600.0, 0, 9, 1.4)  # deeper tier wins over extra epochs
    assert sorted([a, b, c]) == [c, b, a]
    assert bp.rank_key(1245.0, 2, 1, 0.0) < bp.rank_key(1245.0, 2, 1, 1.0)


def test_run_without_wait4(tmp_path, monkeypatch):
    import sys as _sys

    monkeypatch.delattr(bru.os, "wait4", raising=False)
    r = bru._run([_sys.executable, "-c", "raise SystemExit(3)"], tmp_path / "x.log", timeout=30.0)
    assert r["rc"] == 3 and r["peak_rss_mb"] is None and "timed_out" not in r
    r = bru._run(
        [_sys.executable, "-c", "import time; time.sleep(30)"], tmp_path / "y.log", timeout=1.0
    )
    assert r.get("timed_out") is True


def test_parse_utc_and_boundary():
    assert bru.parse_utc("2026-10-17T18:20:00") == bru.parse_utc("2026-10-17T20:20:00+02:00")
    stop = bru.parse_utc("2026-10-17T18:20:00")
    assert not bru.too_late(stop, stop) and bru.too_late(stop + 1, stop)


def test_task_problem():
    assert bru.task_problem(None) == "missing"
    assert bru.task_problem({"arc_diff": {"rc": -9, "timed_out": True}}) == "timed out"
    assert (
        bru.task_problem({"arc_diff": {"rc": -1, "why": "not started: stop time"}}) == "not started"
    )
    assert bru.task_problem({"arc_diff": {"rc": 1}}) == "failed"
    assert bru.task_problem({"arc_diff": {"rc": 0}, "morph": {"rc": 0}}) is None


def test_same_source_verdict():
    import same_source_pairs as ssp

    good = {
        "retrieval": {"recall@10": 0.5, "median_false_matches": 10.0},
        "control": {"decoy_match_rate": 0.2, "positive_match_rate": 0.6},
    }
    assert ssp.verdict(good) == {"P1": True, "P2": True, "pass": True}
    bad = {
        "retrieval": {"recall@10": 0.49, "median_false_matches": 10.5},
        "control": {"decoy_match_rate": 0.1, "positive_match_rate": 0.9},
    }
    assert ssp.verdict(bad) == {"P1": False, "P2": False, "pass": False}
    missing = {"retrieval": {"recall@10": None, "median_false_matches": None}, "control": {}}
    assert ssp.verdict(missing)["pass"] is False
    assert ssp.holdout("1") == (__import__("zlib").crc32(b"1") % 2 == 0)
    assert ssp.auc([1, 2], [3, 4]) == 1.0 and ssp.auc([3], [3]) == 0.5
