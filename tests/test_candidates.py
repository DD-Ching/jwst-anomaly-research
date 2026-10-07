import sqlite3

import numpy as np
import pytest

from jwst_anomaly.candidates import (
    SCHEMA_VERSION,
    STATUSES,
    TRANSITIONS,
    CandidateStore,
    StoreError,
    to_python,
)

CONTEXT = {
    "run_id": "r1",
    "created_utc": "2026-10-07T00:00:00Z",
    "config": {"sha256": "ab" * 32},
    "git": {"commit": "c0ffee", "dirty": False},
}


def _record(uid, rank, sample="s1", **extra):
    return {"source_uid": uid, "sample_id": sample, "rank": rank, "score": 10.0 - rank, **extra}


@pytest.fixture
def store(tmp_path):
    with CandidateStore(tmp_path / "sub" / "c.sqlite") as s:
        s.add_run("r1", CONTEXT, config_name="cfg")
        s.add_candidates(
            "r1",
            [
                _record(
                    "a",
                    1,
                    ra=np.float64(110.8),
                    dec=-73.45,
                    top_features="mag_ref",
                    method_scores={"robust_z": np.float32(3.5), "lof": np.nan},
                    flags={"F200W": {"quality_flag": "ok", "on_edge": np.bool_(False)}},
                    xmatch={"n_matches": np.int64(0), "best_match_sep_arcsec": np.nan},
                    cutouts={"F200W": "s1/cutouts/a_F200W.fits"},
                ),
                _record("b", 2),
                _record("a", 1, sample="s2"),  # same uid in another sample of the same run
            ],
        )
        yield s


def _triage(store, uid="a", **kw):
    return store.set_status(uid, "triaged", note="looked at it", author="tester", **kw)


def _vet(store, uid="a", **kw):
    return store.add_vetting_note(
        uid,
        test="edge_check",
        outcome="pass",
        evidence="cutout s1/cutouts/a_F200W.fits; 0% NaN",
        provenance="observed",
        author="tester",
        **kw,
    )


def test_schema_version_and_reopen(tmp_path, store):
    assert store.schema_version == SCHEMA_VERSION
    path = store.path
    store.close()
    with CandidateStore(path, create=False) as again:
        assert again.schema_version == SCHEMA_VERSION
        assert len(again.list_candidates()) == 3


def test_refuses_newer_schema_and_foreign_db(tmp_path):
    newer = tmp_path / "newer.sqlite"
    con = sqlite3.connect(newer)
    con.execute(f"PRAGMA user_version = {SCHEMA_VERSION + 1}")
    con.close()
    with pytest.raises(StoreError, match="schema version"):
        CandidateStore(newer)
    foreign = tmp_path / "foreign.sqlite"
    con = sqlite3.connect(foreign)
    con.execute("CREATE TABLE x (y)")
    con.commit()
    con.close()
    with pytest.raises(StoreError, match="not a candidate store"):
        CandidateStore(foreign)


def test_create_false_requires_existing(tmp_path):
    with pytest.raises(FileNotFoundError):
        CandidateStore(tmp_path / "missing.sqlite", create=False)
    assert not (tmp_path / "missing.sqlite").exists()


def test_runs(store):
    store.finish_run("r1", "completed", summary={"n_candidates": 3})
    run = store.get_run("r1")
    assert run["status"] == "completed" and run["finished_utc"]
    assert run["context"]["git"]["commit"] == "c0ffee"
    assert run["git_dirty"] is False
    assert run["summary"] == {"n_candidates": 3}
    store.add_run("r2", {**CONTEXT, "run_id": "r2"}, config_name="cfg")
    assert store.latest_run_id() == "r2"
    assert store.latest_run_id(status="completed") == "r1"
    assert [r["run_id"] for r in store.list_runs()] == ["r2", "r1"]
    assert store.list_runs()[1]["n_candidates"] == 3
    assert store.get_run("nope") is None
    with pytest.raises(StoreError, match="unknown run"):
        store.finish_run("nope", "failed")
    with pytest.raises(StoreError, match="run status"):
        store.finish_run("r1", "exploded")


def test_candidate_fields_roundtrip(store):
    cand = store.get_candidate("a", sample_id="s1")
    assert cand["ra"] == 110.8 and cand["rank"] == 1 and cand["status"] == "new"
    assert cand["method_scores"] == {"lof": None, "robust_z": 3.5}
    assert cand["flags"] == {"F200W": {"on_edge": False, "quality_flag": "ok"}}
    assert cand["xmatch"] == {"best_match_sep_arcsec": None, "n_matches": 0}
    assert cand["cutouts"] == {"F200W": "s1/cutouts/a_F200W.fits"}
    assert cand["vetting_notes"] == [] and cand["status_history"] == []


def test_duplicate_candidate_rejected(store):
    with pytest.raises(StoreError, match="cannot add"):
        store.add_candidates("r1", [_record("b", 2)])
    with pytest.raises(StoreError, match="source_uid"):
        store.add_candidates("r1", [{"source_uid": " "}])


def test_list_filters(store):
    assert [c["source_uid"] for c in store.list_candidates(sample_id="s1")] == ["a", "b"]
    _triage(store, "b")
    assert [c["source_uid"] for c in store.list_candidates(status="triaged")] == ["b"]
    assert len(store.list_candidates(status="new", limit=1)) == 1
    assert store.list_candidates(run_id="other") == []
    with pytest.raises(StoreError, match="status must be"):
        store.list_candidates(status="bogus")


def test_ambiguous_uid_needs_sample(store):
    with pytest.raises(StoreError, match="ambiguous"):
        store.get_candidate("a")
    assert store.get_candidate("a", sample_id="s2")["sample_id"] == "s2"
    with pytest.raises(StoreError, match="no candidate"):
        store.get_candidate("zzz")
    with pytest.raises(StoreError, match=r"no candidate 'b' \(sample 's2'\) in any run"):
        store.get_candidate("b", sample_id="s2")
    with pytest.raises(StoreError, match=r"no candidate 'b' \(sample 's2'\) in run 'r1'"):
        store.get_candidate("b", run_id="r1", sample_id="s2")


def test_uid_in_several_runs_needs_run_id(store):
    store.add_run("r2", {**CONTEXT, "run_id": "r2", "created_utc": "2026-10-08T00:00:00Z"})
    store.add_candidates("r2", [_record("b", 1)])
    with pytest.raises(StoreError, match=r"occurs in 2 runs \['r1', 'r2'\]; pass run_id"):
        store.set_status("b", "triaged", note="x", author="t")
    _triage(store, "b", run_id="r2")
    assert store.get_candidate("b", run_id="r1")["status"] == "new"
    assert store.get_candidate("b", run_id="r2")["status"] == "triaged"
    assert store.get_candidate("a", sample_id="s1")["run_id"] == "r1"  # unique: no run needed


def test_concurrent_status_change_detected(store, monkeypatch):
    resolve = store._resolve

    def stale(*args):
        row = resolve(*args)
        with store._conn:  # another writer changes the status after we read it
            store._conn.execute("UPDATE candidates SET status = 'triaged' WHERE source_uid = 'b'")
        return row

    monkeypatch.setattr(store, "_resolve", stale)
    with pytest.raises(StoreError, match="concurrently"):
        store.set_status("b", "triaged", note="x", author="t")
    monkeypatch.undo()
    assert store.get_candidate("b")["status_history"] == []


def test_status_lifecycle(store):
    cand = _triage(store, sample_id="s1")
    assert cand["status"] == "triaged"
    assert cand["status_history"][0]["old_status"] == "new"
    assert cand["status_history"][0]["author"] == "tester"
    with pytest.raises(StoreError, match="vetting note"):
        store.set_status("a", "artifact", note="spike", author="t", sample_id="s1")
    note_id = _vet(store, sample_id="s1")
    assert note_id > 0
    cand = store.set_status("a", "unexplained", note="tests pass", author="t", sample_id="s1")
    cand = store.set_status("a", "followup", note="needs NIRSpec", author="t", sample_id="s1")
    assert [h["new_status"] for h in cand["status_history"]] == [
        "triaged",
        "unexplained",
        "followup",
    ]
    assert cand["vetting_notes"][0]["test"] == "edge_check"
    assert cand["vetting_notes"][0]["provenance"] == "observed"
    cand = store.set_status("a", "triaged", note="re-open", author="t", sample_id="s1")
    assert cand["status"] == "triaged"


@pytest.mark.parametrize("old", STATUSES)
def test_invalid_transitions_rejected(old):
    invalid = set(STATUSES) - TRANSITIONS[old]
    assert old in invalid  # no self-transitions
    assert "new" in invalid  # nothing returns to new


def test_invalid_transition_raises(store):
    with pytest.raises(StoreError, match="invalid transition new -> unexplained"):
        store.set_status("b", "unexplained", note="x", author="t")
    with pytest.raises(StoreError, match="status must be"):
        store.set_status("b", "discovered", note="x", author="t")
    with pytest.raises(StoreError, match="note"):
        store.set_status("b", "triaged", note="  ", author="t")
    assert store.get_candidate("b")["status"] == "new"


@pytest.mark.parametrize(
    "field, value, match",
    [
        ("outcome", "maybe", "outcome"),
        ("provenance", "rumour", "provenance"),
        ("evidence", "", "evidence"),
        ("test", " ", "test"),
        ("author", "", "author"),
    ],
)
def test_vetting_note_validation(store, field, value, match):
    kwargs = {
        "test": "t",
        "outcome": "fail",
        "evidence": "e",
        "provenance": "hypothesis",
        "author": "a",
    }
    kwargs[field] = value
    with pytest.raises(StoreError, match=match):
        store.add_vetting_note("b", **kwargs)


def test_to_python():
    assert to_python(np.float32(1.5)) == 1.5
    assert to_python(np.ma.masked) is None
    assert to_python(float("inf")) is None
    assert to_python(b"x") == "x"
    assert to_python({"a": (np.int64(1), np.array([1.0, np.nan]))}) == {"a": [1, [1.0, None]]}
