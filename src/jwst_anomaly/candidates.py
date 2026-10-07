"""Persistent candidate store with vetting status. Owner: bootstrap unit 6.

One SQLite file (stdlib ``sqlite3``; DECISIONS.md D-007) holds:

* ``runs`` -- one row per pipeline run with its provenance context (JSON) and stage summary;
* ``candidates`` -- top-ranked sources of a run, unique per ``(run_id, sample_id, source_uid)``,
  with score/rank (``model_prediction``), per-method scores, image quality flags (``derived``),
  external cross-match summary (``observed``), cutout paths and a vetting ``status``;
* ``status_history`` -- every status change with a note, author and UTC time;
* ``vetting_notes`` -- one row per vetting test: name, outcome, evidence, provenance label.

The schema version is ``PRAGMA user_version`` (:data:`SCHEMA_VERSION`). A store written by
newer code is refused rather than silently misread.

Status lifecycle (:data:`TRANSITIONS`)::

    new -> triaged -> {artifact, known_object, explained, unexplained}
    unexplained -> followup
    any outcome or followup -> triaged   (re-open when new evidence arrives)
    followup -> {artifact, known_object, explained, unexplained}

Moving to an outcome status or ``followup`` requires at least one vetting note, so every
conclusion lists the tests behind it (docs/methodology.md).
"""

from __future__ import annotations

import json
import math
import sqlite3
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import numpy as np

from jwst_anomaly.provenance import format_utc, utc_now
from jwst_anomaly.schema import Provenance

SCHEMA_VERSION = 1

STATUSES = (
    "new",
    "triaged",
    "artifact",
    "known_object",
    "explained",
    "unexplained",
    "followup",
)
OUTCOME_STATUSES = frozenset({"artifact", "known_object", "explained", "unexplained"})
TRANSITIONS: dict[str, frozenset[str]] = {
    "new": frozenset({"triaged"}),
    "triaged": OUTCOME_STATUSES,
    "artifact": frozenset({"triaged"}),
    "known_object": frozenset({"triaged"}),
    "explained": frozenset({"triaged"}),
    "unexplained": frozenset({"triaged", "followup"}),
    "followup": OUTCOME_STATUSES | {"triaged"},
}
# Statuses that state a conclusion and therefore need recorded vetting tests.
REQUIRES_VETTING = OUTCOME_STATUSES | {"followup"}

# Vetting outcome of one test, from the candidate's point of view:
# pass = the ordinary explanation this test checks is not supported (candidate survives);
# fail = the test supports an ordinary explanation (artifact, known object, ...);
# inconclusive = the test could not decide.
OUTCOMES = ("pass", "fail", "inconclusive")

RUN_STATUSES = ("running", "completed", "failed")

_SCHEMA_V1 = """
CREATE TABLE runs (
    run_id        TEXT PRIMARY KEY,
    created_utc   TEXT NOT NULL,
    finished_utc  TEXT,
    status        TEXT NOT NULL,
    config_name   TEXT,
    config_sha256 TEXT,
    git_commit    TEXT,
    git_dirty     INTEGER,
    context_json  TEXT NOT NULL,
    summary_json  TEXT,
    error         TEXT
);
CREATE TABLE candidates (
    candidate_id       INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id             TEXT NOT NULL REFERENCES runs(run_id) ON DELETE CASCADE,
    sample_id          TEXT NOT NULL DEFAULT '',
    role               TEXT,
    source_uid         TEXT NOT NULL,
    ra                 REAL,
    dec                REAL,
    score              REAL,
    rank               INTEGER,
    top_features       TEXT,
    method_scores_json TEXT,
    flags_json         TEXT,
    xmatch_json        TEXT,
    cutouts_json       TEXT,
    status             TEXT NOT NULL DEFAULT 'new',
    created_utc        TEXT NOT NULL,
    updated_utc        TEXT NOT NULL,
    UNIQUE (run_id, sample_id, source_uid)
);
CREATE INDEX candidates_source_uid ON candidates(source_uid);
CREATE INDEX candidates_status ON candidates(status);
CREATE TABLE status_history (
    change_id    INTEGER PRIMARY KEY AUTOINCREMENT,
    candidate_id INTEGER NOT NULL REFERENCES candidates(candidate_id) ON DELETE CASCADE,
    old_status   TEXT NOT NULL,
    new_status   TEXT NOT NULL,
    note         TEXT NOT NULL,
    author       TEXT NOT NULL,
    changed_utc  TEXT NOT NULL
);
CREATE TABLE vetting_notes (
    note_id      INTEGER PRIMARY KEY AUTOINCREMENT,
    candidate_id INTEGER NOT NULL REFERENCES candidates(candidate_id) ON DELETE CASCADE,
    test         TEXT NOT NULL,
    outcome      TEXT NOT NULL,
    evidence     TEXT NOT NULL,
    provenance   TEXT NOT NULL,
    author       TEXT NOT NULL,
    created_utc  TEXT NOT NULL
);
"""

_JSON_FIELDS = {
    "method_scores_json": "method_scores",
    "flags_json": "flags",
    "xmatch_json": "xmatch",
    "cutouts_json": "cutouts",
    "context_json": "context",
    "summary_json": "summary",
}


class StoreError(Exception):
    """Invalid request to the candidate store (unknown candidate, bad transition, ...)."""


def to_python(value: Any) -> Any:
    """Convert numpy/astropy scalars and containers to JSON-safe Python; NaN and masked -> None."""
    if value is None or value is np.ma.masked:
        return None
    if hasattr(value, "unit") and hasattr(value, "value"):  # astropy Quantity
        value = value.value
    if isinstance(value, Mapping):
        return {str(k): to_python(v) for k, v in value.items()}
    if isinstance(value, list | tuple):
        return [to_python(v) for v in value]
    if isinstance(value, np.ndarray):
        return [to_python(v) for v in value.tolist()]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, bool | int | float | str):
        return value
    return str(value)


def _dumps(value: Any) -> str | None:
    return None if value is None else json.dumps(to_python(value), sort_keys=True)


def _row_dict(row: sqlite3.Row) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key in row.keys():
        value = row[key]
        if key in _JSON_FIELDS:
            out[_JSON_FIELDS[key]] = None if value is None else json.loads(value)
        else:
            out[key] = value
    if "git_dirty" in out and out["git_dirty"] is not None:
        out["git_dirty"] = bool(out["git_dirty"])
    return out


def _require_text(name: str, value: str | None) -> str:
    if value is None or not str(value).strip():
        raise StoreError(f"{name} must be a non-empty string")
    return str(value).strip()


class CandidateStore:
    """SQLite-backed store of runs, ranked candidates and vetting notes.

    ``create=False`` refuses to create a new file (for read/modify commands).
    Usable as a context manager; :meth:`close` releases the connection.
    """

    def __init__(self, path: str | Path, *, create: bool = True) -> None:
        self.path = Path(path)
        if not self.path.exists():
            if not create:
                raise FileNotFoundError(f"no candidate store at {self.path}")
            self.path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(self.path, timeout=30)
        self._conn.row_factory = sqlite3.Row
        self._conn.execute("PRAGMA foreign_keys = ON")
        try:
            self._init_schema()
        except Exception:
            self._conn.close()
            raise

    # -- lifecycle -------------------------------------------------------------------------

    def _init_schema(self) -> None:
        version = self.schema_version
        if version > SCHEMA_VERSION:
            raise StoreError(
                f"{self.path} has schema version {version}; this code supports "
                f"<= {SCHEMA_VERSION}. Upgrade jwst-anomaly."
            )
        if version == 0:
            tables = self._conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
            if tables.fetchall():
                raise StoreError(f"{self.path} is not a candidate store (no schema version)")
            # One transaction: tables and version appear together or not at all.
            self._conn.executescript(
                f"BEGIN;\n{_SCHEMA_V1}\nPRAGMA user_version = {SCHEMA_VERSION};\nCOMMIT;"
            )
        # Future migrations: ``if version < 2: ...`` here, each bumping user_version.

    @property
    def schema_version(self) -> int:
        """``PRAGMA user_version`` of the database file."""
        return int(self._conn.execute("PRAGMA user_version").fetchone()[0])

    def close(self) -> None:
        self._conn.close()

    def __enter__(self) -> CandidateStore:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    # -- runs ------------------------------------------------------------------------------

    def add_run(
        self,
        run_id: str,
        context: Mapping[str, Any],
        *,
        status: str = "running",
        config_name: str | None = None,
    ) -> None:
        """Register a run with its provenance context (see ``provenance.capture_run_context``)."""
        if status not in RUN_STATUSES:
            raise StoreError(f"run status must be one of {RUN_STATUSES}, got {status!r}")
        config = context.get("config") or {}
        git = context.get("git") or {}
        dirty = git.get("dirty")
        with self._conn:
            self._conn.execute(
                "INSERT INTO runs (run_id, created_utc, status, config_name, config_sha256,"
                " git_commit, git_dirty, context_json) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    run_id,
                    context.get("created_utc") or format_utc(utc_now()),
                    status,
                    config_name,
                    config.get("sha256"),
                    git.get("commit"),
                    None if dirty is None else int(bool(dirty)),
                    _dumps(context),
                ),
            )

    def finish_run(
        self,
        run_id: str,
        status: str,
        *,
        summary: Mapping[str, Any] | None = None,
        error: str | None = None,
    ) -> None:
        """Mark a run ``completed`` or ``failed`` and store its stage summary."""
        if status not in RUN_STATUSES:
            raise StoreError(f"run status must be one of {RUN_STATUSES}, got {status!r}")
        with self._conn:
            cur = self._conn.execute(
                "UPDATE runs SET status = ?, finished_utc = ?, summary_json = ?, error = ?"
                " WHERE run_id = ?",
                (status, format_utc(utc_now()), _dumps(summary), error, run_id),
            )
        if cur.rowcount == 0:
            raise StoreError(f"unknown run {run_id!r}")

    def get_run(self, run_id: str) -> dict[str, Any] | None:
        row = self._conn.execute("SELECT * FROM runs WHERE run_id = ?", (run_id,)).fetchone()
        return None if row is None else _row_dict(row)

    def list_runs(self, *, limit: int | None = None) -> list[dict[str, Any]]:
        """Runs, newest first (without the bulky context/summary JSON)."""
        sql = (
            "SELECT run_id, created_utc, finished_utc, status, config_name, config_sha256,"
            " git_commit, git_dirty, error,"
            " (SELECT COUNT(*) FROM candidates c WHERE c.run_id = runs.run_id) AS n_candidates"
            " FROM runs ORDER BY created_utc DESC, rowid DESC"
        )
        params: tuple[Any, ...] = ()
        if limit is not None:
            sql += " LIMIT ?"
            params = (int(limit),)
        return [_row_dict(r) for r in self._conn.execute(sql, params)]

    def latest_run_id(self, *, status: str | None = None) -> str | None:
        sql = "SELECT run_id FROM runs"
        params: tuple[Any, ...] = ()
        if status is not None:
            sql += " WHERE status = ?"
            params = (status,)
        row = self._conn.execute(sql + " ORDER BY created_utc DESC, rowid DESC", params).fetchone()
        return None if row is None else row[0]

    # -- candidates ------------------------------------------------------------------------

    def add_candidates(self, run_id: str, records: Iterable[Mapping[str, Any]]) -> int:
        """Insert candidate records of a run; return the number inserted.

        Each record needs ``source_uid`` and may carry ``sample_id``, ``role``, ``ra``, ``dec``,
        ``score``, ``rank``, ``top_features``, ``method_scores`` (``{method: score}``),
        ``flags`` (``{band: {quality_flag, on_edge, frac_nan}}``), ``xmatch`` (one
        ``schema.XMATCH_COLUMNS`` row as a dict) and ``cutouts`` (``{band: path}``).
        New candidates start in status ``new``.
        """
        now = format_utc(utc_now())
        rows = []
        for rec in records:
            uid = _require_text("source_uid", to_python(rec.get("source_uid")))
            rank = to_python(rec.get("rank"))
            rows.append(
                (
                    run_id,
                    str(rec.get("sample_id") or ""),
                    rec.get("role"),
                    uid,
                    to_python(rec.get("ra")),
                    to_python(rec.get("dec")),
                    to_python(rec.get("score")),
                    None if rank is None else int(rank),
                    to_python(rec.get("top_features")),
                    _dumps(rec.get("method_scores")),
                    _dumps(rec.get("flags")),
                    _dumps(rec.get("xmatch")),
                    _dumps(rec.get("cutouts")),
                    now,
                    now,
                )
            )
        try:
            with self._conn:
                self._conn.executemany(
                    "INSERT INTO candidates (run_id, sample_id, role, source_uid, ra, dec, score,"
                    " rank, top_features, method_scores_json, flags_json, xmatch_json,"
                    " cutouts_json, created_utc, updated_utc)"
                    " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    rows,
                )
        except sqlite3.IntegrityError as exc:
            raise StoreError(f"cannot add candidates to run {run_id!r}: {exc}") from exc
        return len(rows)

    def list_candidates(
        self,
        *,
        run_id: str | None = None,
        status: str | None = None,
        sample_id: str | None = None,
        limit: int | None = None,
    ) -> list[dict[str, Any]]:
        """Candidates ordered by run (newest first), sample and rank."""
        if status is not None and status not in STATUSES:
            raise StoreError(f"status must be one of {STATUSES}, got {status!r}")
        where, params = [], []
        for col, val in (("c.run_id", run_id), ("c.status", status), ("c.sample_id", sample_id)):
            if val is not None:
                where.append(f"{col} = ?")
                params.append(val)
        sql = "SELECT c.* FROM candidates c JOIN runs r ON r.run_id = c.run_id"
        if where:
            sql += " WHERE " + " AND ".join(where)
        sql += " ORDER BY r.created_utc DESC, r.rowid DESC, c.sample_id, c.rank, c.candidate_id"
        if limit is not None:
            sql += " LIMIT ?"
            params.append(int(limit))
        return [_row_dict(r) for r in self._conn.execute(sql, params)]

    def _resolve(self, source_uid: str, run_id: str | None, sample_id: str | None) -> sqlite3.Row:
        """Find exactly one candidate row.

        Without ``run_id`` the uid must occur in a single run: source uids come from per-run
        catalogs and need not denote the same object in another run, so the store never guesses.
        """
        where = f"{source_uid!r}" + (f" (sample {sample_id!r})" if sample_id is not None else "")
        if run_id is None:
            runs = [
                r[0]
                for r in self._conn.execute(
                    "SELECT DISTINCT c.run_id FROM candidates c"
                    " WHERE c.source_uid = ? AND (? IS NULL OR c.sample_id = ?)",
                    (source_uid, sample_id, sample_id),
                )
            ]
            if not runs:
                raise StoreError(f"no candidate {where} in any run")
            if len(runs) > 1:
                raise StoreError(
                    f"candidate {where} occurs in {len(runs)} runs {sorted(runs)}; "
                    "pass run_id (CLI: --run)"
                )
            run_id = runs[0]
        rows = self._conn.execute(
            "SELECT * FROM candidates WHERE run_id = ? AND source_uid = ?"
            " AND (? IS NULL OR sample_id = ?)",
            (run_id, source_uid, sample_id, sample_id),
        ).fetchall()
        if not rows:
            raise StoreError(f"no candidate {where} in run {run_id!r}")
        if len(rows) > 1:
            samples = sorted(r["sample_id"] for r in rows)
            raise StoreError(
                f"source_uid {source_uid!r} is ambiguous in run {run_id!r} "
                f"(samples {samples}); pass sample_id (CLI: --sample)"
            )
        return rows[0]

    def get_candidate(
        self, source_uid: str, *, run_id: str | None = None, sample_id: str | None = None
    ) -> dict[str, Any]:
        """One candidate with its vetting notes and status history."""
        row = self._resolve(source_uid, run_id, sample_id)
        cand = _row_dict(row)
        cid = row["candidate_id"]
        cand["vetting_notes"] = [
            _row_dict(r)
            for r in self._conn.execute(
                "SELECT * FROM vetting_notes WHERE candidate_id = ? ORDER BY note_id", (cid,)
            )
        ]
        cand["status_history"] = [
            _row_dict(r)
            for r in self._conn.execute(
                "SELECT * FROM status_history WHERE candidate_id = ? ORDER BY change_id", (cid,)
            )
        ]
        return cand

    def set_status(
        self,
        source_uid: str,
        status: str,
        *,
        note: str,
        author: str,
        run_id: str | None = None,
        sample_id: str | None = None,
    ) -> dict[str, Any]:
        """Move a candidate along :data:`TRANSITIONS`; return the updated candidate."""
        if status not in STATUSES:
            raise StoreError(f"status must be one of {STATUSES}, got {status!r}")
        note = _require_text("note", note)
        author = _require_text("author", author)
        row = self._resolve(source_uid, run_id, sample_id)
        old = row["status"]
        if status not in TRANSITIONS[old]:
            allowed = sorted(TRANSITIONS[old])
            raise StoreError(f"invalid transition {old} -> {status}; allowed from {old}: {allowed}")
        cid = row["candidate_id"]
        if status in REQUIRES_VETTING:
            n_notes = self._conn.execute(
                "SELECT COUNT(*) FROM vetting_notes WHERE candidate_id = ?", (cid,)
            ).fetchone()[0]
            if n_notes == 0:
                raise StoreError(
                    f"status {status!r} states a conclusion; add at least one vetting note "
                    "(add_vetting_note / `candidates add-vetting`) first"
                )
        now = format_utc(utc_now())
        with self._conn:
            cur = self._conn.execute(
                "UPDATE candidates SET status = ?, updated_utc = ?"
                " WHERE candidate_id = ? AND status = ?",
                (status, now, cid, old),
            )
            if cur.rowcount != 1:  # someone else changed it since we read it
                raise StoreError(f"status of {source_uid!r} changed concurrently; retry")
            self._conn.execute(
                "INSERT INTO status_history (candidate_id, old_status, new_status, note, author,"
                " changed_utc) VALUES (?, ?, ?, ?, ?, ?)",
                (cid, old, status, note, author, now),
            )
        return self.get_candidate(source_uid, run_id=row["run_id"], sample_id=row["sample_id"])

    def add_vetting_note(
        self,
        source_uid: str,
        *,
        test: str,
        outcome: str,
        evidence: str,
        provenance: str,
        author: str,
        run_id: str | None = None,
        sample_id: str | None = None,
    ) -> int:
        """Record one vetting test for a candidate; return the note id.

        ``outcome`` is one of :data:`OUTCOMES`; ``provenance`` is a ``schema.Provenance``
        value describing the evidence (e.g. ``observed`` for an archive image,
        ``hypothesis`` for a proposed interpretation).
        """
        if outcome not in OUTCOMES:
            raise StoreError(f"outcome must be one of {OUTCOMES}, got {outcome!r}")
        valid_prov = [p.value for p in Provenance]
        if provenance not in valid_prov:
            raise StoreError(f"provenance must be one of {valid_prov}, got {provenance!r}")
        test = _require_text("test", test)
        evidence = _require_text("evidence", evidence)
        author = _require_text("author", author)
        row = self._resolve(source_uid, run_id, sample_id)
        now = format_utc(utc_now())
        with self._conn:
            cur = self._conn.execute(
                "INSERT INTO vetting_notes (candidate_id, test, outcome, evidence, provenance,"
                " author, created_utc) VALUES (?, ?, ?, ?, ?, ?, ?)",
                (row["candidate_id"], test, outcome, evidence, provenance, author, now),
            )
            self._conn.execute(
                "UPDATE candidates SET updated_utc = ? WHERE candidate_id = ?",
                (now, row["candidate_id"]),
            )
        return int(cur.lastrowid)
