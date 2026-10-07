"""``jwst-anomaly`` command-line entry point. Owner: bootstrap unit 6.

::

    jwst-anomaly run --config configs/reference_sample.yaml [--sample ID ...]
    jwst-anomaly runs list
    jwst-anomaly candidates list [--run ID] [--status S] [--sample ID] [--json]
    jwst-anomaly candidates show <source_uid> [--run ID] [--sample ID]
    jwst-anomaly candidates set-status <source_uid> <status> --note TEXT [--run ID]
    jwst-anomaly candidates add-vetting <source_uid> --test NAME --outcome pass|fail|inconclusive
        --evidence TEXT --provenance LABEL [--run ID]

All commands take ``--db PATH`` (default ``$JWST_ANOMALY_OUTPUTS/candidates.sqlite``).
``--author`` defaults to ``$JWST_ANOMALY_AUTHOR``, then the login name.
Exit codes: 0 success, 1 failure (failed run, unknown candidate, invalid transition, ...),
2 usage error.
"""

from __future__ import annotations

import argparse
import getpass
import json
import logging
import os
import sqlite3
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from jwst_anomaly import __version__, paths, pipeline
from jwst_anomaly.candidates import OUTCOMES, STATUSES, CandidateStore, StoreError
from jwst_anomaly.schema import Provenance

AUTHOR_ENV = "JWST_ANOMALY_AUTHOR"

EXIT_OK, EXIT_FAILURE, EXIT_USAGE = 0, 1, 2


def _default_author() -> str:
    author = os.environ.get(AUTHOR_ENV)
    if author:
        return author
    try:
        return getpass.getuser()
    except Exception:  # no login name available (e.g. some containers)
        return "unknown"


def _add_db(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--db",
        type=Path,
        default=None,
        help="candidate store (default: <outputs>/candidates.sqlite)",
    )


def _add_target(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("source_uid")
    parser.add_argument(
        "--run", dest="run_id", help="run id (needed if the uid is in several runs)"
    )
    parser.add_argument("--sample", dest="sample_id", help="sample id, if the uid is ambiguous")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="jwst-anomaly",
        description="Rank unusual sources in public JWST data and track their vetting. "
        "Anomaly scores are rankings, not evidence of new physics.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    sub = parser.add_subparsers(dest="command", required=True)

    p_run = sub.add_parser("run", help="run the pipeline from a YAML config")
    p_run.add_argument("--config", type=Path, required=True)
    p_run.add_argument("--sample", action="append", dest="samples", help="only this sample id")
    p_run.add_argument("--outputs", type=Path, default=None, help="outputs directory override")
    _add_db(p_run)

    p_runs = sub.add_parser("runs", help="inspect recorded runs")
    runs_sub = p_runs.add_subparsers(dest="action", required=True)
    p_rl = runs_sub.add_parser("list", help="list runs, newest first")
    p_rl.add_argument("--limit", type=int, default=None)
    p_rl.add_argument("--json", action="store_true")
    _add_db(p_rl)

    p_cand = sub.add_parser("candidates", help="inspect and vet candidates")
    cand_sub = p_cand.add_subparsers(dest="action", required=True)

    p_list = cand_sub.add_parser("list", help="list candidates")
    p_list.add_argument("--run", dest="run_id")
    p_list.add_argument("--status", choices=STATUSES)
    p_list.add_argument("--sample", dest="sample_id")
    p_list.add_argument("--limit", type=int, default=None)
    p_list.add_argument("--json", action="store_true")
    _add_db(p_list)

    p_show = cand_sub.add_parser("show", help="one candidate with vetting notes and history")
    _add_target(p_show)
    _add_db(p_show)

    p_set = cand_sub.add_parser("set-status", help="move a candidate along the status lifecycle")
    _add_target(p_set)
    p_set.add_argument("status", choices=STATUSES)
    p_set.add_argument("--note", required=True, help="why the status changes")
    p_set.add_argument("--author", default=None)
    _add_db(p_set)

    p_vet = cand_sub.add_parser("add-vetting", help="record one vetting test")
    _add_target(p_vet)
    p_vet.add_argument("--test", required=True, help="e.g. edge_check, diffraction_spike, simbad")
    p_vet.add_argument(
        "--outcome",
        required=True,
        choices=OUTCOMES,
        help="pass = the ordinary explanation tested is not supported; fail = it is",
    )
    p_vet.add_argument("--evidence", required=True, help="what was looked at (paths, URLs, values)")
    p_vet.add_argument("--provenance", required=True, choices=[p.value for p in Provenance])
    p_vet.add_argument("--author", default=None)
    _add_db(p_vet)
    return parser


def _db_path(args: argparse.Namespace) -> Path:
    return args.db if args.db is not None else pipeline.default_db_path()


def _print_table(rows: list[dict[str, Any]], columns: Sequence[str]) -> None:
    cells = [["" if r.get(c) is None else str(r.get(c)) for c in columns] for r in rows]
    widths = [max([len(c), *(len(row[i]) for row in cells)]) for i, c in enumerate(columns)]
    print("  ".join(c.upper().ljust(w) for c, w in zip(columns, widths, strict=True)).rstrip())
    for row in cells:
        print("  ".join(v.ljust(w) for v, w in zip(row, widths, strict=True)).rstrip())


def _cmd_run(args: argparse.Namespace) -> int:
    try:
        run_id = pipeline.run(
            args.config, outputs_dir=args.outputs, db_path=args.db, samples=args.samples
        )
    except pipeline.PipelineError as exc:
        print(f"error: run {exc.run_id} failed: {exc}", file=sys.stderr)
        if exc.report_path is not None:
            print(f"report: {exc.report_path}", file=sys.stderr)
        return EXIT_FAILURE
    outputs = args.outputs if args.outputs is not None else paths.outputs_dir()
    db = args.db if args.db is not None else pipeline.default_db_path(outputs)
    with CandidateStore(db, create=False) as store:
        run = store.get_run(run_id) or {}
    n = (run.get("summary") or {}).get("n_candidates", 0)
    print(f"run {run_id} completed: {n} candidates stored in {db}")
    print(f"report: {Path(outputs) / 'runs' / run_id / 'report.md'}")
    return EXIT_OK


def _cmd_runs(args: argparse.Namespace) -> int:
    with CandidateStore(_db_path(args), create=False) as store:
        runs = store.list_runs(limit=args.limit)
    if args.json:
        print(json.dumps(runs, indent=2))
        return EXIT_OK
    for r in runs:
        r["git_commit"] = (r.get("git_commit") or "")[:10]
    _print_table(runs, ("run_id", "status", "config_name", "git_commit", "n_candidates", "error"))
    return EXIT_OK


def _cmd_candidates(args: argparse.Namespace) -> int:
    with CandidateStore(_db_path(args), create=False) as store:
        if args.action == "list":
            rows = store.list_candidates(
                run_id=args.run_id, status=args.status, sample_id=args.sample_id, limit=args.limit
            )
            if args.json:
                print(json.dumps(rows, indent=2))
                return EXIT_OK
            for r in rows:
                r["score"] = None if r["score"] is None else f"{r['score']:.4g}"
            _print_table(rows, ("run_id", "sample_id", "rank", "source_uid", "score", "status"))
            return EXIT_OK
        if args.action == "show":
            cand = store.get_candidate(
                args.source_uid, run_id=args.run_id, sample_id=args.sample_id
            )
            print(json.dumps(cand, indent=2))
            return EXIT_OK
        author = args.author or _default_author()
        if args.action == "set-status":
            cand = store.set_status(
                args.source_uid,
                args.status,
                note=args.note,
                author=author,
                run_id=args.run_id,
                sample_id=args.sample_id,
            )
            print(
                f"{cand['source_uid']} (run {cand['run_id']}, sample {cand['sample_id']}): "
                f"status {cand['status_history'][-1]['old_status']} -> {cand['status']}"
            )
            return EXIT_OK
        note_id = store.add_vetting_note(
            args.source_uid,
            test=args.test,
            outcome=args.outcome,
            evidence=args.evidence,
            provenance=args.provenance,
            author=author,
            run_id=args.run_id,
            sample_id=args.sample_id,
        )
        print(f"vetting note {note_id} added to {args.source_uid}: {args.test} -> {args.outcome}")
        return EXIT_OK


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``jwst-anomaly`` console script."""
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:  # --help/--version (0) or usage error (2)
        return exc.code if isinstance(exc.code, int) else EXIT_USAGE
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )
    handlers = {"run": _cmd_run, "runs": _cmd_runs, "candidates": _cmd_candidates}
    try:
        return handlers[args.command](args)
    except (pipeline.ConfigError, StoreError, sqlite3.Error, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_FAILURE


if __name__ == "__main__":
    raise SystemExit(main())
