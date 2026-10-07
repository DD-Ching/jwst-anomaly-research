"""Run provenance capture (code version, config, environment). Owner: bootstrap unit 6.

Everything here uses the standard library (``subprocess`` for git, ``hashlib``,
``importlib.metadata``, ``platform``); see DECISIONS.md D-007.
"""

from __future__ import annotations

import hashlib
import os
import platform
import re
import subprocess
import sys
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Any

from jwst_anomaly import __version__, paths

# Distributions whose versions can change results. Missing ones are recorded as None.
KEY_PACKAGES = (
    "jwst-anomaly",
    "astropy",
    "astroquery",
    "numpy",
    "scipy",
    "pandas",
    "pyarrow",
    "scikit-learn",
    "pyyaml",
    "matplotlib",
    "fsspec",
    "s3fs",
    "photutils",
    "jwst",
)

# Run ids sort chronologically: UTC timestamp, then 8 hex chars of a hash of the run inputs.
RUN_ID_PATTERN = re.compile(r"^\d{8}T\d{6}Z-[0-9a-f]{8}$")

_GIT_TIMEOUT_S = 10


def utc_now() -> datetime:
    """Current time as an aware UTC datetime."""
    return datetime.now(UTC)


def format_utc(dt: datetime) -> str:
    """ISO 8601 UTC timestamp with second precision and a ``Z`` suffix."""
    return dt.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def make_run_id(created: datetime, *parts: str) -> str:
    """Run id ``YYYYMMDDTHHMMSSZ-<8 hex>``; the hash covers the full timestamp and ``parts``."""
    created = created.astimezone(UTC)
    digest = hashlib.sha256("|".join([created.isoformat(), *parts]).encode()).hexdigest()
    return f"{created:%Y%m%dT%H%M%SZ}-{digest[:8]}"


def file_sha256(path: str | Path) -> str:
    """Hex sha256 of a file's bytes."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _git(repo_dir: Path, *args: str) -> str | None:
    try:
        out = subprocess.run(
            ["git", "-C", str(repo_dir), *args],
            capture_output=True,
            encoding="utf-8",  # git emits UTF-8 regardless of the console code page
            errors="replace",
            timeout=_GIT_TIMEOUT_S,
            check=True,
        )
    except (OSError, ValueError, subprocess.SubprocessError):
        return None
    return (out.stdout or "").strip()


_MANIFESTS = "data/manifests"
_NO_GIT = {
    "available": False,
    "commit": None,
    "branch": None,
    "dirty": None,
    "manifests_dirty": None,
}


def git_state(repo_dir: str | Path | None = None) -> dict[str, Any]:
    """Commit, branch and dirty flags of the code checkout; all None outside a git work tree.

    ``repo_dir`` (default: the repository containing this package) must be the work-tree root,
    so a package installed somewhere inside an unrelated checkout is not credited with that
    checkout's commit. ``dirty`` covers code and config (untracked, non-ignored files included:
    they may be code that ran); acquisition manifests are reported as ``manifests_dirty``.
    """
    repo = Path(repo_dir) if repo_dir is not None else paths.repo_root()
    top = _git(repo, "rev-parse", "--show-toplevel")
    if not top or Path(top).resolve() != repo.resolve():
        return dict(_NO_GIT)
    commit = _git(repo, "rev-parse", "HEAD")
    if not commit:
        return dict(_NO_GIT)
    code = _git(repo, "status", "--porcelain", "--", ".", f":(exclude){_MANIFESTS}")
    manifests = _git(repo, "status", "--porcelain", "--", _MANIFESTS)
    return {
        "available": True,
        "commit": commit,
        "branch": _git(repo, "rev-parse", "--abbrev-ref", "HEAD"),
        "dirty": None if code is None else bool(code),
        "manifests_dirty": None if manifests is None else bool(manifests),
    }


def package_versions(names: tuple[str, ...] = KEY_PACKAGES) -> dict[str, str | None]:
    """Installed distribution versions via ``importlib.metadata``; None if not installed."""
    versions: dict[str, str | None] = {}
    for name in names:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def capture_run_context(
    config_path: str | Path | None = None,
    *,
    repo_dir: str | Path | None = None,
    now: datetime | None = None,
) -> dict[str, Any]:
    """Return git commit/dirty state, config hash, package versions and UTC time for a run.

    The result is JSON-serializable. ``run_id`` follows :data:`RUN_ID_PATTERN`.
    ``repo_dir`` (default: the repository containing this package) and ``now`` exist for tests.
    """
    created = (now or utc_now()).astimezone(UTC)
    config: dict[str, Any] | None = None
    if config_path is not None:
        cpath = Path(config_path)
        config = {
            "path": cpath.resolve().as_posix(),
            "sha256": file_sha256(cpath),
            "size": cpath.stat().st_size,
        }
    git = git_state(repo_dir)
    run_id = make_run_id(
        created,
        config["sha256"] if config else "",
        git["commit"] or "",
        str(os.getpid()),
    )
    return {
        "run_id": run_id,
        "created_utc": format_utc(created),
        "config": config,
        "git": git,
        "jwst_anomaly_version": __version__,
        "python": {
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "platform": platform.platform(),
        },
        "packages": package_versions(),
        "argv": list(sys.argv),
    }
