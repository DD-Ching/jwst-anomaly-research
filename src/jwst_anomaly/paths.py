"""Filesystem locations. Data never lives in git; manifests do."""

from __future__ import annotations

import os
from pathlib import Path

DATA_ENV = "JWST_ANOMALY_DATA"
OUTPUTS_ENV = "JWST_ANOMALY_OUTPUTS"


def repo_root() -> Path:
    """Repository root (valid for the editable/src-layout install used by this project)."""
    return Path(__file__).resolve().parents[2]


def data_root() -> Path:
    """Root for downloaded and derived data: ``$JWST_ANOMALY_DATA`` or ``<repo>/data``."""
    env = os.environ.get(DATA_ENV)
    return Path(env).expanduser() if env else repo_root() / "data"


def cache_dir() -> Path:
    """Downloaded archive products (gitignored)."""
    return data_root() / "cache"


def manifests_dir() -> Path:
    """Tracked, machine-readable acquisition manifests. Always inside the repo."""
    return repo_root() / "data" / "manifests"


def outputs_dir() -> Path:
    """Run outputs (scores, cutouts, reports, candidate DB).

    ``$JWST_ANOMALY_OUTPUTS`` or ``<repo>/outputs``.
    """
    env = os.environ.get(OUTPUTS_ENV)
    return Path(env).expanduser() if env else repo_root() / "outputs"
