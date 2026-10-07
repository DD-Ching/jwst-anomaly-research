"""Shared pytest configuration.

Tests marked ``@pytest.mark.network`` hit live services (MAST, CDS, NED, ...). They are
skipped by default so CI stays offline and deterministic; enable them with
``pytest --run-network`` or ``JWST_ANOMALY_NETWORK=1``.
"""

from __future__ import annotations

import os

import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--run-network", action="store_true", default=False, help="run tests that need internet"
    )


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    if config.getoption("--run-network") or os.environ.get("JWST_ANOMALY_NETWORK") == "1":
        return
    skip = pytest.mark.skip(reason="network test; use --run-network or JWST_ANOMALY_NETWORK=1")
    for item in items:
        if "network" in item.keywords:
            item.add_marker(skip)
