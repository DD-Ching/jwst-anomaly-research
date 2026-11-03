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


#: Per-test wall-clock cap for network tests, in seconds (ASSUMPTION: a healthy live test takes
#: well under 5 min; the slowest seen, the MAST product listing, took about 2.5 min). A hung
#: service then fails one named test instead of stalling the whole job until the workflow's
#: job timeout cancels it with no test summary (#374, run 41476680422). Override with
#: ``JWST_ANOMALY_NETWORK_TIMEOUT`` (seconds; 0 disables).
NETWORK_TEST_TIMEOUT_S = 300.0


def network_test_timeout_s() -> float:
    raw = os.environ.get("JWST_ANOMALY_NETWORK_TIMEOUT", "")
    if not raw.strip():
        return NETWORK_TEST_TIMEOUT_S
    value = float(raw)
    if value < 0:
        raise ValueError(f"JWST_ANOMALY_NETWORK_TIMEOUT must be >= 0, got {raw!r}")
    return value


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    if config.getoption("--run-network") or os.environ.get("JWST_ANOMALY_NETWORK") == "1":
        timeout_s = network_test_timeout_s()
        if timeout_s > 0 and config.pluginmanager.hasplugin("timeout"):
            for item in items:
                if "network" in item.keywords and item.get_closest_marker("timeout") is None:
                    item.add_marker(pytest.mark.timeout(timeout_s))
        return
    skip = pytest.mark.skip(reason="network test; use --run-network or JWST_ANOMALY_NETWORK=1")
    for item in items:
        if "network" in item.keywords:
            item.add_marker(skip)
