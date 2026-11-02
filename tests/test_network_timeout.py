"""The per-test timeout that conftest.py gives network tests (#374)."""

from __future__ import annotations

import pytest
from conftest import NETWORK_TEST_TIMEOUT_S, network_test_timeout_s


def test_default_timeout(monkeypatch):
    monkeypatch.delenv("JWST_ANOMALY_NETWORK_TIMEOUT", raising=False)
    assert network_test_timeout_s() == NETWORK_TEST_TIMEOUT_S


@pytest.mark.parametrize(("raw", "expected"), [("120", 120.0), ("0", 0.0), (" ", 300.0)])
def test_timeout_override(monkeypatch, raw, expected):
    monkeypatch.setenv("JWST_ANOMALY_NETWORK_TIMEOUT", raw)
    assert network_test_timeout_s() == expected


def test_negative_timeout_rejected(monkeypatch):
    monkeypatch.setenv("JWST_ANOMALY_NETWORK_TIMEOUT", "-1")
    with pytest.raises(ValueError, match=">= 0"):
        network_test_timeout_s()
