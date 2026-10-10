"""Offline tests for the durable derived-data store (S-1)."""

from __future__ import annotations

import hashlib

import pytest
from astropy.table import Table

from jwst_anomaly import derived_store as ds


def _manifest(tmp_path, rows):
    t = Table(
        rows=rows, names=("name", "tag", "sha256", "size", "sources", "code_commit", "created")
    )
    p = tmp_path / "derived_data.ecsv"
    t.write(p, format="ascii.ecsv")
    return p


class _Resp:
    def __init__(self, content, status=200):
        self.content, self.status_code = content, status


def test_pinned_takes_the_newest_tag(tmp_path):
    m = _manifest(
        tmp_path,
        [
            ("a.npz", "derived-data-20261001", "x" * 64, 1, "u", "c", "d"),
            ("a.npz", "derived-data-20261010", "y" * 64, 1, "u", "c", "d"),
        ],
    )
    assert ds.pinned("a.npz", m)["tag"] == "derived-data-20261010"
    assert ds.pinned("b.npz", m) is None


def test_fetch_verifies_sha256(tmp_path, monkeypatch):
    good = b"reduced maps"
    sha = hashlib.sha256(good).hexdigest()
    m = _manifest(tmp_path, [("a.npz", "derived-data-20261010", sha, len(good), "u", "c", "d")])
    import requests

    monkeypatch.setattr(requests, "get", lambda *a, **k: _Resp(good))
    dest = tmp_path / "out" / "a.npz"
    assert ds.fetch("a.npz", dest, m) == dest
    assert dest.read_bytes() == good
    # an existing matching file is reused without a request
    monkeypatch.setattr(requests, "get", lambda *a, **k: pytest.fail("no request expected"))
    assert ds.fetch("a.npz", dest, m) == dest


def test_fetch_refuses_wrong_bytes_and_errors(tmp_path, monkeypatch):
    m = _manifest(tmp_path, [("a.npz", "derived-data-20261010", "0" * 64, 3, "u", "c", "d")])
    import requests

    monkeypatch.setattr(requests, "get", lambda *a, **k: _Resp(b"bad"))
    dest = tmp_path / "a.npz"
    assert ds.fetch("a.npz", dest, m) is None
    assert not dest.exists()

    def boom(*a, **k):
        raise requests.ConnectionError("blocked")

    monkeypatch.setattr(requests, "get", boom)
    assert ds.fetch("a.npz", dest, m) is None


def test_publish_tags_the_code_commit_and_pins_only_after_success(tmp_path, monkeypatch):
    import subprocess
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import derived_publish as P

    f = tmp_path / "a.npz"
    f.write_bytes(b"maps")
    monkeypatch.setattr(ds, "MANIFEST", tmp_path / "derived_data.ecsv")
    calls = []

    def run(cmd, **kw):
        calls.append(cmd)
        if cmd[0] == "git":
            return subprocess.CompletedProcess(cmd, 0, stdout="abc123\n")
        if fail:
            raise subprocess.CalledProcessError(1, cmd)
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(P.subprocess, "run", run)
    fail = True
    with pytest.raises(subprocess.CalledProcessError):
        P.main([str(f), "--sources", "u", "--tag", "derived-data-20261010"])
    assert not ds.MANIFEST.exists()  # a refused release leaves no pin
    fail = False
    assert P.main([str(f), "--sources", "u", "--tag", "derived-data-20261010"]) == 0
    gh = calls[-1]
    assert gh[gh.index("--target") + 1] == "abc123"
    assert ds.pinned("a.npz", ds.MANIFEST)["code_commit"] == "abc123"


def test_manifest_row_refuses_large_files(tmp_path, monkeypatch):
    f = tmp_path / "big.bin"
    f.write_bytes(b"0" * 10)
    monkeypatch.setattr(ds, "MAX_BYTES", 5)
    with pytest.raises(ValueError):
        ds.manifest_row(f, "t", "s", "c", "d")
    monkeypatch.setattr(ds, "MAX_BYTES", 50)
    assert ds.manifest_row(f, "t", "s", "c", "d")["size"] == 10
