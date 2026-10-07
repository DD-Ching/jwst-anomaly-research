import hashlib
import json
import subprocess
from datetime import UTC, datetime

from jwst_anomaly import __version__, paths, provenance


def test_capture_run_context(tmp_path):
    config = tmp_path / "c.yaml"
    config.write_bytes(b"name: x\n")
    now = datetime(2026, 10, 7, 12, 34, 56, 789, tzinfo=UTC)
    ctx = provenance.capture_run_context(config, now=now)

    json.dumps(ctx)  # JSON-serializable as stored
    assert provenance.RUN_ID_PATTERN.match(ctx["run_id"])
    assert ctx["run_id"].startswith("20261007T123456Z-")
    assert ctx["created_utc"] == "2026-10-07T12:34:56Z"
    assert ctx["config"]["sha256"] == hashlib.sha256(b"name: x\n").hexdigest()
    assert ctx["config"]["size"] == 8
    assert ctx["jwst_anomaly_version"] == __version__
    assert ctx["packages"]["astropy"]
    assert ctx["packages"]["jwst-anomaly"] == __version__
    assert ctx["python"]["version"].count(".") == 2
    assert set(ctx["git"]) == {"available", "commit", "branch", "dirty", "manifests_dirty"}


def test_git_state_of_this_checkout():
    git = provenance.git_state()
    if git["available"]:  # CI and dev checkouts; tolerated otherwise
        assert len(git["commit"]) == 40
        assert isinstance(git["dirty"], bool)


def test_non_git_directory_tolerated(tmp_path):
    ctx = provenance.capture_run_context(None, repo_dir=tmp_path)
    assert ctx["config"] is None
    assert ctx["git"]["available"] is False
    assert {ctx["git"][k] for k in ("commit", "branch", "dirty", "manifests_dirty")} == {None}
    assert provenance.RUN_ID_PATTERN.match(ctx["run_id"])


def test_subdirectory_of_checkout_not_credited():
    """A package dir inside some other checkout must not report that checkout's commit."""
    assert provenance.git_state(paths.repo_root() / "src")["available"] is False


def test_undecodable_git_output_tolerated(monkeypatch):
    class Done:
        stdout = None  # what subprocess leaves when its reader thread fails to decode

    monkeypatch.setattr(subprocess, "run", lambda *a, **k: Done())
    assert provenance.git_state()["available"] is False


def test_missing_git_binary_tolerated(monkeypatch):
    def no_git(*args, **kwargs):
        raise FileNotFoundError("git")

    monkeypatch.setattr(subprocess, "run", no_git)
    assert provenance.git_state()["commit"] is None


def test_git_timeout_tolerated(monkeypatch):
    def slow(*args, **kwargs):
        raise subprocess.TimeoutExpired(cmd="git", timeout=1)

    monkeypatch.setattr(subprocess, "run", slow)
    assert provenance.git_state()["available"] is False


def test_run_ids_unique_and_sortable():
    t1 = datetime(2026, 10, 7, 0, 0, 0, 1, tzinfo=UTC)
    t2 = datetime(2026, 10, 7, 0, 0, 0, 2, tzinfo=UTC)
    a, b = provenance.make_run_id(t1, "cfg"), provenance.make_run_id(t2, "cfg")
    assert a != b and a[:16] == b[:16]
    assert provenance.make_run_id(t1, "cfg") == a  # deterministic given inputs
    assert provenance.make_run_id(datetime(2026, 10, 8, tzinfo=UTC)) > a


def test_package_versions_missing():
    assert provenance.package_versions(("surely-not-installed-pkg",)) == {
        "surely-not-installed-pkg": None
    }
