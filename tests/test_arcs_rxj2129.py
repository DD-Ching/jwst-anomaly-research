"""Tests for scripts/arcs_rxj2129.py: zip-member staging, the validate gate and the frozen outcome
(local zip and synthetic summaries; no network)."""

from __future__ import annotations

import hashlib
import importlib.util
import io
import json
import sys
import zipfile
from pathlib import Path

import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("arcs_rxj2129", _DIR / "arcs_rxj2129.py")
ar = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ar)


def _zip(tmp_path: Path) -> tuple[Path, bytes]:
    data = b"SIMPLE  =" + bytes(range(256)) * 50
    path = tmp_path / "models.zip"
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("M/LENSTOOL/alphaX_DLSDS1.fits", data)
        z.writestr("M/other.fits", b"x")
    return path, data


def _local(monkeypatch, tmp_path, path: Path) -> None:
    monkeypatch.setattr(ar, "_RangeFile", lambda url: io.FileIO(path, "r"))
    monkeypatch.setattr(ar.paths, "cache_dir", lambda: tmp_path / "cache")


def test_stage_members_writes_verified_member_under_fetch_name(monkeypatch, tmp_path):
    path, data = _zip(tmp_path)
    _local(monkeypatch, tmp_path, path)
    sha = hashlib.sha256(data).hexdigest()
    ar.stage_members({"alpha_x": ("https://x/zip/content#M/LENSTOOL/alphaX_DLSDS1.fits", sha)})
    target = tmp_path / "cache" / "external" / f"{sha[:12]}_alphaX_DLSDS1.fits"
    assert target.read_bytes() == data
    assert not list((tmp_path / "cache" / "external").glob("*.part"))


def test_stage_members_refuses_a_checksum_mismatch(monkeypatch, tmp_path):
    path, _ = _zip(tmp_path)
    _local(monkeypatch, tmp_path, path)
    with pytest.raises(ValueError, match="sha256 mismatch"):
        ar.stage_members({"alpha_x": ("https://x#M/LENSTOOL/alphaX_DLSDS1.fits", "0" * 64)})
    assert not list((tmp_path / "cache" / "external").iterdir())


def test_failed_validate_leaves_an_empty_control_and_the_second_model_decides(
    monkeypatch, tmp_path
):
    def fake_main(argv):
        d = tmp_path / argv[1]
        d.mkdir(parents=True, exist_ok=True)
        rel = 0.5 if argv[1] == ar.PRIMARY else 1e-5
        (d / "validate.json").write_text(json.dumps({"kappa_map": {"median_rel_diff": rel}}))

    monkeypatch.setattr(ar.lc, "main", fake_main)
    second = {
        "model": ar.SECOND,
        "strong_shear_background": {"n": 30, "aligned": 25, "anti": 1},
        "anti_expected_shape_noise_background": {"p_ge_observed": 0.4},
    }
    monkeypatch.setattr(ar, "_arcs", lambda model, out, cat, pz: second)
    primary = ar._run(ar.PRIMARY, tmp_path, Path("c"), Path("p"))
    assert primary["validate_passes"] is False
    assert primary["strong_shear_background"] == {"n": 0}
    res = ar.outcome(primary, ar._run(ar.SECOND, tmp_path, Path("c"), Path("p")))
    assert res == {
        "deciding_model": ar.SECOND,
        "outcome": "does not replicate",
        "p_ge_observed": 0.4,
    }
