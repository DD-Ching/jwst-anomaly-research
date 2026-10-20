"""Tests for scripts/arcs_rxc0018.py: model registration and wiring (no network)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_model_registered_as_planned():
    lc = _load("lens_consistency")
    ar = _load("arcs_rxc0018")
    spec = lc.MODELS[ar.PRIMARY]
    assert spec["kind"] == "maps"
    assert spec["z_lens"] == 0.546  # Zitrin et al. 2011
    assert set(spec["files"]) == {"alpha_x", "alpha_y", "kappa_map", "mag_z2"}
    for url, sha in spec["files"].values():
        assert url.startswith(
            "https://archive.stsci.edu/hlsps/relics/rxc0018p16/models/lenstool/v1/"
        )
        assert len(sha) == 64


def test_main_passes_the_one_model_as_primary_and_second(tmp_path, monkeypatch):
    ar = _load("arcs_rxc0018")
    calls = {}
    monkeypatch.setattr(ar, "decide", lambda *a: calls.setdefault("decide", a))
    assert ar.main(["--catalog", "c.ecsv", "--photoz", "z.fits", "--out", str(tmp_path)]) == 0
    primary, second, out, catalog, photoz, plan = calls["decide"]
    assert primary == second == "rxc0018-relics-lenstool"
    assert (out, catalog, photoz) == (tmp_path, Path("c.ecsv"), Path("z.fits"))
    assert plan == "results/arcs/rxc0018-relics/PLAN.md"


def test_decide_runs_a_repeated_model_once(tmp_path, monkeypatch):
    rx = _load("arcs_rxj2129")
    ran = []

    def fake_run(model, out, catalog, photoz):
        ran.append(model)
        return {"model": model, "validate_passes": False, "strong_shear_background": {"n": 0}}

    monkeypatch.setattr(rx, "_run", fake_run)
    monkeypatch.setattr(rx, "sha256_file", lambda p: "0" * 64)
    res = rx.decide("m", "m", tmp_path, Path("c"), Path("z"), "plan")
    assert ran == ["m"]
    assert res["outcome"] == "no result"
