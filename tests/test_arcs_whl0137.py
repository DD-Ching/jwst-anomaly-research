"""Tests for scripts/arcs_whl0137.py: the registered RELICS map models (no network)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("arcs_whl0137", _DIR / "arcs_whl0137.py")
aw = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(aw)

import lens_consistency as lc  # noqa: E402


def test_both_models_are_pinned_map_models_at_the_cluster_redshift():
    for name in (aw.PRIMARY, aw.SECOND):
        m = lc.MODELS[name]
        assert m["kind"] == "maps"
        assert m["z_lens"] == 0.566
        assert {"alpha_x", "alpha_y", "kappa_map"} <= set(m["files"])
        for url, sha in m["files"].values():
            assert url.startswith("https://archive.stsci.edu/hlsps/relics/whl0137m08/")
            assert len(sha) == 64
        for key in m.get("mag_maps", {}).values():
            assert key in m["files"]


def test_second_model_is_glafic():
    assert "glafic" in lc.MODELS[aw.SECOND]["files"]["alpha_x"][0]
