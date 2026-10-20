"""Offline tests for scripts/m1149_sharon.py (frozen gates, transcribed image list; no network)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import lens_consistency as lc  # noqa: E402
import m1149_sharon as ms  # noqa: E402

GOLD = ROOT / "results/lens_consistency/macs1149_sharon/sharon_gold_images.ecsv"


def _maps(kappa: float, mu: float) -> dict:
    return {"kappa": {"median_rel_diff": kappa}, "mu_z2": {"median_rel_diff": mu}}


def _ip(systems: list[str], d: list[float]) -> Table:
    return Table({"system": systems, "dtheta_arcsec": np.array(d, float)})


def test_gold_images_table():
    t = ms.gold_images(GOLD)
    assert len(t) == 23
    assert sorted(set(t["system"]), key=float) == ["1", "2", "3", "4", "5", "9", "13", "110"]
    assert np.all(t["z"] > 0.9)
    assert t.meta["provenance"] == "observed"


def test_gates_pass():
    g = ms.gates(_maps(0.01, 0.02), _ip(["1", "1", "2", "2"], [0.2, 0.3, 0.1, 0.4]))
    assert g["G1"]["passes"] and g["G2"]["passes"] and g["pinned_for_use"]


def test_gates_unsolved_system_fails_g2():
    g = ms.gates(_maps(0.01, 0.02), _ip(["1", "1", "9", "9"], [0.2, 0.3, np.nan, np.nan]))
    assert g["G1"]["passes"]
    assert not g["G2"]["passes"]
    assert g["G2"]["systems_reproduced"] == {"1": True, "9": False}
    assert not g["pinned_for_use"]


def test_gates_far_match_and_rms_and_maps():
    # one image beyond the 3" match radius fails G2 even with a small rms elsewhere
    assert not ms.gates(_maps(0.0, 0.0), _ip(["1"] * 3, [0.1, 0.1, 3.2]))["G2"]["passes"]
    # rms above 1" fails G2 although every image is within 3"
    assert not ms.gates(_maps(0.0, 0.0), _ip(["1"] * 2, [1.2, 1.2]))["G2"]["passes"]
    # a magnification map difference above 5 % fails G1
    g = ms.gates(_maps(0.0, 0.06), _ip(["1"] * 2, [0.1, 0.1]))
    assert not g["G1"]["passes"] and not g["pinned_for_use"]


def test_model_registered():
    spec = lc.MODELS["macs1149-sharon"]
    assert lc.is_map_model("macs1149-sharon")
    assert set(spec["files"]) == {"alpha_x", "alpha_y", "kappa_map", "mag_z2"}
    assert spec["z_lens"] == lc.MODELS["macs1149-cats"]["z_lens"]
