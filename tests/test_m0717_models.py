"""Offline tests for scripts/m0717_models.py (registered models, transcribed image list; no network)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import lens_consistency as lc  # noqa: E402
import m0717_models as mm  # noqa: E402
import m1149_sharon as ms  # noqa: E402

GOLD = ROOT / "results/lens_consistency/macs0717_models/sharon_gold_images.ecsv"


def test_gold_images_table():
    t = ms.gold_images(GOLD)
    assert len(t) == 26
    assert sorted(set(t["system"]), key=float) == ["1", "3", "4", "6", "12", "13", "14", "19"]
    assert np.all(t["z"] > 1.7)
    assert t.meta["provenance"] == "observed"
    # one redshift per system
    for s in set(t["system"]):
        assert len(set(t["z"][t["system"] == s])) == 1


def test_models_registered():
    for name in mm.GATED:
        spec = lc.MODELS[name]
        assert lc.is_map_model(name)
        assert set(spec["files"]) == {"alpha_x", "alpha_y", "kappa_map", "mag_z2"}
        assert spec["z_lens"] == lc.MODELS[mm.COMPARISON]["z_lens"]
    shas = [sha for name in mm.GATED for _, sha in lc.MODELS[name]["files"].values()]
    assert len(set(shas)) == len(shas) == 8


def test_rms_ignores_unsolved():
    r, n = mm.rms(Table({"dtheta_arcsec": np.array([3.0, 4.0, np.nan])}))
    assert n == 2 and np.isclose(r, np.sqrt(12.5))
    assert mm.rms(Table({"dtheta_arcsec": np.array([np.nan])})) == (None, 0)
