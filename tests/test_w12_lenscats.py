"""Tests for scripts/w12_lenscats.py helpers (tiny synthetic tables, offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("w12_lenscats", _DIR / "w12_lenscats.py")
w12 = importlib.util.module_from_spec(_spec)
sys.modules["w12_lenscats"] = w12  # dataclasses look their module up
_spec.loader.exec_module(w12)


def _row(name, entries="", refs="", system_type=""):
    t = Table(
        {
            "name": [name],
            "entries": [entries],
            "refs": [refs],
            "system_type": [system_type],
            "simbad": [""],
        }
    )
    return t[0]


def test_selection_class():
    assert w12.selection_class(_row("SPT0418-47")) == "submm/radio"
    assert w12.selection_class(_row("J0904", refs="Negrello et al. (2016)")) == "submm/radio"
    assert (
        w12.selection_class(
            _row("J0011-0845", refs="https://research.ast.cam.ac.uk/lensedquasars/")
        )
        == "quasar"
    )
    assert w12.selection_class(_row("020000+100000", system_type="GQ")) == "quasar"
    assert w12.selection_class(_row("DESI-049.7700-49.3639", refs="Storfer et al.")) == "galaxy"
    # case-sensitive survey tokens: "Others" is not HerS, "class" is not CLASS
    assert w12.selection_class(_row("X", refs="Others et al.; a class of lenses")) == "galaxy"
    assert w12.is_submm_radio(_row("HATLASJ121301.5-004922"))


def test_simbad_galaxy():
    assert w12.is_simbad_galaxy("[ASK2018] DES J0407-5006 G", "QSO")  # lens-galaxy component
    assert w12.is_simbad_galaxy("SDSS J1", "Galaxy")
    assert not w12.is_simbad_galaxy("QSO J0158-4325", "QSO")


def test_limits(tmp_path):
    n = 10
    systems = Table(
        {
            "system_id": [f"L{i}" for i in range(n)],
            "covered": [True] * 9 + [False],
            "detectable_typical": [True] * 8 + [False, True],
            "detectable": [True, True] + [False] * 8,
            "selection": ["quasar", "galaxy"] + ["galaxy"] * 8,
        }
    )
    systems.write(tmp_path / "systems.ecsv")
    vetted = Table(
        {
            "system_id": ["L2", "L3"],
            "verdict": [
                "point sources within 3'': lens blended or PSF-typed",  # undecided
                "galaxy bright enough for the lens within 5''",  # decisive
            ],
        }
    )
    lim = w12.limits(tmp_path, vetted, w12.Params())
    assert lim["unexplained"] == 0 and lim["covered"] == 9 and lim["decisive"] == 8
    assert lim["typical_all_N"] == 7  # L0, L1, L3..L7
    assert lim["typical_all_f95"] == pytest.approx(2.996 / 7)
    assert lim["typical_lens_light_independent_N"] == 1
    assert lim["conservative_all_N"] == 2
    vetted["verdict"][1] = "unexplained by catalogue tests"
    lim = w12.limits(tmp_path, vetted, w12.Params())
    assert lim["unexplained"] == 1 and lim["typical_all_f95"] is None
    assert np.isfinite(lim["typical_all_N"])
