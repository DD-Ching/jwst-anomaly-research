"""Tests for scripts/w12_lenscats.py helpers (tiny synthetic tables, offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest
from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("w12_lenscats", _DIR / "w12_lenscats.py")
w12 = importlib.util.module_from_spec(_spec)
sys.modules["w12_lenscats"] = w12  # dataclasses look their module up
_spec.loader.exec_module(w12)

D = 1 / 3600


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
    assert w12.selection_class(_row("SPT0418-47")) == "submm"
    assert w12.selection_class(_row("J0904", refs="Negrello et al. (2016)")) == "submm"
    assert w12.selection_class(_row("B1555+375", refs="lensedquasars")) == "radio"
    assert w12.selection_class(_row("MG0751+2716")) == "radio"
    lemon = "https://research.ast.cam.ac.uk/lensedquasars/"
    assert w12.selection_class(_row("J0011-0845", refs=lemon)) == "quasar"
    assert w12.selection_class(_row("020000+100000", system_type="GQ")) == "quasar"
    assert w12.selection_class(_row("DESI-049.7700-49.3639", refs="Storfer et al.")) == "galaxy"
    # case-sensitive survey tokens: "Others" is not HerS, "class" is not CLASS
    assert w12.selection_class(_row("X", refs="Others et al.; a class of lenses")) == "galaxy"


def test_simbad_galaxy():
    assert w12.is_simbad_galaxy("[ASK2018] DES J0407-5006 G", "QSO")  # lens-galaxy component
    assert w12.is_simbad_galaxy("SDSS J1", "Galaxy")
    assert not w12.is_simbad_galaxy("QSO J0158-4325", "QSO")


def test_box_terms_wrap_and_pole():
    assert len(w12.box_terms(180.0, 10.0, 5.0)) == 1
    lo = w12.box_terms(0.0005, 0.0, 5.0)
    assert len(lo) == 2 and "AND 360.0000000" in lo[1]
    hi = w12.box_terms(359.9995, 0.0, 5.0)
    assert len(hi) == 2 and "ra BETWEEN 0.0000000" in hi[1]
    assert w12.box_terms(10.0, 89.999, 5.0) == ["(dec BETWEEN 89.9976111 AND 90.0000000)"]


def test_poisson95():
    assert w12.poisson95(0) == pytest.approx(2.996, abs=1e-3)
    assert w12.poisson95(1) == pytest.approx(4.744, abs=1e-3)


def _sources(rows):
    return Table(rows=rows, names=("ra", "dec", "type", "mag_z", "maskbits"))


def _system(sel):
    return Table(
        {"system_id": ["q"], "ra": [10.0], "dec": [0.0], "selection": [sel], "defl_mag_max": [21.0]}
    )


IMAGES = [(10.0 - 1.2 * D, 0.0, "PSF", 19.0, 0), (10.0 + 1.2 * D, 0.0, "PSF", 19.5, 0)]
LENS = [(10.0 + 0.1 * D, 0.0, "DEV", 20.0, 0)]


def test_quasar_test_can_fail_and_injection_recovers():
    p = w12.Params()
    sysq = _system("quasar")
    with_lens = w12.deflector_test(sysq, _sources(IMAGES + LENS), p)
    assert with_lens["test_status"][0] == "deflector"
    assert with_lens["image_sep"][0] == pytest.approx(2.4, abs=0.01)
    assert w12.deflector_test(sysq, _sources(IMAGES), p)["test_status"][0] == "none"
    faint = _sources(IMAGES + [(10.0, 0.0, "DEV", 22.5, 0)])
    assert w12.deflector_test(sysq, faint, p)["test_status"][0] == "none"  # too faint
    close = [(10.0 - 0.6 * D, 0.0, "PSF", 19.0, 0), (10.0 + 0.6 * D, 0.0, "PSF", 19.5, 0)]
    assert w12.deflector_test(sysq, _sources(close), p)["test_status"][0] == "too close"
    assert w12.deflector_test(sysq, _sources(IMAGES[:1]), p)["test_status"][0] == "blended"
    # galaxy-selected systems are insensitive whatever the imaging shows
    sysg = _system("galaxy")
    assert w12.deflector_test(sysg, _sources(IMAGES), p)["test_status"][0] == "insensitive"
    res = sysq.copy()
    res["test_status"] = with_lens["test_status"]
    inj = w12.inject_dark(res, _sources(IMAGES + LENS), p)
    assert inj["n"] == 1 and inj["recovery"] == 1.0


def test_radio_test():
    p = w12.Params()
    sysr = _system("radio")
    gal = [(10.0 + 0.3 * D, 0.0, "DEV", 20.0, 0)]
    assert w12.deflector_test(sysr, _sources(gal), p)["test_status"][0] == "deflector"
    far = [(10.0 + 4.0 * D, 0.0, "DEV", 20.0, 0)]
    assert w12.deflector_test(sysr, _sources(far), p)["test_status"][0] == "none"
