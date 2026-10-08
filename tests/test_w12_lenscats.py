"""Tests for scripts/w12_lenscats.py helpers (tiny synthetic tables, offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest
from astropy.table import Table

from jwst_anomaly import lenscats

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
    """rows: (ra, dec, type, mag_g, mag_z, maskbits)."""
    return Table(rows=rows, names=("ra", "dec", "type", "mag_g", "mag_z", "maskbits"))


def _system(sel):
    return Table(
        {"system_id": ["q"], "ra": [10.0], "dec": [0.0], "selection": [sel], "defl_mag_max": [21.0]}
    )


IMAGES = [(10.0 - 1.2 * D, 0.0, "PSF", 19.3, 19.0, 0), (10.0 + 1.2 * D, 0.0, "PSF", 19.8, 19.5, 0)]
LENS = [(10.0 + 0.1 * D, 0.0, "DEV", 22.0, 20.0, 0)]


def _test(sel, rows, p=None):
    p = p or w12.Params()
    t, src = _system(sel), _sources(rows)
    images = lenscats.pair_images(t, src, p.image_radius)
    return w12.deflector_test(t, src, p, images)["test_status"][0]


def test_quasar_test_can_fail():
    assert _test("quasar", IMAGES + LENS) == "deflector"
    assert _test("quasar", IMAGES) == "none"  # a dark lens can fail the test
    assert _test("quasar", IMAGES + [(10.0, 0.0, "DEV", 24.0, 22.5, 0)]) == "faint galaxy"
    close = [
        (10.0 - 0.6 * D, 0.0, "PSF", 19.3, 19.0, 0),
        (10.0 + 0.6 * D, 0.0, "PSF", 19.8, 19.5, 0),
    ]
    assert _test("quasar", close) == "too close"
    assert _test("quasar", IMAGES[:1]) == "blended"
    assert _test("galaxy", IMAGES) == "insensitive"


def test_quad_images_are_not_the_deflector():
    # 3rd and 4th images: PSF, colour like the pair, inside the pair circle
    quad = IMAGES + [
        (10.0, 0.9 * D, "PSF", 20.3, 20.0, 0),
        (10.0, -0.9 * D, "PSF", 20.4, 20.1, 0),
    ]
    assert _test("quasar", quad) == "none"
    # a red compact source typed PSF is a lens candidate, not an image
    red = IMAGES + [(10.0, 0.2 * D, "PSF", 22.5, 20.0, 0)]
    assert _test("quasar", red) == "deflector"


def test_pair_circle_beyond_image_radius():
    # wide pair (4.8'') centred 1.5'' north of the catalogue position: a lens 3.7'' from the
    # position (outside image_radius + 0.5) but inside the pair circle is still found
    p = w12.Params()
    imgs = [
        (10.0 - 2.4 * D, 1.5 * D, "PSF", 19.3, 19.0, 0),
        (10.0 + 2.4 * D, 1.5 * D, "PSF", 19.8, 19.5, 0),
    ]
    lens = [(10.0, 3.7 * D, "DEV", 22.0, 20.0, 0)]
    assert _test("quasar", imgs + lens, p) == "deflector"


def test_radio_test():
    gal = [(10.0 + 0.3 * D, 0.0, "DEV", 22.0, 20.0, 0)]
    assert _test("radio", gal) == "deflector"
    assert _test("radio", [(10.0 + 0.3 * D, 0.0, "DEV", 24.0, 22.5, 0)]) == "faint galaxy"
    assert _test("radio", [(10.0 + 4.0 * D, 0.0, "DEV", 22.0, 20.0, 0)]) == "none"


def test_empty_inputs_keep_schema():
    t = _system("quasar")[:0]
    src = _sources(IMAGES)
    img = lenscats.pair_images(t, src)
    assert list(img.colnames) == list(lenscats.PAIR_COLUMNS) and len(img) == 0
    out = lenscats.quasar_pair_test(t, src, [], img)
    assert list(out.colnames) == list(lenscats.PAIR_TEST_COLUMNS) and len(out) == 0
    res = w12.deflector_test(t, src, w12.Params(), img)
    assert "test_status" in res.colnames and len(res) == 0
    one = _system("quasar")
    none_src = _sources(IMAGES)[:0]
    assert lenscats.pair_images(one, none_src)["n_images"][0] == 0
