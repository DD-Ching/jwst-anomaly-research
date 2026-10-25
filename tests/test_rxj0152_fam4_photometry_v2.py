"""Frozen rule functions of scripts/rxj0152_fam4_photometry_v2.py (PLAN.md rules 6-10)."""

import importlib.util
from pathlib import Path

import pytest
from astropy.coordinates import SkyCoord

_PATH = Path(__file__).resolve().parents[1] / "scripts" / "rxj0152_fam4_photometry_v2.py"
_spec = importlib.util.spec_from_file_location("rxj0152_fam4_photometry_v2", _PATH)
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)

REF_SED = {b: float(i + 1) for i, b in enumerate(m.BANDS)}
SIG = {b: 0.01 for b in m.BANDS}


def test_flux_window():
    lo, hi = m.flux_window(100.0)
    assert lo == pytest.approx(100.0 / 7.482 / 3)
    assert hi == pytest.approx(3 * 100.0 / 5.325)


def test_colour_match_scaled_copy_and_offsets():
    assert m.colour_match({b: v / 6 for b, v in REF_SED.items()}, REF_SED, SIG)
    off = {b: v / 6 for b, v in REF_SED.items()}
    off["f090w"] *= 10**0.31
    assert not m.colour_match(off, REF_SED, SIG)
    off["f090w"] = REF_SED["f090w"] / 6 * 10**0.29
    assert m.colour_match(off, REF_SED, SIG)
    neg = {b: v / 6 for b, v in REF_SED.items()}
    neg["f444w"] = -1.0
    assert not m.colour_match(neg, REF_SED, SIG)


def test_colour_match_skips_undetectable_bands_and_needs_three():
    sig = {b: 1e6 for b in m.BANDS}  # nothing detectable
    f = {b: v / 6 for b, v in REF_SED.items()}
    assert not m.colour_match(f, REF_SED, sig)
    sig = dict(sig, f090w=0.01, f115w=0.01, f150w=0.01)
    bad = dict(f, f444w=-5.0)  # undetectable band ignored
    assert m.colour_match(bad, REF_SED, sig)
    assert not m.colour_match(dict(f, f277w=0.0), REF_SED, SIG)


def test_chance_p():
    lam, p = m.chance_p(62, 61.5, 2501.3, 3)
    assert lam == pytest.approx(62 * 61.5 / (2501.3 - 61.5))
    assert 0 < p < 1
    assert m.chance_p(0, 61.5, 2501.3, 1)[1] == 0.0
    assert m.chance_p(10, 61.5, 2501.3, 0)[1] == 1.0


def test_v_pass():
    assert m.v_pass([True] * 3, [True] * 24 + [False] * 6)
    assert not m.v_pass([True] * 3, [True] * 23 + [False] * 7)
    assert not m.v_pass([True, False, True], [True] * 30)


def test_outcomes():
    assert m.outcome(False, True, True, 0, 1.0) == "O3_undetermined"
    assert m.outcome(True, False, True, 2, 0.01) == "O3_undetermined"
    assert m.outcome(True, True, True, 1, 0.05) == "O1_counterpart"
    assert m.outcome(True, True, False, 1, 0.2) == "O1c_chance_level"
    assert m.outcome(True, True, True, 0, 1.0) == "O2_missing"
    assert m.outcome(True, True, False, 0, 1.0) == "O3_undetermined"
    assert m.known_case_pass(6) and not m.known_case_pass(5)


def test_offsets_and_region():
    c = m.PREDICTED["glafic"]
    offs = m.offset_positions(c)
    assert len(offs) == 9
    for p in offs:
        assert SkyCoord(*p, unit="deg").separation(
            SkyCoord(*c, unit="deg")
        ).arcsec == pytest.approx(1.5)
        assert m.in_region(p)
    assert not m.in_region((28.18344, -13.95682))  # 4.1 is far from every prediction
