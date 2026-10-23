import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w2_desi_evidence as ev  # noqa: E402


def _side(gaia_id, d_gaia, d_tr, g):
    return {
        "gaia_nearest": {"source_id": gaia_id, "d": d_gaia},
        "tractor_nearest": {"d": d_tr, "type": "PSF", "mag_g": g, "mag_r": "", "mag_z": ""},
    }


def test_flags_distinct_sources_and_close_pair():
    rec = {
        "a": _side("1", 0.05, 0.05, "19.0"),
        "b": _side("2", 0.10, 0.10, "19.5"),
        "tractor": [{}, {}],
        "sep_desi": 1.07,
        "sep_cat": 1.02,
    }
    f = ev.flags(rec)
    assert f["both_match_distinct_gaia"] and f["close_le_2arcsec"]
    assert f["dg_tractor"] == 0.5 and f["sep_desi_minus_cat"] == 0.05


def test_flags_one_fibre_without_its_own_source():
    # J2351-like: fibre a is 1.94" from the only bright Gaia source, which is fibre b's match
    rec = {
        "a": _side("7", 1.94, 0.07, "23.45"),
        "b": _side("7", 0.11, 0.11, "19.57"),
        "tractor": [{}, {}, {}],
        "sep_desi": 2.05,
        "sep_cat": 1.86,
    }
    f = ev.flags(rec)
    assert f["gaia_a"] is None and f["gaia_b"] == "7"
    assert not f["both_match_distinct_gaia"]
    assert f["dg_tractor"] == 3.88 and not f["close_le_2arcsec"]


def test_flags_ignore_non_finite_magnitudes():
    rec = {
        "a": _side("1", 0.05, 0.05, "NaN"),
        "b": _side("2", 0.10, 0.10, "19.5"),
        "tractor": [{}, {}],
        "sep_desi": 2.5,
        "sep_cat": 2.5,
    }
    assert ev.flags(rec)["dg_tractor"] is None
