import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import heldout_synthesis as hs  # noqa: E402


def _t(rows):
    return Table(
        rows=rows,
        names=["image_id", "system", "type", "z", "dtheta_arcsec", "n_predicted", "inside"],
    )


def test_offsets_drops_outside_and_marks_no_prediction_inf():
    t = _t(
        [
            ("1a", "1", "N", 2.0, 0.5, 3, True),
            ("1b", "1", "N", 2.0, np.nan, 0, True),
            ("2a", "2", "N", 2.0, np.nan, 0, False),
        ]
    )
    off = hs.offsets(t)
    assert set(off) == {"1a", "1b"}
    assert off["1a"] == ("1", 0.5) and np.isinf(off["1b"][1])


def test_model_stats_counts_inf_in_percentiles_not_rms():
    off = {f"{i}a": (str(i), d) for i, d in enumerate([0.1, 0.2, 0.3, np.inf])}
    s = hs.model_stats(off)
    assert s["n"] == 4 and s["n_no_prediction"] == 1
    assert np.isclose(s["rms_arcsec"], np.sqrt((0.01 + 0.04 + 0.09) / 3))
    assert s["median_arcsec"] == 0.2 and np.isinf(s["p95_arcsec"]) and s["frac_over_3"] == 0.25


def test_isolated_needs_two_good_family_members():
    off = {
        "1a": ("1", 4.0),
        "1b": ("1", 0.5),
        "1c": ("1", 0.9),
        "2a": ("2", 4.0),
        "2b": ("2", 0.5),
        "3a": ("3", 4.0),
        "3b": ("3", 0.5),
        "3c": ("3", 1.2),
        "4a": ("4", np.inf),
        "4b": ("4", 0.1),
        "4c": ("4", 0.2),
    }
    assert hs.isolated_misses(off) == ["1a", "4a"]


def test_shared_and_flags():
    per = {"A": {"1a": ("1", 4.0), "2a": ("2", 3.5)}, "B": {"1a": ("1", np.inf), "2a": ("2", 1.0)}}
    sh = hs.shared_misses(per)
    assert sh == {"1a": ["A", "B"]}
    assert hs.flags(sh, {"A": [], "B": ["1a"]}) == ["1a"]
    assert hs.flags(sh, {"A": ["2a"], "B": []}) == []


def test_field_floor_takes_lowest_rms_model():
    stats = {
        "A": {"n": 5, "rms_arcsec": 2.0, "p95_arcsec": 4.0},
        "B": {"n": 5, "rms_arcsec": 1.0, "p95_arcsec": 2.5},
    }
    assert hs.field_floor(stats) == {"model": "B", "p95_arcsec": 2.5, "rms_arcsec": 1.0}


def test_stored_tables_exist():
    for _, _, rel, _ in hs.TABLES:
        assert (Path(__file__).resolve().parents[1] / hs.ROOT / rel).exists(), rel


def test_seen_inclusive_isolated_uses_seen_members():
    off = {"6c": ("6", 3.05), "7c": ("7", 4.0)}
    seen = {"6a": 0.03, "6b": 0.05, "7a": 0.02}
    fam = lambda i: i.rstrip("abc")  # noqa: E731
    assert hs.seen_inclusive_isolated(off, seen, fam) == ["6c"]
