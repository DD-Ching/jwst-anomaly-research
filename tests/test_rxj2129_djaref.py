"""Offline checks of the RXJ2129 image-list transcription (PLAN_1.md inputs)."""

import hashlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import rxj2129_djaref as rx  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

# sha256 of "paper_id,ra,dec;..." after the row-by-row check against the TeX table (2026-10-23)
COORD_SHA256 = "a77d5206a0a2d6e695e929842a5a71f5b2d9e32e743959369408ff1e86ec7d60"


def test_image_table_counts_and_families():
    t = rx.image_table()
    assert len(t) == 31
    assert len(set(t["system"])) == 10
    assert all(t["z_spec"])
    # ids parse back to their family, so back-traces group correctly
    assert [lensmodel.image_family(i) for i in t["image_id"]] == list(t["system"])
    assert len(set(zip(t["ra"], t["dec"], strict=True))) == 31  # no duplicate positions


def test_image_table_redshifts():
    t = rx.image_table()
    z = dict(zip(t["system"], t["z"], strict=False))
    assert z["G"] == 9.51
    assert z["S3-1"] == z["S3-2"] == z["3"] == 1.5194
    assert np.all((t["dec"] > 0.083) & (t["dec"] < 0.095))
    assert t.meta["provenance"] == "observed"


def test_image_coordinates_pinned():
    """The transcribed positions, checked row by row against the arXiv TeX table on 2026-10-23."""
    t = rx.image_table()
    text = ";".join(
        f"{p},{r:.7f},{d:.7f}" for p, r, d in zip(t["paper_id"], t["ra"], t["dec"], strict=True)
    )
    assert hashlib.sha256(text.encode()).hexdigest() == COORD_SHA256


def _merge(systems, offsets_arcsec, classes, search=1.0):
    ra0, dec0 = 322.4165, 0.0905
    ra = [ra0 + dx / 3600.0 for dx, _ in offsets_arcsec]
    dec = [dec0 + dy / 3600.0 for _, dy in offsets_arcsec]
    return rx.merge_knot_targets(systems, ra, dec, [search] * len(systems), classes)


def test_knot_group_joins_overlapping_targets_of_different_families():
    # S3-1 absent overlaps a recovered system-3 target: one region, detected -> no flag
    out = _merge(["S3-1", "3"], [(0.0, 0.0), (1.5, 0.0)], ["absent", "recovered"])
    assert out["flags"] == []
    assert out["knot_sets"][0]["class"] == "recovered"
    assert out["knot_sets"][0]["systems"] == ["S3-1", "3"]


def test_knot_group_all_absent_flags_once():
    out = _merge(["S3-1", "S3-2", "3"], [(0, 0), (1.5, 0), (3.0, 0)], ["absent"] * 3)
    assert len(out["flags"]) == 1 and out["flags"][0].startswith("S3-1@")  # transitive join
    assert out["knot_sets"][0]["rows"] == [0, 1, 2]


def test_knot_group_order_ambiguous_before_absent():
    out = _merge(["3", "S3-2"], [(0, 0), (0.5, 0)], ["absent", "ambiguous"])
    assert out["flags"] == [] and out["knot_sets"][0]["class"] == "ambiguous"


def test_knot_group_leaves_other_targets_alone():
    # outside the group, same family, or circles not overlapping: per-row rule 6
    out = _merge(
        ["4", "3", "S3-1", "S3-1", "2"],
        [(0, 0), (0.5, 0), (10, 0), (10.5, 0), (20, 0)],
        ["absent", "recovered", "absent", "absent", "recovered"],
    )
    assert out["knot_sets"] == []
    assert len(out["flags"]) == 3
    # separation just beyond the sum of the radii does not join
    out = _merge(["3", "S3-1"], [(0, 0), (2.01, 0)], ["absent", "recovered"])
    assert out["knot_sets"] == [] and len(out["flags"]) == 1
