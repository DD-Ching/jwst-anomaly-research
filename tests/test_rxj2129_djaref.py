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
