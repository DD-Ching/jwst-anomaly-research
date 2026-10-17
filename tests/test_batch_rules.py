"""The batch pass rule and class rule (scripts/batch_run.py, scripts/batch_review.py), offline."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import batch_pairs as bp  # noqa: E402
import batch_review as br  # noqa: E402
import batch_run as bru  # noqa: E402


def test_box_passes_boundaries():
    assert bru.box_passes(7, 10, 2, 0.70, 2)  # exactly at both limits
    assert not bru.box_passes(6, 10, 0, 0.70, 2)  # vanish fraction below
    assert not bru.box_passes(10, 10, 3, 0.70, 2)  # too many false hits
    assert not bru.box_passes(0, 0, 0, 0.70, 2)  # no vanish sources: undeterminable


def test_classify():
    assert br.classify(True, {"F150W": -30, "F444W": -10}, None) == "artifact (edge)"
    assert br.classify(False, {"F150W": -2.9, "F444W": 2.0}, None) == "no significant change"
    assert br.classify(False, {"F150W": -27, "F444W": -13}, None) == "real change"
    assert br.classify(False, {"F150W": -5, "F444W": 4}, None) == "undetermined"
    assert br.classify(False, {"F150W": -5}, None).startswith("undetermined (photometry")
    assert br.classify(False, {"F150W": -27, "F444W": -13}, "veto: spike").startswith("artifact")
    assert (
        br.classify(False, {"F150W": 1, "F444W": 1}, "looked at: fine") == "no significant change"
    )


def test_wilson():
    assert bru.wilson(0, 0) is None
    lo, hi = bru.wilson(62, 71)
    assert lo < 62 / 71 < hi


def test_touching_boxes_do_not_overlap():
    import math

    ra, dec = 64.04, -24.07
    step = bp.BOX / 3600 / math.cos(math.radians(dec))
    neighbour = (round(ra + step, 6), round(dec, 6))  # stored rounded, as in boxes.json
    assert not bp.overlaps(ra, dec, [neighbour])
    assert not bp.overlaps(ra - step, dec, [(round(ra, 6), round(dec, 6))])
    assert bp.overlaps(ra + step / 2, dec, [neighbour])
    assert bp.overlaps(ra, dec, [(ra, dec)])
