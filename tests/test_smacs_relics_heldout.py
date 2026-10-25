import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import mars_hff_heldout as mo  # noqa: E402
import smacs_relics_heldout as sr  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def test_inputs_are_the_frozen_ones():
    assert set(sr.MODELS) == {"glafic_v2", "zitrin-ltm-gauss_v2", "lenstool_v1"}
    seen = Table.read(ROOT / sr.SEEN)
    assert len(seen) == 17 and {"image_id", "ra", "dec"} <= set(seen.colnames)
    fam = {str(i).split(".")[0] for i in seen["image_id"]}
    assert fam == {"1", "2", "3", "4", "5"}  # the five HST-era families of Golubchik+22


def _mahler_like(seen: Table, dra: float, ddec: float) -> Table:
    """The seen rows moved by (dra, ddec) arcsec, plus one far image of a new family."""
    cos = np.cos(np.deg2rad(np.asarray(seen["dec"], float)))
    t = Table(
        {
            "image_id": [str(i).rstrip("ab") for i in seen["image_id"]] + ["7.1"],
            "ra": list(np.asarray(seen["ra"], float) + dra / 3600 / cos) + [110.79],
            "dec": list(np.asarray(seen["dec"], float) + ddec / 3600) + [-73.449],
            "z": [1.5] * len(seen) + [5.17],
        }
    )
    return t


def test_frame_offset_and_seen_split_on_the_golubchik_list():
    seen = Table.read(ROOT / sr.SEEN)
    t = _mahler_like(seen, 0.3, -0.2)
    off = mo.frame_offset(t, seen)
    assert off["applied"]
    s = mo.split_seen(mo.shift(t, off), seen)
    assert s["seen"][:-1].all() and not s["seen"][-1]  # the new family is held out
