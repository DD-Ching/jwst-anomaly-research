import sys
from pathlib import Path

from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import m1149_hst_heldout as hh  # noqa: E402
import mars_hff_heldout as mo  # noqa: E402
import mars_m1149_heldout as m1  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def test_inputs_are_the_frozen_ones():
    assert hh.MODELS == ("macs1149-cats", "macs1149-sharon")
    assert len(Table.read(ROOT / m1.IMAGES)) == 18
    gold = Table.read(ROOT / hh.SHARON_GOLD)
    assert len(gold) == 23 and {"ra", "dec", "image_id"} <= set(gold.colnames)


def test_sharon_seen_list_is_the_gold_list():
    assert len(hh.seen_list("macs1149-sharon", {}, None)) == 23


def test_schuldt_images_are_not_in_the_sharon_gold_list():
    t = Table.read(ROOT / m1.IMAGES)
    gold = Table.read(ROOT / hh.SHARON_GOLD)
    s = mo.split_seen(t, gold)
    assert not s["seen"].any()  # held out by construction: no Sharon gold row within 0.5"
    assert not mo.frame_offset(t, gold)["applied"]


def test_posthoc_groups_rms_and_max(tmp_path):
    rows = Table(
        {
            "image_id": ["401.1", "401.2", "202.1"],
            "system": ["401", "401", "202"],
            "dtheta_arcsec": [3.0, 4.0, float("nan")],
        }
    )
    p = tmp_path / "t.ecsv"
    rows.write(p, format="ascii.ecsv")
    g = hh.posthoc({"m": p})["groups"]
    assert g["401"]["m"] == {"n": 2, "rms_arcsec": 3.536, "max_arcsec": 4.0}
    assert g["202-204"]["m"] == {"n": 0}
