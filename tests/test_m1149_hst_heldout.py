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
