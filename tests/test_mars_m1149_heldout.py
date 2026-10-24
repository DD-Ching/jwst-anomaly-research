import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import mars_hff_heldout as mo  # noqa: E402
import mars_hff_inputs as mh  # noqa: E402
import mars_m1149_heldout as m1  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def _list() -> Table:
    return Table.read(ROOT / m1.IMAGES)


def test_schuldt_list_matches_the_published_table():
    t = _list()
    assert len(t) == 18 and t["z_spec"].all() and set(t["qf"]) <= {3, 9}
    fams = [s for s in np.unique([mo.lensmodel.image_family(i) for i in t["image_id"]])]
    assert sorted(fams) == ["201", "202", "203", "204", "301", "401", "SXa"]
    assert np.all((t["z"] > 1.4) & (t["z"] < 6.0))
    assert m1.FIELD in mh.FIELDS


def test_all_held_out_gives_type_n_anchors_and_no_anchor_for_sxa():
    t = _list()
    far = Table({"ra": [0.0], "dec": [0.0]})  # no MARS row near: every image is held out
    s = mo.split_seen(t, far)
    assert not s["seen"].any()
    a = mo.anchors(s)
    assert "SXa" not in a and set(a) == {"201", "202", "203", "204", "301", "401"}
    assert {k for k, _ in a.values()} == {"N"}
    assert sum(int(m.sum()) for _, m in a.values()) == 17


def test_one_seen_image_makes_type_m():
    t = _list()
    i = int(np.where(t["image_id"] == "301.2")[0][0])
    mars = Table({"ra": [t["ra"][i]], "dec": [t["dec"][i]]})
    a = mo.anchors(mo.split_seen(t, mars))
    assert a["301"][0] == "M" and a["201"][0] == "N"
    assert not mo.frame_offset(t, mars)["applied"]  # 1 pair < 5: O = 0
