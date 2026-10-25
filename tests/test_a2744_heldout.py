import sys
from pathlib import Path

from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import a2744_heldout as ah  # noqa: E402
import mars_hff_inputs as mh  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def test_inputs_are_the_frozen_ones():
    t = Table.read(ROOT / ah.IMAGES)
    assert len(t) > 100 and bool(t["z_spec"].all())
    assert (t["z"] > 0).all()
    assert ah.MODEL == "abell2744-bergamini23"
    assert ah.MARS_FIELD in mh.FIELDS


def test_knots_are_separate_families():
    assert lensmodel.image_family("1.1a") == "1.1"
    assert lensmodel.image_family("1.2a") == "1.2"
    assert lensmodel.image_family("83b") == "83"


def test_shared_over_3():
    blocks = {"a": {"over_3arcsec": ["1.1a", "2.1b"]}, "b": {"over_3arcsec": ["2.1b", "9a"]}}
    assert ah.shared_over_3(blocks) == ["2.1b"]
    assert ah.shared_over_3({"a": {"over_3arcsec": []}}) == []
