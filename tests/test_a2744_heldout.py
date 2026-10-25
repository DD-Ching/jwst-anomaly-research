import sys
from pathlib import Path

from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import a2744_heldout as ah  # noqa: E402
import mars_hff_inputs as mh  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def test_inputs_are_the_frozen_ones():
    t = Table.read(ROOT / ah.IMAGES)
    assert len(t) > 100 and bool(t["z_spec"].all())
    assert (t["z"] > 0).all()
    assert ah.MODEL == "abell2744-bergamini23"
    assert ah.MARS_FIELD in mh.FIELDS


def test_shared_over_3():
    blocks = {"a": {"over_3arcsec": ["1.1a", "2.1b"]}, "b": {"over_3arcsec": ["2.1b", "9a"]}}
    assert ah.shared_over_3(blocks) == ["2.1b"]
    assert ah.shared_over_3({"a": {"over_3arcsec": []}}) == []


def test_posthoc_split(tmp_path):
    imgs = Table(
        {
            "image_id": ["80.1a", "70a", "70b"],
            "ra": [3.5, 3.5, 3.5],
            "dec": [-30.376, -30.40, -30.41],
        }
    )
    imgs.write(tmp_path / "i.ecsv", format="ascii.ecsv")
    r = Table(
        {
            "image_id": ["80.1a", "70a", "70b"],
            "system": ["80.1", "70", "70"],
            "z": [3.67, 2.39, 2.39],
            "dtheta_arcsec": [4.0, 0.3, 0.4],
            "inside": [True, True, True],
        }
    )
    r.write(tmp_path / f"heldout_{ah.MODEL}.ecsv", format="ascii.ecsv")
    p = ah.posthoc(tmp_path, tmp_path / "i.ecsv")
    assert p["north"] == {"n": 1, "rms_arcsec": 4.0, "within_3arcsec": 0.0, "systems": ["80"]}
    assert p["north_z_ge_3"]["n"] == 1 and p["north_z_lt_3"]["n"] == 0
    assert p["core"]["n"] == 2 and p["core"]["within_3arcsec"] == 1.0
    assert p["without_system_80"]["rms_arcsec"] == 0.354


def test_posthoc_counts_no_prediction_as_a_miss(tmp_path):
    imgs = Table({"image_id": ["72a", "72b"], "ra": [3.55, 3.55], "dec": [-30.375, -30.3751]})
    imgs.write(tmp_path / "i.ecsv", format="ascii.ecsv")
    r = Table(
        {
            "image_id": ["72a", "72b"],
            "system": ["72", "72"],
            "z": [3.75, 3.75],
            "dtheta_arcsec": [1.0, float("inf")],
            "inside": [True, True],
        }
    )
    r.write(tmp_path / f"heldout_{ah.MODEL}.ecsv", format="ascii.ecsv")
    p = ah.posthoc(tmp_path, tmp_path / "i.ecsv")
    assert p["north_z_ge_3"] == {
        "n": 2,
        "rms_arcsec": 1.0,
        "within_3arcsec": 0.5,
        "systems": ["72"],
    }
    assert p["min_sep_north_z_ge_3_to_z_lt_3_arcsec"] is None


def test_mars_vs_bergamini(tmp_path):
    Table({"image_id": ["86a", "62.1a"], "dtheta_arcsec": [1.08, 0.5]}).write(
        tmp_path / f"heldout_{ah.MODEL}.ecsv", format="ascii.ecsv"
    )
    Table(
        {
            "image_id": ["86a", "83a", "62.1a"],
            "dtheta_arcsec": [4.86, 51.7, 0.8],
            "inside": [True, True, False],
        }
    ).write(tmp_path / "heldout_A2744.ecsv", format="ascii.ecsv")
    assert ah.mars_vs_bergamini(tmp_path) == [
        {"image_id": "86a", "bergamini_arcsec": 1.08, "mars_arcsec": 4.86}
    ]
