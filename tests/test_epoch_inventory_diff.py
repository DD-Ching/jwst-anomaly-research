import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from epoch_inventory_diff import diff_inventories  # noqa: E402


def _inv(epochs, centre=None):
    f = {"epochs": epochs}
    if centre:
        f["centre_deg"] = centre
    return {"meta": {}, "fields": {"x": f}}


def test_added_removed_newly_public_and_centre():
    old = _inv(
        {"1-1": {"public": True}, "2-1": {"public": False}, "3-1": {"public": True}}, [10.0, 0.0]
    )
    new = _inv(
        {"1-1": {"public": True}, "2-1": {"public": True}, "4-1": {"public": False}},
        [10.0, 1 / 3600],
    )
    d = diff_inventories(old, new)["x"]
    assert d["added"] == ["4-1"]
    assert d["removed"] == ["3-1"]
    assert d["newly_public"] == ["2-1"]
    assert d["centre_shift_arcsec"] == 1.0


def test_unchanged_field_is_empty():
    inv = _inv({"1-1": {"public": True}}, [1.0, 1.0])
    assert diff_inventories(inv, inv)["x"] == {"added": [], "removed": [], "newly_public": []}


def test_centre_shift_wraps_at_ra_zero():
    old = _inv({}, [359.99995, 0.0])
    new = _inv({}, [0.00001, 0.0])
    assert "centre_shift_arcsec" not in diff_inventories(old, new)["x"]
