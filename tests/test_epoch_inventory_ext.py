import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import epoch_inventory_ext as ex  # noqa: E402


def _poly(ra0, dec0, half_deg):
    pts = [
        (ra0 - half_deg, dec0 - half_deg),
        (ra0 + half_deg, dec0 - half_deg),
        (ra0 + half_deg, dec0 + half_deg),
        (ra0 - half_deg, dec0 + half_deg),
    ]
    return "POLYGON ICRS " + " ".join(f"{a} {b}" for a, b in pts)


def test_classify_thresholds():
    assert ex.classify({"F200W": {"a_ov": 3.0, "a_core": 1.01}}) == "core"
    assert ex.classify({"F200W": {"a_ov": 3.0, "a_core": 1.0}}) == "overlap"  # strict >
    assert ex.classify({"F200W": {"a_ov": 0.0, "a_core": 0.0}}) == "none"
    assert ex.classify({}) == "none"
    # "some shared filter" is an or over filters
    per = {"F150W": {"a_ov": 0.5, "a_core": 0.2}, "F444W": {"a_ov": 4.0, "a_core": 2.9}}
    assert ex.classify(per) == "core"


def test_pair_areas_core_and_offset():
    centre = (10.0, 0.0)
    rows = [
        # epochs A and B: 2' x 2' boxes on the centre -> core disc fully covered (pi arcmin^2)
        {
            "obs_id": "jw00001-o001_t001_nircam_clear-f200w",
            "filters": "F200W",
            "s_region": _poly(10.0, 0.0, 1 / 60),
        },
        {
            "obs_id": "jw00002-o001_t001_nircam_clear-f200w",
            "filters": "F200W",
            "s_region": _poly(10.0, 0.0, 1 / 60),
        },
        # epoch C: 1' box 3' east -> no core
        {
            "obs_id": "jw00003-o001_t001_nircam_clear-f200w",
            "filters": "F200W",
            "s_region": _poly(10.0 + 3 / 60, 0.0, 0.5 / 60),
        },
    ]
    pairs = [
        {"early": "00001-001", "late": "00002-001", "shared_filters": ["F200W", "F444W"]},
        {"early": "00001-001", "late": "00003-001", "shared_filters": ["F200W"]},
    ]
    out, epoch_core = ex.pair_areas(pairs, rows, centre)
    assert out[0]["class"] == "core"
    assert abs(out[0]["per_filter"]["F200W"]["a_core"] - 3.1416) < 0.05
    assert out[0]["missing_filters"] == ["F444W"]
    assert out[1]["class"] == "none"
    assert epoch_core["00003-001:F200W"] == 0.0


def test_used_excludes_inventory_dirs(tmp_path):
    (tmp_path / "grid").mkdir()
    (tmp_path / "grid" / "RESULT.md").write_text("pair jw01176-o211 -> jw01176-o212")
    for d in ex.INVENTORY_DIRS:
        (tmp_path / d).mkdir()
        (tmp_path / d / "inventory.json").write_text('{"e": "jw09999-o001"}')
    (tmp_path / "notes.md").write_text("jw02282-o010")
    assert ex.used_excluding_inventories(tmp_path) == {
        ("01176", "211"),
        ("01176", "212"),
        ("02282", "010"),
    }
