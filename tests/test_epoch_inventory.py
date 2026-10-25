import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import epoch_inventory as ei  # noqa: E402


def _row(obs_id, filt, t_min, release=0.0):
    return {
        "obs_id": obs_id,
        "filters": filt,
        "t_min": t_min,
        "t_obs_release": release,
        "dataRights": "PUBLIC",
    }


def test_epoch_key_and_filter():
    assert ei.epoch_key("jw02282-o010_t001_nircam_clear-f200w") == ("02282", "010")
    assert ei.epoch_key("hst_123") is None
    assert ei.filter_name("CLEAR;F200W") == "F200W"
    assert ei.filter_name("F405N;F444W") == "F405N"


def test_new_pair_rule():
    rows = [
        _row("jw00001-o001_t001_nircam_clear-f150w", "F150W", 100.0),
        _row("jw00001-o001_t001_nircam_clear-f444w", "F444W", 100.0),
        _row("jw00002-o001_t001_nircam_clear-f150w", "F150W", 120.0),  # dt 20 d: fails (c)
        _row(
            "jw00003-o001_t001_nircam_clear-f200w", "F200W", 300.0
        ),  # no shared filter with 1 / 2: fails (b)
        _row("jw00004-o001_t001_nircam_clear-f444w", "F444W", 400.0),
        _row(
            "jw00005-o001_t001_nircam_clear-f444w", "F444W", 500.0, release=1e9
        ),  # not public: fails (a)
    ]
    ep = ei.build_epochs(rows)
    assert ep[("00001", "001")]["filters"] == {"F150W", "F444W"}
    pairs = ei.new_pairs(ep, used=set(), now_mjd=1000.0)
    assert [(p["early"], p["late"]) for p in pairs] == [("00001-001", "00004-001")]
    # (d): both epochs already used -> not new
    assert ei.new_pairs(ep, used={("00001", "001"), ("00004", "001")}, now_mjd=1000.0) == []
    assert len(ei.new_pairs(ep, used={("00001", "001")}, now_mjd=1000.0)) == 1


def test_used_epochs(tmp_path):
    (tmp_path / "a.md").write_text("pair jw01176-o211 vs jw01176213 and jw02282-o010_t001")
    assert ei.used_epochs(tmp_path) == {("01176", "211"), ("01176", "213"), ("02282", "010")}
