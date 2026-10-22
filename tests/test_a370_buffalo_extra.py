"""Tests of the Abell 370 BUFFALO extra-image script: constants and the input pinning (PLAN.md)."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import a370_buffalo_extra as m  # noqa: E402
import m0417_extra as mx  # noqa: E402
import m0417_images as mi  # noqa: E402


def test_constants_name_the_planned_field_model_and_tile():
    assert (m.FIELD, m.MODEL, m.OBS, m.TILE, m.GRID_STEP_K) == (
        "abell370",
        "abell370-buffalo-a-maps",
        "o002",
        "t001",
        0.25,
    )
    assert m.IMAGES == "buffalo_sl_gold.ecsv"
    assert m.OUT == Path("results/lens_consistency/abell370_buffalo_extra")


def test_refuses_to_run_until_the_canucs_files_are_pinned(tmp_path, monkeypatch):
    monkeypatch.setattr(m, "CAT_SHA", None)
    with pytest.raises(SystemExit, match="pin CAT_SHA"):
        m.check_inputs(tmp_path / "c", tmp_path / "i")


def test_refuses_files_with_other_content(tmp_path, monkeypatch):
    f = tmp_path / "f"
    f.write_bytes(b"x")
    monkeypatch.setattr(m, "CAT_SHA", "0" * 64)
    monkeypatch.setattr(m, "IMAGE_SHA", "0" * 64)
    with pytest.raises(SystemExit, match="not the CANUCS"):
        m.check_inputs(f, f)


def test_register_points_the_reused_stages_at_this_field(monkeypatch):
    monkeypatch.setattr(mi, "FIELDS", dict(mi.FIELDS))
    monkeypatch.setattr(mx, "FIELDS", dict(mx.FIELDS))
    monkeypatch.setattr(mx, "IMAGES", mx.IMAGES)
    monkeypatch.setattr(mx, "OUT", mx.OUT)
    m.register()
    assert mi.FIELDS["abell370"] == {"images": m.IMAGES, "models": (m.MODEL,)}
    assert mx.FIELDS["abell370"] == (m.IMAGES, m.MODEL, "o002", "t001")
    assert (mx.IMAGES, mx.OUT) == (m.OUT, m.FRAME)
