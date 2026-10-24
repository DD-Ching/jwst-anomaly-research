"""Tests of the Abell S1063 extra-image script: constants, rule 1 threshold, rule 2 repeat check
and the model hook (PLAN.md)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bullet_extra as be  # noqa: E402
import lens_consistency as lc  # noqa: E402
import m0417_extra as mx  # noqa: E402
import s1063_extra as m  # noqa: E402

STORED = {"n_images": 36, "rms_arcsec": 0.58, "matched_fraction": 1.0, "passes": True}


def test_constants_name_the_planned_field_model_tile_and_reference():
    assert (m.FIELD, m.MODEL, m.OBS, m.TILE) == ("abells1063", "abells1063-wslap", "o001", "t001")
    assert m.S3_3293.format("o001", "t001", "i2d.fits").endswith(
        "jw03293/L3/t/o001/jw03293-o001_t001_nircam_clear-f200w_i2d.fits"
    )
    assert m.DJA_PHOT.endswith("/v7/abells1063-grizli-v7.0-fix_phot.fits")
    assert (m.F_MIN_N, m.GRID_STEP_K, m.K_REPEAT_TOL_ARCSEC) == (18, 0.25, 1e-6)
    assert (be.F_MATCH_ARCSEC, be.F_MAX_MEDIAN_ARCSEC) == (0.5, 0.1)


def test_spec_images_are_the_36_gate_k_images():
    t = m.spec_images_table(m.sw.spec_images())
    assert len(t) == 36 and "system" not in t.colnames
    assert len({str(i).split(".")[0] for i in t["image_id"]}) == 12


def test_k_repeat_needs_the_stored_numbers_and_a_pass():
    assert m.k_repeat(dict(STORED), STORED)["passes"]
    assert not m.k_repeat(dict(STORED, rms_arcsec=0.581), STORED)["reproduced"]
    bad = dict(STORED, passes=False)
    r = m.k_repeat(bad, dict(bad))
    assert r["reproduced"] and not r["passes"]


def test_frame_check_f_uses_the_lower_count():
    ra0, dec0 = 342.18, -44.53
    k = np.arange(20)
    ra, dec = ra0 + 0.001 * (k % 5), dec0 + 0.001 * (k // 5)
    assert be.frame_check_f(ra, dec, ra, dec, min_n=m.F_MIN_N)["passes"]
    assert not be.frame_check_f(ra[:17], dec[:17], ra, dec, min_n=m.F_MIN_N)["passes"]
    shifted = be.frame_check_f(ra, dec + 0.15 / 3600, ra, dec, min_n=m.F_MIN_N)
    assert shifted["n_matched"] == 20 and not shifted["passes"]


def test_register_points_the_frame_tables_at_go3293_tile():
    m.register()
    assert mx.FIELDS[m.FIELD] == (m.IMAGES_SPEC, m.MODEL, "o001", "t001")
    assert mx.IMAGES == m.OUT and mx.OUT == m.FRAME and mx.S3 == m.S3_3293
    mx.FIELDS.pop(m.FIELD)


def test_run_field_uses_a_given_model_without_registry_or_offset(monkeypatch):
    seen = {}

    def no_lookup(name):
        raise AssertionError("registry lookup with a given model")

    class Stop(Exception):
        pass

    def backtrace(model, images, _):
        seen["model"] = model
        raise Stop

    monkeypatch.setattr(lc, "load_model", no_lookup)
    monkeypatch.setattr(lc, "apply_frame_offset", no_lookup)
    monkeypatch.setattr(mx, "load_images", lambda field: None)
    monkeypatch.setattr(mx.lensmodel, "backtrace_images", backtrace)
    monkeypatch.setitem(mx.FIELDS, "x", ("f", "name", "o", "t"))
    sentinel = object()
    try:
        mx.run_field("x", Path("."), model=sentinel)
    except Stop:
        pass
    assert seen["model"] is sentinel
