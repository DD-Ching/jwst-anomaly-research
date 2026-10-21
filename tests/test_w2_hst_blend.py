"""Offline tests of the line-C HST rules (scripts/w2_hst_blend_chain.py); synthetic data only."""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.modeling import models
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w2_hst_blend as cli  # noqa: E402
import w2_hst_blend_chain as ch  # noqa: E402
import w12_niq_hst as hst  # noqa: E402

P = ch.Params()
PIX = 0.04


def _classes():
    rows = [
        ("L1", "blended", True),
        ("L2", "too close", False),
        ("L3", "blended", False),
        ("L4", "too close", True),
        ("L5", "none", False),  # decided by LS: never selected
        ("L6", "blended", False),
    ]
    return Table(
        {
            "system_id": [r[0] for r in rows],
            "name": [f"N{r[0]}" for r in rows],
            "ra": np.arange(len(rows), dtype=float),
            "dec": np.zeros(len(rows)),
            "selection": ["quasar"] * len(rows),
            "status": [r[1] for r in rows],
            "ls_sep": np.full(len(rows), np.nan),
            "theta_e": np.full(len(rows), np.nan),
            "z_lens": np.full(len(rows), np.nan),
            "z_source": np.full(len(rows), np.nan),
            "lens_z_known": [r[2] for r in rows],
        }
    )


def _coverage():
    # (system, instrument, filter, exptime, programme, contains)
    rows = [
        ("L1", "ACS/WFC", "F814W", 2000.0, "100", True),
        ("L2", "WFC3/UVIS", "F814W", 1500.0, "100", True),
        ("L3", "ACS/WFC", "F814W", 3000.0, "200", True),
        ("L4", "ACS/WFC", "F814W", 5000.0, "300", True),  # control from another programme
        ("L5", "ACS/WFC", "F814W", 9000.0, "100", True),
        ("L6", "ACS/WFC", "F814W", 900.0, "100", True),  # too shallow
        ("L6", "WFC3/IR", "F814W", 5000.0, "100", True),  # wrong instrument
        ("L2", "ACS/WFC", "F814W", 9000.0, "100", False),  # does not contain the position
    ]
    return Table(
        {
            "system_id": [r[0] for r in rows],
            "obs_id": [f"j{i:03d}aa{i:02d}0" for i in range(len(rows))],  # one visit per row
            "instrument": [r[1] for r in rows],
            "filters": [r[2] for r in rows],
            "t_exptime": [r[3] for r in rows],
            "proposal_id": [r[4] for r in rows],
            "contains": [r[5] for r in rows],
        }
    )


def test_select_systems_rule():
    s = ch.select_systems(_classes(), _coverage(), P)
    tg = s[s["role"] == "target"]
    ct = s[s["role"] == "control"]
    assert list(tg["system_id"]) == ["L3", "L2"]  # deepest first; L6 too shallow; L5 decided
    assert list(tg["t814_s"]) == [3000.0, 1500.0]  # non-containing exposure not counted
    assert list(ct["system_id"]) == ["L1"]  # L4 shares no programme with a target


def test_select_systems_caps_and_ties():
    p = ch.Params(max_targets=1)
    c = _classes()
    cov = _coverage()
    cov["t_exptime"][1] = 3000.0  # L2 ties with L3: system_id order decides
    s = ch.select_systems(c, cov, p)
    assert list(s[s["role"] == "target"]["system_id"]) == ["L2"]


def test_classify_cases():
    c = ch.classify
    assert c(1, np.nan, np.nan, np.nan, 1.0, P).startswith("unresolved")
    assert c(2, 0.5, 10, np.nan, 1.0, P) == "too compact for the core mask"
    assert c(2, 1.5, 10, 2.2, 1.0, P) == "pair mismatch"
    assert c(4, 1.5, 10, 2.2, 1.0, P) == "light between images"  # quads are not pair-checked
    assert c(2, 1.5, 3.0, np.nan, 1.0, P) == "not detected between images"
    assert c(2, 1.5, -6.0, np.nan, 1.0, P).startswith("fit inconclusive")
    assert c(2, 1.5, 6.0, np.nan, np.nan, P) == "noise not estimable"
    assert c(2, 1.5, 5.0, 1.6, 1.0, P) == "light between images"


def test_control_gate():
    light, none = ch.DECIDED
    nine = [light] * 9 + [none]
    assert ch.control_gate(nine, P)["passed"]
    assert not ch.control_gate([light] * 8 + [none] * 2, P)["passed"]
    assert not ch.control_gate([light] * 4, P)["passed"]  # too few decided controls
    g = ch.control_gate([light] * 5 + ["unresolved (< 2 point sources)"] * 5, P)
    assert g["passed"] and g["n_decided"] == 5 and g["efficiency"] == 1.0


def _scene(lens_mag=None, reff=0.3, seed=1, quad=False):
    rng = np.random.default_rng(seed)
    n = 120
    yy, xx = np.indices((n, n))
    zp = 25.9
    sigma = 0.01
    img = rng.normal(0.0, sigma, (n, n))
    pts = [(45.0, 60.0, 50.0), (75.0, 60.0, 25.0)]
    if quad:
        pts = [(45.0, 60.0, 50.0), (75.0, 60.0, 30.0), (60.0, 45.0, 20.0), (60.0, 76.0, 15.0)]
    for x, y, a in pts:
        img += models.Moffat2D(a, x, y, gamma=1.6, alpha=2.8)(xx, yy)
    if lens_mag is not None:
        lx, ly = (60.0, 60.0)
        psf = ch.moffat_kernel(1.6, 2.8)
        g = ch.sersic_unit(img.shape, lx, ly, reff / PIX, 4, 0.2, psf=psf)
        img += 10 ** (-0.4 * (lens_mag - zp)) * g
    return img, sigma, zp, pts


def _noise_fn(sigma):
    return lambda npix: sigma


def test_find_point_sources_skips_extended_lens():
    img, sigma, _, pts = _scene(lens_mag=20.0, reff=0.4)
    pk = ch.find_point_sources(img, (60, 60), PIX, sigma, P)
    assert len(pk) == 2
    for x, y, _ in pts:
        assert min(np.hypot(x - a, y - b) for a, b in pk) <= 1.5


@pytest.mark.parametrize("quad", [False, True])
def test_measure_detects_lens_and_not_empty(quad):
    valid = np.ones((120, 120), bool)
    img, sigma, zp, _ = _scene(lens_mag=None, quad=quad)
    base = ch.measure(img, valid, (60, 60), PIX, P, sigma, _noise_fn(sigma))
    assert base["n_img"] == (4 if quad else 2)
    assert base["status"] == "not detected between images"
    img2, _, _, _ = _scene(lens_mag=21.0, quad=quad)
    m = ch.measure(img2, valid, (60, 60), PIX, P, sigma, _noise_fn(sigma))
    assert m["status"] == "light between images"
    assert ch.recovered(base, m, P)


def test_remove_halos_n_equals_d064_for_a_pair():
    rng = np.random.default_rng(3)
    res = rng.normal(0, 1, (60, 60))
    pos = [(20.3, 30.1), (41.7, 29.4)]
    fit = SimpleNamespace(
        x_0_1=SimpleNamespace(value=pos[0][0]), y_0_1=SimpleNamespace(value=pos[0][1]),
        x_0_2=SimpleNamespace(value=pos[1][0]), y_0_2=SimpleNamespace(value=pos[1][1]),
    )  # fmt: skip
    assert np.allclose(ch.remove_halos_n(res, pos), hst.remove_halos(res, fit))


def test_lens_position_sis_and_centroid():
    x, y = ch.lens_position([(0.0, 0.0), (10.0, 0.0)], [2.0, 1.0])
    assert x == pytest.approx(10 - 10 * 0.5 / 1.5) and y == 0.0
    assert ch.lens_position([(0, 0), (2, 0), (0, 2), (2, 2)], [1, 1, 1, 1]) == (1.0, 1.0)


def test_point_in_region():
    reg = "POLYGON 10.0 -1.0 10.1 -1.0 10.1 -0.9 10.0 -0.9"
    assert cli.point_in_region(10.05, -0.95, reg)
    assert not cli.point_in_region(10.15, -0.95, reg)
    wrap = "POLYGON 359.95 0.0 0.05 0.0 0.05 0.1 359.95 0.1"
    assert cli.point_in_region(0.01, 0.05, wrap)


def test_predicted_f814w_monotonic():
    a = ch.predicted_f814w(0.5, 0.5, 2.0, P)["m_typ"]
    b = ch.predicted_f814w(1.0, 0.5, 2.0, P)["m_typ"]
    assert b < a  # a larger Einstein radius needs a more luminous lens
    assert ch.predicted_f814w(1.0, 0.5, 2.0, P)["m_faint"] == pytest.approx(b + 2 * 0.88)


def test_run_helpers():
    import w2_hst_blend_run as run

    t = Table({"a": np.ma.MaskedArray([1.0, 0.0], mask=[False, True])})
    assert np.isnan(run.fcol(t, "a")[1]) and run.fcol(t, "a")[0] == 1.0  # masked is not 0
    yy, xx = np.indices((801, 801))
    num = models.Moffat2D(2.0, 400, 400, gamma=3.0, alpha=3.0)(xx, yy).sum()
    assert run.moffat_flux(2.0, 3.0, 3.0) == pytest.approx(num, rel=1e-3)


def test_visit_exposure_counts_each_exposure_once():
    # SDSSJ0246-0825 rows in hst_coverage.ecsv: HAP-MVM skycells, HAP-SVM visits and pipeline
    # association products of the same exposures (v1 summed them: 67232 s)
    rows = [
        ("hst_skycell-p1064x06y18_wfc3_uvis_f814w_all", "WFC3/UVIS", 19360.0),
        ("hst_skycell-p1064x06y18_acs_wfc_f814w_all", "ACS/WFC", 2288.0),
        ("hst_skycell-p1153x06y01_acs_wfc_f814w_all", "ACS/WFC", 2288.0),
        ("hst_skycell-p1153x06y01_wfc3_uvis_f814w_all", "WFC3/UVIS", 19360.0),
        ("j8oi08040", "ACS/WFC", 1144.0),
        ("hst_14254_53_wfc3_uvis_f814w_icud53", "WFC3/UVIS", 9680.0),
        ("hst_9744_08_acs_wfc_f814w_j8oi08", "ACS/WFC", 2288.0),
        ("j8oi08030", "ACS/WFC", 1144.0),
        ("hst_14254_02_wfc3_uvis_f814w_icud02", "WFC3/UVIS", 9680.0),
    ]
    o, i, t = zip(*rows, strict=True)
    assert ch.visit_exposure(o, i, t) == 2288.0 + 2 * 9680.0
    # a visit with pipeline products only is summed over its associations
    assert (
        ch.visit_exposure(["idxo03020", "idxo03030"], ["WFC3/UVIS"] * 2, [700.0, 728.0]) == 1428.0
    )
    # an association and a separate exposure of the same visit: the association only
    o = ["jf5626010", "jf5626ryq"]
    assert ch.visit_exposure(o, ["ACS/WFC"] * 2, [674.0, 5.0]) == 674.0
    # individual exposures only (no association, no SVM): summed
    assert ch.visit_exposure(["jf5626ryq", "jf5626rzq"], ["ACS/WFC"] * 2, [300.0, 300.0]) == 600.0
    # multi-visit mosaics alone give no depth
    assert ch.visit_exposure(["hst_skycell-p1_acs_wfc_f814w_all"], ["ACS/WFC"], [674.0]) == 0.0


def test_merge_rows_keeps_other_systems(tmp_path):
    import w2_hst_blend_run as run

    f = tmp_path / "m.ecsv"
    prev = [
        {"system_id": "A", "role": "control", "status": "old"},
        {"system_id": "B", "role": "control", "status": "keep"},
        {"system_id": "C", "role": "control", "status": "deselected"},
        {"system_id": "T", "role": "target", "status": "other role"},
    ]
    Table(rows=prev).write(f)
    t = run.merge_rows(f, [{"system_id": "A", "role": "control", "status": "new"}], "control",
                       {"A", "B"})  # fmt: skip
    got = {(str(r["system_id"]), str(r["status"])) for r in t}
    assert got == {("A", "new"), ("B", "keep"), ("T", "other role")}


def test_role_summary_empty_role():
    import w2_hst_blend_run as run

    meas = Table(
        rows=[{"system_id": "A", "role": "control", "status": "light between images",
               "m_lim_reff0.3": 22.0, "max_sep_arcsec": 1.5}]
    )  # fmt: skip
    s = run.role_summary(meas, "target")  # the gate-failed path: no target rows
    assert s["target_status_counts"] == {} and s["target_max_sep_range"] is None
    assert s["target_m_lim_reff0.3"]["n"] == 0
    assert run.role_summary(meas, "control")["control_m_lim_reff0.3"]["median"] == 22.0


def test_hapcut_failure_is_not_no_data(monkeypatch, tmp_path):
    import warnings

    import w2_hst_blend_run as run
    from astroquery.mast import Hapcut

    monkeypatch.setattr(run.time, "sleep", lambda s: None)
    c = run.SkyCoord(10.0, 0.0, unit="deg")

    def boom(*a, **k):
        raise ConnectionError("proxy reset")

    monkeypatch.setattr(Hapcut, "download_cutouts", boom)
    t, st, msg = run.hapcut_with_retry(c, 250, tmp_path)
    assert t is None and st == run.FAILED and "proxy reset" in msg

    def empty(*a, **k):
        from astroquery.exceptions import NoResultsWarning

        warnings.warn("Missing HAP files for input target.", NoResultsWarning, stacklevel=1)
        return Table(names=("Local Path",), dtype=(str,))

    monkeypatch.setattr(Hapcut, "download_cutouts", empty)
    t, st, msg = run.hapcut_with_retry(c, 250, tmp_path)
    assert st == run.NO_HAP and "Missing HAP files" in msg


def test_fp_reference_uses_max_system():
    import w2_hst_blend_run as run

    t = Table(
        {"system_id": ["A", "A", "B", "B"], "psf": ["epsf", "moffat", "epsf", "epsf"],
         "inject_mag": [np.nan, np.nan, np.nan, 20.0], "resid_snr": [3.0, 50.0, 9.0, 99.0]}
    )  # fmt: skip
    assert run.fp_reference(t) == "B"


def test_method_validation_owner_rule():
    ok_gate = {"passed": True, "efficiency": 0.9}
    v = ch.method_validation(ok_gate, {"epsf_no_lens_max_snr": 13.0}, P)
    assert not v["validated"] and v["target_interpretation"].startswith("not interpreted")
    assert ch.method_validation(ok_gate, {"epsf_no_lens_max_snr": 4.9}, P)["validated"]
    assert ch.method_validation(ok_gate, {"epsf_no_lens_max_snr": None}, P)["verdict"].startswith(
        "incomplete"
    )
    assert (
        ch.method_validation({"passed": False}, {"epsf_no_lens_max_snr": 1.0}, P)["verdict"]
        == "failed (control gate)"
    )


def test_no_lens_max_all_nan_and_json_safe():
    import json

    import w2_hst_blend_run as run

    t = Table({"system_id": ["A", "B"], "psf": ["epsf", "epsf"], "inject_mag": [np.nan, np.nan],
               "resid_snr": [np.nan, np.nan]})  # fmt: skip
    assert run.no_lens_max(t) == (None, None)
    assert run.fp_reference(t) is None
    assert json.dumps({"x": run.finite_or_none(np.nan)}, allow_nan=False) == '{"x": null}'
    t["resid_snr"][1] = 13.0
    assert run.no_lens_max(t) == (13.0, "B")


def test_reuse_check_and_label_normalisation():
    import w2_hst_blend_run as run

    row = {"hst_image": "a.fits", "img_pos_px": "1.00,2.00;3.00,4.00", "img_amps": "1;2",
           "psf_gamma": 1.5, "psf_alpha": 2.0, "sigma_pix": 0.005, "zp_ab": 25.1}  # fmt: skip
    assert run.reuse_ok(row, dict(row))
    assert not run.reuse_ok(row, dict(row, img_pos_px="1.00,2.01;3.00,4.00"))
    assert not run.reuse_ok(row, dict(row, sigma_pix=0.0051))
    assert not run.reuse_ok(row, dict(row, hst_image="b.fits"))
    assert run.normalise_status("none between images") == "not detected between images"
    assert run.normalise_status("light between images") == "light between images"


def test_apply_verdict_idempotent_and_restores():
    import w2_hst_blend_run as run

    light = "light between images"
    meas = Table(
        rows=[
            {"system_id": "C", "role": "control", "status": light},
            {"system_id": "T", "role": "target", "status": light},
        ]
    )
    fail = {
        "validated": False,
        "target_interpretation": "not interpreted (method validation failed)",
    }
    ok = {"validated": True, "target_interpretation": "interpreted under PLAN.md"}
    a = run.apply_verdict(meas, fail)
    b = run.apply_verdict(a, fail)
    assert (
        list(b["status"])
        == list(a["status"])
        == ["light between images", fail["target_interpretation"]]
    )
    c = run.apply_verdict(b, ok)
    assert list(c["status"]) == ["light between images", "light between images"]
    # a row added later by measure (no stored chain status) is filled from its status
    from astropy.table import vstack

    d = vstack(
        [
            a,
            Table(
                rows=[{"system_id": "U", "role": "target", "status": "not detected between images"}]
            ),
        ]
    )
    e = run.apply_verdict(d, ok)
    assert list(e[run.RAW])[-1] == "not detected between images"


def test_download_status_helpers(tmp_path):
    import w2_hst_blend_run as run

    st = {"A": {"status": run.NO_HAP}, "B": {"status": run.FAILED}, "C": {"status": "in progress"}}
    assert run.missing_cutout_status(st, "A") == "no HAP cutout"
    assert run.missing_cutout_status(st, "B") == run.FAILED
    assert run.missing_cutout_status(st, "C") == run.NOT_ATTEMPTED  # interrupted attempt
    assert run.missing_cutout_status(st, "Z") == run.NOT_ATTEMPTED  # never attempted
    d = tmp_path / "C"
    d.mkdir()
    (d / "x.fits").write_bytes(b"")
    assert not run.cutout_complete(st, "C", d)  # partial set from an interrupted run
    assert run.cutout_complete({"C": {"status": "ok"}}, "C", d)


def test_merge_rows_no_inject_drops_stale_rows(tmp_path):
    import w2_hst_blend_run as run

    f = tmp_path / "inj.ecsv"
    Table(rows=[{"system_id": "A", "role": "control", "inject_mag": 20.0},
                {"system_id": "B", "role": "control", "inject_mag": 20.0}]).write(f)  # fmt: skip
    t = run.merge_rows(f, [], "control", {"A", "B"}, {"A"})  # A re-measured with --no-inject
    assert list(t["system_id"]) == ["B"]


def test_needs_download_no_hap_final(tmp_path):
    import w2_hst_blend_run as run

    d = tmp_path / "A"
    assert run.needs_download({}, "A", d)  # never attempted
    assert not run.needs_download({"A": {"status": run.NO_HAP}}, "A", d)  # final answer
    assert run.needs_download({"A": {"status": run.NO_HAP}}, "A", d, force=True)
    assert run.needs_download({"A": {"status": run.FAILED}}, "A", d)
    d.mkdir()
    (d / "x.fits").write_bytes(b"")
    assert run.needs_download({"A": {"status": "in progress"}}, "A", d)  # partial set
    assert not run.needs_download({"A": {"status": "ok"}}, "A", d)


def test_run_system_refuses_partial_cutout(tmp_path, monkeypatch):
    import json

    import w2_hst_blend_run as run

    monkeypatch.setattr(run, "RES", tmp_path)
    (tmp_path / run.DOWNLOAD_STATUS).write_text(json.dumps({"S1": {"status": "in progress"}}))
    d = tmp_path / "cut" / "S1"
    d.mkdir(parents=True)
    (d / "partial.fits").write_bytes(b"")  # would fail to load if it were used
    r = {"ra": 10.0, "dec": 0.0, "system_id": "S1", "name": "N", "role": "target",
         "status_ls": "blended", "t814_s": 1428.0}  # fmt: skip
    row, panel, inj = run.run_system(r, tmp_path / "cut", P)
    assert row["status"] == run.NOT_ATTEMPTED and panel is None and inj == []


def test_psf_one_redoes_partial_wide_dir(tmp_path, monkeypatch):
    import json

    import w2_hst_blend_run as run
    import w12_niq_epsf as ep

    monkeypatch.setattr(run, "RES", tmp_path)
    (tmp_path / run.DOWNLOAD_STATUS).write_text(json.dumps({"wide_S1": {"status": run.FAILED}}))
    d = tmp_path / "wide_S1"
    d.mkdir()
    (d / "stale.fits").write_bytes(b"")
    monkeypatch.setattr(run, "hapcut_with_retry", lambda c, size, path: (None, run.NO_HAP, "none"))
    r = {"ra": 10.0, "dec": 0.0, "name": "N"}
    rows, ok = run.psf_one("S1", r, {}, d, tmp_path / "log", P, ep.Params())
    assert not ok and rows[0]["status"] == "no wide F814W cutout"
    assert not (d / "stale.fits").exists()  # the partial set was removed, not reused
    st = json.loads((tmp_path / run.DOWNLOAD_STATUS).read_text())
    assert st["wide_S1"]["status"] == run.NO_HAP


def test_no_lens_max_ignores_inf():
    import w2_hst_blend_run as run

    t = Table({"system_id": ["A", "B"], "psf": ["epsf", "epsf"], "inject_mag": [np.nan, np.nan],
               "resid_snr": [np.inf, 13.0]})  # fmt: skip
    assert run.no_lens_max(t) == (13.0, "B")
    assert run.finite_or_none(np.inf) is None


def test_pick_reuse_refusal_and_missing_column():
    import w2_hst_blend_run as run

    cur = {"hst_image": "a.fits", "img_pos_px": "1,2;3,4", "img_amps": "1;2", "psf_gamma": 1.5,
           "psf_alpha": 2.0, "sigma_pix": 0.005, "zp_ab": 25.1}  # fmt: skip
    psf = Table({"system_id": ["S"], "psf": ["epsf"], "status": ["none between images"]})
    bad = Table(rows=[dict(cur, system_id="S", sigma_pix=0.006)])
    good = Table(rows=[dict(cur, system_id="S")])
    rows, mism = run.pick_reuse("S", cur, [(psf, bad), (psf, good)])
    assert rows is not None and not mism  # a later source accepted: no refusal
    assert rows[0]["status"] == "not detected between images"
    rows, mism = run.pick_reuse("S", cur, [(psf, bad)])
    assert rows is None and mism
    nocol = Table(rows=[{k: v for k, v in dict(cur, system_id="S").items() if k != "zp_ab"}])
    assert not run.reuse_ok(nocol[0], cur)  # missing input column: mismatch, recompute
