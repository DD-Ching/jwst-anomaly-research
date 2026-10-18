"""Tests for scripts/inject_shear.py (synthetic inputs only)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
sys.path.insert(0, str(Path(__file__).resolve().parent))
_spec = importlib.util.spec_from_file_location("inject_shear", _DIR / "inject_shear.py")
ish = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ish)

from test_inject_radial import SCALE, _model, _shapes  # noqa: E402


def _field(n=900, seed=0):
    rng = np.random.default_rng(seed)
    model = _model()
    x, y = rng.uniform(-40, 40, (2, n))
    shapes = _shapes(
        model,
        x,
        y,
        pa_obs=rng.uniform(0, 180, n),
        ellipticity=rng.uniform(0.1, 0.5, n),
        semimajor_px=rng.uniform(3.0, 6.0, n),
    )
    args = SimpleNamespace(
        **(
            ish.SHEAR_DEFAULTS
            | {"max_radius": 40.0, "recover_tol": 2.0, "n_random": 40, "responsivity": 1.0}
        )
    )
    return model, shapes, args


def test_shear_defaults_come_from_the_screen():
    assert ish.SHEAR_DEFAULTS == ish.es.shear_defaults()
    assert (
        ish.SHEAR_DEFAULTS["filter"] == "schirmer" and ish.SHEAR_DEFAULTS["aperture_arcsec"] == 10
    )


def test_shear_field_spec_overrides_only_the_shear_inputs():
    spec = ish.field_spec("abells1063")
    assert spec["photoz"][0] == "eazy" and spec["spike_stars"] and len(spec["centre"]) == 2
    assert spec["catalog"] == ish.ir.FIELDS["abells1063"]["catalog"]
    assert ish.ir.FIELDS["abells1063"]["photoz"] is None  # the radial entry is unchanged
    assert ish.field_spec("macs0416") == ish.ir.FIELDS["macs0416"]


def test_combine_checks_photoz_against_the_given_specs(capsys):
    eff = {"2e+12": {"efficiency": 0.1, "theta_e_zs2_arcsec": 3.0}}
    summ = {
        "field": "abells1063",
        "photoz": "zout.fits",
        "screened_area_deg2": 1e-3,
        "efficiency": eff,
    }
    out = ish.ir.combine([summ], [2e12], {"abells1063": ish.field_spec("abells1063")})
    assert out["fields"] == ["abells1063"]
    assert "differs from its field spec" not in capsys.readouterr().out
    ish.ir.combine([summ], [2e12])  # against the radial FIELDS entry (no photo-z) it warns
    assert "differs from its field spec" in capsys.readouterr().out


def test_injected_ellipticity_keeps_r_of_the_measured_change_only():
    e_src = np.array([0.1 + 0.05j, np.nan])
    raw_src = np.array([0.3 + 0.0j, np.nan])  # measured, cluster shear included
    raw_img = np.array([0.5 + 0.1j, 0.2j])
    out = ish.injected_ellipticity(e_src, raw_img, raw_src, 0.5)
    # resolved source: the measured change only; an unresolved source is never painted (NaN)
    np.testing.assert_allclose(out[0], 0.1 + 0.05j + 0.5 * (0.2 + 0.1j))
    assert np.isnan(out[1])
    # no lens change: the source's corrected shape comes back unchanged, whatever its cluster g
    same = ish.injected_ellipticity(e_src[:1], raw_src[:1], raw_src[:1], 0.45)
    np.testing.assert_allclose(same, e_src[:1])
    # |e| is capped below 1 (a large inner-image change on an already elliptical source)
    big = ish.injected_ellipticity([0.8 + 0j], [0.95 + 0j], [0.1 + 0j], 0.45)
    np.testing.assert_allclose(big, [0.99 + 0j])


def test_only_measurable_rows_are_lensed():
    model, shapes, args = _field()
    shapes["semimajor_px"][:50] = 0.5  # unresolved after deconvolving a 1 px PSF
    inj = ish.ShearInjector(model, shapes, args, psf_sigma=1.0)
    assert not inj.background[:50].any()
    assert inj.background[50:].sum() > 0.9 * (len(shapes) - 50)


def test_a_user_responsivity_outside_its_range_is_refused():
    model, shapes, _ = _field()
    for bad in (0.0, -0.3, 2.0):
        with pytest.raises(ValueError, match="outside"):
            ish.es.shear_sources(model, shapes, 1.0, 10.0, responsivity=bad)


def test_a_massive_w1_lens_is_recovered_and_no_lens_is_not():
    model, shapes, args = _field()
    inj = ish.ShearInjector(model, shapes, args, psf_sigma=1.0)
    assert inj.r == 1.0  # synthetic shapes carry no cluster shear to measure R from
    inj.draw_null(np.random.default_rng(1))
    big = np.where(inj.background, 6.0, np.nan)  # theta_E 6": dozens of radial images
    hit = inj.trial(15.0, -10.0, big, SCALE)
    assert hit["n_lensed"] > 5 and hit["n_image_used"] > 5
    assert hit["recovered"] and hit["s_best"] > np.max(inj.null)
    none = np.full(len(shapes), np.nan)  # no lensable rows: nothing painted
    miss = inj.trial(15.0, -10.0, none, SCALE)
    assert miss["n_lensed"] == 0 and miss["s_best"] < hit["s_best"]


def test_painted_images_on_a_spike_axis_are_vetoed_like_real_segments():
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy.table import Table

    model, shapes, args = _field()
    star = SkyCoord(shapes["ra"][0], shapes["dec"][0], unit="deg")
    shapes["mag"][0], shapes["is_extended"][0] = 16.0, False  # spikes reach ~19"
    for i, sep in enumerate(np.arange(2.0, 10.0), start=1):  # real segments along PA 0 (north)
        c = star.directional_offset_by(0.0 * u.deg, sep * u.arcsec)
        shapes["ra"][i], shapes["dec"][i], shapes["pa_obs"][i] = c.ra.deg, c.dec.deg, 0.0
        shapes["ellipticity"][i] = 0.6
    inj = ish.ShearInjector(model, shapes, args, psf_sigma=1.0)
    assert inj.base["spike"][1:9].all() and len(inj.spike_ref) >= 8
    # painted: on the north axis and aligned (vetoed); 30 deg off every axis, pointing at the star
    pas = [0.0, 30.0, 0.0]
    pos = [star.directional_offset_by(pa * u.deg, 6.5 * u.arcsec) for pa in (180.0, 30.0)]
    pos.append(star.directional_offset_by(0.0 * u.deg, 40.0 * u.arcsec))  # beyond the spike
    img = Table(
        {
            "ra": [c.ra.deg for c in pos],
            "dec": [c.dec.deg for c in pos],
            "pa_obs": pas,
            "label": [10**6, 10**6 + 1, 10**6 + 2],
        }
    )
    np.testing.assert_array_equal(inj.spike_veto(img), [True, False, False])


def test_setting_reads_missing_options_as_legacy_values():
    old = {"assumptions": {"screen": {"aperture_arcsec": 10.0, "max_radius": 100.0}}}
    new = {
        "assumptions": {
            "screen": {
                "aperture_arcsec": 10.0,
                "max_radius": 60.0,
                "psf_chi": [0.0, 0.0],
                "min_coverage": 0.0,
            },
            "spike_stars": 25,
        }
    }
    assert ish.setting(old) == ish.setting(new)
    cut = {"assumptions": {"screen": {**new["assumptions"]["screen"], "min_coverage": 0.5}}}
    assert ish.setting(old) != ish.setting(cut)


def test_streak_veto_blanks_rows_and_painted_images():
    model, shapes, args = _field()
    ra0 = float(np.median(shapes["ra"]))

    def streak(ra, dec):  # a north-south "trail" through the field centre
        return np.abs(np.asarray(ra) - ra0) * 3600.0 * np.cos(np.radians(-30.0)) < 4.0

    on = streak(shapes["ra"], shapes["dec"])
    assert 20 < on.sum() < len(shapes) // 2
    plain = ish.ShearInjector(model, shapes, args, psf_sigma=1.0)
    inj = ish.ShearInjector(model, shapes, args, psf_sigma=1.0, streak=streak)
    c = inj.base["counts"]
    assert inj.n_streak == int(on.sum())
    n_used = int(np.isfinite(inj.base["e"]).sum())
    assert c["n_veto_used"] > 0 and n_used == c["n_used_before_veto"] - c["n_veto_used"]
    assert not np.isfinite(inj.base["e"][on]).any()  # vetoed rows leave the statistic
    assert not inj.background[on].any()  # and are never lensed
    np.testing.assert_array_equal(inj.background[~on], plain.background[~on])
    with pytest.raises(ValueError, match="bool array"):
        ish.es.shear_screen(model, shapes, args, 1.0, None, on.astype(int))


def test_null_lens_model_has_no_shear_and_the_lens_model_frame():
    ref = _model()
    null = ish.ir.NullLensModel(ref.ra0, ref.dec0)
    ra, dec = ref.to_sky(np.array([12.0, -30.0]), np.array([5.0, 40.0]))
    np.testing.assert_allclose(null.to_frame(ra, dec), ref.to_frame(ra, dec))
    np.testing.assert_allclose(null.to_sky(*null.to_frame(ra, dec)), (ra, dec))
    g, kappa = ish.es.cluster_reduced_shear(null, ra, dec, np.array([1.0, 3.0]))
    assert np.all(g == 0) and np.all(kappa == 0) and null.z_lens == 0.0


def test_a_massive_w1_lens_is_recovered_without_a_lens_model():
    model, shapes, args = _field()
    null = ish.ir.NullLensModel(model.ra0, model.dec0)
    args.responsivity = 0.42
    inj = ish.ShearInjector(null, shapes, args, psf_sigma=1.0)
    assert inj.r == 0.42
    np.testing.assert_allclose(inj.e[inj.use], inj.eps[inj.use])  # nothing removed
    inj.draw_null(np.random.default_rng(1))
    hit = inj.trial(15.0, -10.0, np.where(inj.background, 6.0, np.nan), SCALE)
    assert hit["recovered"]


def test_ncf_fields_need_their_fixed_responsivity_and_a_lens_redshift(tmp_path):
    import json

    spec = ish.field_spec("macs0416-ncf")
    assert spec["model"] is None and spec["responsivity"] == "ncf"
    assert spec["catalog"] == ish.ncf_catalog("o023_t002")
    assert spec["catalog"].endswith("jw01208-o023_t002_nircam_clear-f150w_cat.ecsv")
    gate = tmp_path / "gate_r.json"
    gate.write_text(json.dumps({"gate_r": {"passed": True, "adopted_r": 0.43}}))
    assert ish.fixed_responsivity("macs0416-ncf", spec, 0.3, gate) == 0.43
    assert ish.fixed_responsivity("macs0416", ish.field_spec("macs0416"), None, gate) is None
    with pytest.raises(ValueError, match="lens-z"):
        ish.fixed_responsivity("macs0416-ncf", spec, None, gate)
    gate.write_text(json.dumps({"gate_r": {"passed": False, "adopted_r": None}}))
    with pytest.raises(ValueError, match="gate R"):
        ish.fixed_responsivity("macs0416-ncf", spec, 0.3, gate)
    with pytest.raises(ValueError, match="gate R"):
        ish.fixed_responsivity("macs0416-ncf", spec, 0.3, tmp_path / "missing.json")
    assert ish.setting(  # the fixed R is per field kind, not a screen setting to compare
        {"assumptions": {"screen": {"max_radius": 1.0}, "fixed_responsivity": 0.43}}
    ) == ish.setting({"assumptions": {"screen": {"max_radius": 2.0}}})


def test_gate_r_rule():
    _ld = importlib.util.spec_from_file_location("w1_ncf_gate_r", _DIR / "w1_ncf_gate_r.py")
    g = importlib.util.module_from_spec(_ld)
    _ld.loader.exec_module(g)
    r200 = {"a": 0.49, "b": 0.60}
    ok = g.gate_r({"a": 0.42, "b": 0.61}, r200)
    assert ok["passed"] and ok["adopted_r"] == 0.42
    assert not g.gate_r({"a": 0.30, "b": 0.61}, r200)["passed"]  # 0.19 > 0.15 in one field
    bad = g.gate_r({"a": None, "b": 0.61}, r200)  # an unusable fit fails the gate
    assert not bad["passed"] and bad["adopted_r"] is None
    assert not g.gate_r({"a": 0.42}, r200)["passed"]  # a field without a fit fails, no KeyError


def test_goodsn_fields_use_the_frozen_r_and_the_shared_dja_zout():
    import json

    names = sorted(ish.GOODSN_FIELDS)
    assert len(names) == 6
    for name in names:
        spec = ish.field_spec(name)
        assert spec["model"] is None and spec["responsivity"] == ish.GOODSN_R
        assert spec["catalog"].endswith("_nircam_clear-f200w_cat.ecsv")
        assert "jw01181-" in spec["catalog"]
        assert spec["photoz"] == ("eazy", ish._GOODSN_ZOUT)
        assert ish.fixed_responsivity(name, spec, 0.3) == pytest.approx(0.414)
        with pytest.raises(ValueError, match="lens-z"):
            ish.fixed_responsivity(name, spec, None)
    # the frozen R is the Abell 370 F200W fit, rounded (results/w1_shear_goodsn/PLAN.md)
    a370 = _DIR.parent / "results" / "w1_volume" / "summaries" / "z0.3_abell370.json"
    assert round(json.loads(a370.read_text())["screen"]["responsivity"], 3) == ish.GOODSN_R


def test_gate_p_rule():
    assert ish.passes_gate_p({"n_used": 100, "n_used_photoz": 50})
    assert not ish.passes_gate_p({"n_used": 100, "n_used_photoz": 49})
    assert not ish.passes_gate_p({"n_used": 0, "n_used_photoz": 0})
    assert not ish.passes_gate_p({"n_used": 100, "n_used_photoz": None})
    assert not ish.passes_gate_p({"n_used_before_veto": 100})  # older summaries: no gate data


def test_fixed_r_sensitivity_runs_are_refused_by_w1_volume(tmp_path):
    import json

    _lv = importlib.util.spec_from_file_location("w1_volume", _DIR / "w1_volume.py")
    wv = importlib.util.module_from_spec(_lv)
    _lv.loader.exec_module(wv)
    d = tmp_path / "z0.3" / "f"
    d.mkdir(parents=True)
    summ = {"assumptions": {"w1_lens_z": 0.3, "fixed_r_override": 0.231}}
    (d / "summary.json").write_text(json.dumps(summ))
    with pytest.raises(SystemExit, match="sensitivity"):
        wv.main(
            ["--runs", str(tmp_path), "--fields", "f", "--z", "0.3", "--out", str(tmp_path / "o")]
        )
    with pytest.raises(ValueError, match="'ncf' or a number"):
        ish.fixed_responsivity("x", {"responsivity": "NCF"}, 0.3)


def test_ceers_fields_use_the_frozen_r_and_the_shared_dja_zout():
    names = sorted(ish.CEERS_FIELDS)
    assert len(names) == 10 and "ceers-o052" not in names  # o052 overlaps o002 / o064 (PLAN.md)
    for name in names:
        spec = ish.field_spec(name)
        assert spec["model"] is None and spec["responsivity"] == ish.GOODSN_R
        assert "jw01345-" in spec["catalog"]
        assert spec["catalog"].endswith("_nircam_clear-f200w_cat.ecsv")
        assert spec["photoz"] == ("eazy", ish._CEERS_ZOUT)
        assert ish.fixed_responsivity(name, spec, 0.3) == pytest.approx(0.414)


def test_primer_tiles_are_fields_with_a_box_and_the_stacked_zout():
    names = sorted(ish.PRIMER_FIELDS)
    assert sum(n.startswith("primer-uds-") for n in names) == 12  # PLAN.md: 12 UDS + 8 COSMOS
    assert sum(n.startswith("primer-cosmos-") for n in names) == 8
    for name in names:
        spec = ish.field_spec(name)
        assert spec["model"] is None and spec["box_half"] == 140.0
        assert spec["responsivity"] == ish.GOODSN_R and spec["max_radius"] == 200.0
        assert "jw01837-" in spec["catalog"]
        survey = name.split("-")[1]
        assert spec["photoz"] == ("eazy", ish.PRIMER_ZOUT.format(survey))
        assert np.hypot(spec["box_half"], spec["box_half"]) < spec["max_radius"]
