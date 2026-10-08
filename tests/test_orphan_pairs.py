"""Tests for scripts/orphan_pairs.py (synthetic catalogue only)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("orphan_pairs", _DIR / "orphan_pairs.py")
op = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(op)

BANDS = ("F435W", "F606W", "F814W", "F090W", "F115W", "F150W", "F200W", "F277W", "F356W", "F444W")
WAVE = np.array([0.435, 0.606, 0.814, 0.90, 1.15, 1.50, 2.00, 2.77, 3.56, 4.44])
RA0, DEC0 = 64.04, -24.07
Z_CLUSTER = 0.4


def _sed(slope: float, bump: float, amp: float) -> np.ndarray:
    """A smooth toy SED (nJy): power law plus a Gaussian bump in log wavelength."""
    lw = np.log(WAVE)
    return amp * WAVE**slope * (1 + 2.0 * np.exp(-0.5 * ((lw - bump) / 0.3) ** 2))


def synthetic_catalogue(seed: int = 3) -> tuple[Table, dict[str, int]]:
    """Random background sources, one injected identical-SED pair and one random pair."""
    rng = np.random.default_rng(seed)
    n = 150
    x = rng.uniform(-100, 100, n)  # arcsec East
    y = rng.uniform(-100, 100, n)
    seds = [_sed(rng.uniform(-2, 2), rng.uniform(-0.5, 1.5), rng.uniform(300, 3000)) for _ in x]
    # injected: one SED at two places 1.2" apart (fluxes 1 : 0.5), far from everything else
    x = np.r_[x, 150.0, 150.8, -150.0, -149.2]
    y = np.r_[y, 150.0, 150.9, -150.0, -149.1]
    base = _sed(0.5, 0.8, 1500.0)
    seds += [base, 0.5 * base, _sed(-1.5, 0.0, 1500.0), _sed(1.5, 1.4, 1500.0)]
    flux = np.array(seds)
    err = 0.02 * flux + 5.0
    flux = flux + rng.normal(0, 1, flux.shape) * err
    m = len(x)
    cat = Table(
        {
            "SOURCE": np.arange(1, m + 1),
            "RA": RA0 + x / 3600.0 / np.cos(np.deg2rad(DEC0)),
            "DEC": DEC0 + y / 3600.0,
        }
    )
    px, py = -x / op.PIXSCALE, y / op.PIXSCALE  # pixel frame: x towards West
    cat["X"], cat["Y"] = px, py
    cat["X_MIN"], cat["X_MAX"], cat["Y_MIN"], cat["Y_MAX"] = px - 3, px + 3, py - 3, py + 3
    cat["KRON_RADIUS"], cat["A"], cat["PA"] = np.full(m, 1.0), np.full(m, 2.0), np.zeros(m)
    cat["AREA_ISO"] = np.full(m, 20.0)
    for col in ("FLAG_DEBLEND", "FLAG_BCG", "FLAG_POINTSRC"):
        cat[col] = np.zeros(m, bool)
    cat["USE_PHOT_APER03"] = np.ones(m, bool)
    for col, z in (("Z_ML", 2.0), ("Z025", 1.6), ("Z160", 1.8), ("Z840", 2.2), ("Z975", 2.4)):
        cat[col] = np.full(m, z)
    cat["Z_SPEC"], cat["MU"] = np.full(m, -99.0), np.ones(m)
    for k, b in enumerate(BANDS):
        cat[f"FLUX_COLOR03_TOTAL_{b}"] = flux[:, k]
        cat[f"FLUXERR_COLOR03_TOTAL_{b}"] = err[:, k]
    return cat, {"injected": (n + 1, n + 2), "random": (n + 3, n + 4)}


def test_sed_chi2_recovers_the_scale_of_an_identical_sed():
    f = _sed(0.3, 0.5, 1000.0)[None, :]
    e = 0.01 * f
    chi2, nb, scale, n_agree = op.sed_chi2(f, e, 0.25 * f, 0.25 * e)
    assert nb[0] == len(BANDS) and n_agree[0] == len(BANDS)
    assert chi2[0] < 1e-6 and abs(scale[0] - 4.0) < 1e-6


def test_search_finds_the_injected_pair_and_rejects_the_random_one():
    cat, ids = synthetic_catalogue()
    images = Table({"image_id": [], "ra": [], "dec": []}, dtype=[str, float, float])
    pairs, summary = op.search(cat, Z_CLUSTER, images, n_shift=3)
    found = {tuple(sorted((int(a), int(b)))) for a, b in pairs["id_a", "id_b"]}
    assert ids["injected"] in found
    assert ids["random"] not in found
    sel = [tuple(sorted((int(a), int(b)))) == ids["injected"] for a, b in pairs["id_a", "id_b"]]
    row = pairs[sel]
    assert row["pair_class"][0] == "orphan"
    assert summary["n_selected"] == len(cat)
    ranked = op.rank_orphans(pairs)
    assert ranked["rank"][0] == 1


def test_published_image_and_visible_lens_reclassify_the_pair():
    cat, ids = synthetic_catalogue()
    a = cat[cat["SOURCE"] == ids["injected"][0]]
    images = Table({"image_id": ["1.1"], "ra": [a["RA"][0]], "dec": [a["DEC"][0]]})
    pairs, _ = op.search(cat, Z_CLUSTER, images, n_shift=1)
    sel = [tuple(sorted((int(x), int(y)))) == ids["injected"] for x, y in pairs["id_a", "id_b"]]
    assert pairs["pair_class"][sel][0] == "published"

    # a faint catalogued galaxy (any redshift, not selected) at the midpoint: a visible lens
    b = cat[cat["SOURCE"] == ids["injected"][1]]
    lens = cat[:1].copy()
    lens["SOURCE"], lens["USE_PHOT_APER03"] = 9999, False
    lens["RA"], lens["DEC"] = 0.5 * (a["RA"][0] + b["RA"][0]), 0.5 * (a["DEC"][0] + b["DEC"][0])
    from astropy.table import vstack

    empty = Table({"image_id": [], "ra": [], "dec": []}, dtype=[str, float, float])
    pairs, _ = op.search(vstack([cat, lens]), Z_CLUSTER, empty, n_shift=1)
    sel = [tuple(sorted((int(x), int(y)))) == ids["injected"] for x, y in pairs["id_a", "id_b"]]
    assert pairs["pair_class"][sel][0] == "visible_lens"


def test_dark_lens_numbers_scale_with_separation():
    t = op.dark_lens_numbers([1.0, 2.0], 0.4, [2.0, 2.0])
    assert np.allclose(t["theta_e"], [0.5, 1.0])
    assert np.isclose(t["mass_e_msun"][1] / t["mass_e_msun"][0], 4.0)
    assert 1e10 < t["mass_e_msun"][0] < 1e12  # theta_E = 0.5" at z_l = 0.4: ~1e11 Msun


def test_lens_rules_ignore_fragments_of_a_member_and_report_only_qualifying_lenses():
    cat, ids = synthetic_catalogue()
    a = cat[cat["SOURCE"] == ids["injected"][0]]
    # a fragment 0.1" from member a: not a lens (not clear of both members)
    frag = cat[:1].copy()
    frag["SOURCE"], frag["USE_PHOT_APER03"] = 9998, False
    frag["RA"], frag["DEC"] = a["RA"][0], a["DEC"][0] + 0.1 / 3600.0
    from astropy.table import vstack

    empty = Table({"image_id": [], "ra": [], "dec": []}, dtype=[str, float, float])
    pairs, _ = op.search(vstack([cat, frag]), Z_CLUSTER, empty, n_shift=1)
    sel = [tuple(sorted((int(x), int(y)))) == ids["injected"] for x, y in pairs["id_a", "id_b"]]
    row = pairs[sel][0]
    assert row["pair_class"] == "orphan"
    assert row["lens_source"] == 0


def test_select_sources_drops_invalid_long_wavelength_photometry():
    cat, ids = synthetic_catalogue()
    k = int(np.flatnonzero(cat["SOURCE"] == ids["injected"][0])[0])
    assert op.select_sources(cat, Z_CLUSTER)[k]
    cat["FLUXERR_COLOR03_TOTAL_F444W"][k] = -5.0  # flagged invalid: no summed S/N
    assert not op.select_sources(cat, Z_CLUSTER)[k]


def test_summary_is_strict_json():
    import json

    raw = {"a": float("nan"), "b": [np.float64(np.inf), 1.0], "c": np.array([np.nan, 2.0])}
    raw["d"] = np.ma.masked
    text = json.dumps(op._finite(raw), allow_nan=False)
    assert json.loads(text) == {"a": None, "b": [None, 1.0], "c": [None, 2.0], "d": None}


# ------------------------------------------------------------------- column abstraction (D-051)


def dja_from_canucs(cat: Table) -> tuple[Table, Table]:
    """The same synthetic catalogue in DJA grizli layout (uJy apertures) plus its zout."""
    import astropy.units as u

    n = len(cat)
    d = Table({"id": np.asarray(cat["SOURCE"]), "ra": cat["RA"], "dec": cat["DEC"]})
    d["x"], d["y"] = cat["X"], cat["Y"]
    for a, b in (("xmin", "X_MIN"), ("xmax", "X_MAX"), ("ymin", "Y_MIN"), ("ymax", "Y_MAX")):
        d[a] = cat[b]
    # half-light radius (px) that gives the CANUCS Kron aperture radius under DJA_KRON_PER_R50
    d["flux_radius"] = op.aperture_radius(cat) / (op.DJA_KRON_PER_R50 * 0.04)
    d["area_iso"], d["theta_image"] = cat["AREA_ISO"], np.deg2rad(cat["PA"])
    d["flag"] = np.where(cat["FLAG_DEBLEND"], 1, 0)
    d["flag_aper_0"] = np.zeros(n, int)
    for b in BANDS:
        lo = b.lower()
        d[f"{lo}_flux_aper_0"] = np.asarray(cat[f"FLUX_COLOR03_TOTAL_{b}"]) / 1000.0 * u.uJy
        d[f"{lo}_fluxerr_aper_0"] = np.asarray(cat[f"FLUXERR_COLOR03_TOTAL_{b}"]) / 1000.0 * u.uJy
        d[f"{lo}_flag_aper_0"] = np.zeros(n, int)
    # a MIRI band and a "u" duplicate: both must be left out of the SED
    for extra in ("f770w", "f444wu"):
        d[f"{extra}_flux_aper_0"] = np.ones(n) * u.uJy
        d[f"{extra}_fluxerr_aper_0"] = np.ones(n) * u.uJy
    d.meta.update(ASEC_0=0.36, APER_0=9.0, KRONFACT=2.5)
    z = Table({"id": np.asarray(cat["SOURCE"])})
    z["z_spec"], z["z_ml"] = np.asarray(cat["Z_SPEC"]), np.asarray(cat["Z_ML"])
    for a, b in (("z025", "Z025"), ("z160", "Z160"), ("z840", "Z840")):
        z[a] = np.asarray(cat[b])
    return d, z


def test_canucs_and_dja_layouts_give_the_same_standard_table_and_pairs():
    cat, ids = synthetic_catalogue()
    std_c = op.canucs_standard(cat)
    std_d = op.dja_standard(*dja_from_canucs(cat))
    assert std_c.meta["bands"] == std_d.meta["bands"] == list(BANDS)
    for col in ("flux", "err", "ra", "dec", "z_low", "z16", "z84", "ap_radius", "pa_deg"):
        assert np.allclose(std_c[col], std_d[col], equal_nan=True), col
    assert np.array_equal(std_c["deblend"], std_d["deblend"])
    images = Table({"image_id": [], "ra": [], "dec": []}, dtype=[str, float, float])
    pc, _ = op.search(std_c, Z_CLUSTER, images, n_shift=2)
    pd, _ = op.search(std_d, Z_CLUSTER, images, n_shift=2)
    assert list(pc["id_a"]) == list(pd["id_a"])
    assert list(pc["pair_class"]) == list(pd["pair_class"])
    # as_standard is idempotent and the raw CANUCS table converts implicitly
    assert op.as_standard(std_c) is std_c
    assert np.array_equal(op.select_sources(cat, Z_CLUSTER), op.select_sources(std_c, Z_CLUSTER))


def test_dja_flags_invalidate_a_band_and_misaligned_zout_is_refused():
    import pytest

    cat, _ = synthetic_catalogue()
    d, z = dja_from_canucs(cat)
    d["f444w_flag_aper_0"][0] = 0x10  # APER_TRUNC
    std = op.dja_standard(d, z)
    assert np.isnan(std["flux"][0, BANDS.index("F444W")])
    with pytest.raises(ValueError, match="row-aligned"):
        op.dja_standard(d, z[::-1])


def test_snr_bands_fall_back_when_f356w_is_missing():
    cat, ids = synthetic_catalogue()
    cat.remove_columns(["FLUX_COLOR03_TOTAL_F356W", "FLUXERR_COLOR03_TOTAL_F356W"])
    assert op.snr_bands_of(cat) == ["F277W", "F444W"]
    k = int(np.flatnonzero(cat["SOURCE"] == ids["injected"][0])[0])
    assert op.select_sources(cat, Z_CLUSTER)[k]


def test_deep_field_lens_check_needs_no_model():
    cat, _ = synthetic_catalogue()
    images = Table({"image_id": [], "ra": [], "dec": []}, dtype=[str, float, float])
    pairs, summary = op.search(cat, op.Z_LENS_REF, images, n_shift=1)
    top = op.lens_check(op.rank_orphans(pairs), None, op.Z_LENS_REF)
    assert len(top) and not np.any(top["cluster_explains"])
    assert all(c in top.meta["hypothesis_columns"] for c in ("mass_e_msun", "sis_sigma"))
    keys = ("null_d_near_zmatched", "null_e_conditioned", "null_f_near_conditioned")
    for key in (*keys, "zoverlap_match_fraction_by_sep"):
        assert key in summary


def test_fetch_tar_member_streams_and_verifies(tmp_path):
    import hashlib
    import io
    import tarfile

    import pytest

    payload = b"zout bytes" * 100
    arch = tmp_path / "a.tar.gz"
    with tarfile.open(arch, "w:gz") as tar:
        for name, data in (("big.h5", b"x" * 5000), ("cat.zout.fits", payload)):
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    spec = {
        "url": arch.as_uri(),
        "sha256": hashlib.sha256(arch.read_bytes()).hexdigest(),
        "member": "cat.zout.fits",
        "member_sha256": hashlib.sha256(payload).hexdigest(),
        "max_bytes": 10**6,
    }
    out = op.fetch_tar_member(spec, cache_dir=tmp_path / "cache")
    assert out.read_bytes() == payload
    assert op.fetch_tar_member(spec, cache_dir=tmp_path / "cache") == out  # cached
    bad = dict(spec, sha256="0" * 64)
    with pytest.raises(ValueError, match="archive sha256"):
        op.fetch_tar_member(bad, cache_dir=tmp_path / "cache2")
    assert not list((tmp_path / "cache2").glob("*"))  # nothing kept on failure


def test_in_region_handles_polygons():
    poly = "POLYGON 10.0 -1.0 10.1 -1.0 10.1 -0.9 10.0 -0.9"
    assert op._in_region(poly, 10.05, -0.95)
    assert not op._in_region(poly, 10.2, -0.95)
    assert op._in_region("POLYGON 359.95 0.0 0.05 0.0 0.05 0.1 359.95 0.1", 0.0, 0.05)


def test_footprint_covers_the_sources_and_excludes_empty_sky():
    cat, _ = synthetic_catalogue()
    gx, gy, area, _ = op.searched_footprint(cat)
    assert 0 < area < 300.0 * 300.0  # sources span ~300" x 300" but sparsely
    assert len(gx) == int(area)


def test_in_region_accepts_a_frame_token():
    sq = "10.0 -1.0 10.0 1.0 12.0 1.0 12.0 -1.0"
    assert op._in_region(f"POLYGON {sq}", 11.0, 0.0)
    assert op._in_region(f"POLYGON ICRS {sq}", 11.0, 0.0)
    assert not op._in_region(f"POLYGON ICRS {sq}", 13.0, 0.0)


def test_pair_cells_are_symmetric_and_nan_colour_has_its_own_bin():
    cat, _ = synthetic_catalogue()
    sub = op.as_standard(cat)
    src = op.source_cells(sub)
    n = len(sub)
    i, j = np.arange(n - 1), np.arange(1, n)
    fwd = op.pair_cells(src, Table({"i": i, "j": j}))
    rev = op.pair_cells(src, Table({"i": j, "j": i}))
    assert np.array_equal(fwd, rev)
    assert list(op.colour_bin(np.array([np.nan, -0.1, 0.1, 0.4, 0.9, np.inf]))) == [
        4,
        0,
        1,
        2,
        3,
        4,
    ]
    # a pair with a NaN-colour member no longer shares the cell of a 0-0.3 colour member
    snr, rad, _ = src
    nan_cb = np.array([1, 4, 1])  # rows 0, 2: colour 0-0.3; row 1: NaN
    src3 = (snr[:3] * 0 + 50.0, rad[:3] * 0 + 0.5, nan_cb)
    cells = op.pair_cells(src3, Table({"i": [0, 1], "j": [2, 2]}))
    assert cells[0] != cells[1]


def test_conditioned_null_uses_cell_rates_and_counts_fallbacks(monkeypatch):
    monkeypatch.setattr(op, "MIN_CELL_REF", 1)
    ref = Table({"z_overlap": [True, True, True, False], "match": [True, False, True, True]})
    key_ref = np.array([1, 1, 2, 3])
    m = Table({"z_overlap": [True, True, True, False], "pair_class": ["orphan"] * 4})
    key_c = np.array([1, 2, 3, 1])  # cell 3 has no z-overlapping reference pair
    out = op.conditioned_null(key_c, m, key_ref, ref, {"orphan": 2})
    # 0.5 (cell 1) + 1.0 (cell 2) + 2/3 (global fallback); the non-overlapping pair adds nothing
    assert np.isclose(out["expected_by_class"]["orphan"], 0.5 + 1.0 + 2 / 3)
    assert out["n_fallback"] == 1 and out["n_ref_zoverlap"] == 3
    # cells with fewer than MIN_CELL_REF reference pairs fall back to the global rate
    monkeypatch.setattr(op, "MIN_CELL_REF", 2)
    out = op.conditioned_null(key_c, m, key_ref, ref, {"orphan": 2})
    assert np.isclose(out["expected_by_class"]["orphan"], 0.5 + 2 / 3 + 2 / 3)
    assert out["n_fallback"] == 2


def test_conditioned_null_without_reference_is_undefined():
    import pytest

    ref = Table({"z_overlap": [False], "match": [True]})
    m = Table({"z_overlap": [True], "pair_class": ["orphan"]})
    with pytest.warns(UserWarning, match="no z-overlapping"):
        out = op.conditioned_null(np.array([1]), m, np.array([1]), ref, {"orphan": 1})
    assert np.isnan(out["expected_by_class"]["orphan"])
