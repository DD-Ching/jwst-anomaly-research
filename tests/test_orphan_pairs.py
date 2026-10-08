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
