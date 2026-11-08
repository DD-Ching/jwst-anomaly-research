"""`arcs` anti null with an empirical blank-field shape prior (exploratory).

Plan: results/arcs/shape-prior/PLAN.md (frozen before the prior was read). The Gaussian σ_ε 0.26 null of
``lens_consistency.anti_probability`` is replaced by |ε| resampled from the CEERS F200W pipeline catalogue
(same cuts as the `arcs` runs except the ellipticity cut), with uniform phases (ASSUMPTION). The
lensing map, the e ≥ 0.5 and 60° anti cuts and the g_min rounding are those of ``anti_probability``.
The catalogue is read from S3, pinned by sha256 in ``data/manifests/reference_sample.ecsv``, and
deleted after reading unless ``--keep``.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import tempfile
from functools import cache
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack
from scipy.stats import poisson

sys.path.insert(0, str(Path(__file__).resolve().parent))
import arcs_bands as ab  # noqa: E402
import arcs_null_calibration as nc  # noqa: E402
import arcs_pooled as ap  # noqa: E402

from jwst_anomaly.photometry import fetch_catalog  # noqa: E402

lc = ap.lc
OUT = ap.ROOT / "shape-prior"
PRIOR_FILE = "jw01345-o001_t021_nircam_clear-f200w_cat.ecsv"
MANIFEST = "reference_sample"
# frozen (PLAN.md; ASSUMPTIONS)
ALPHA = 0.05
GATE_P = 0.05
ISO_TOL = 0.01
OPEN_ROWS = {"WHL0137": 2194, "RXC J0600": 581, "MACS0308": 2157}
N_DRAW = 200_000
SEED = 1


def eps_from_e(e: np.ndarray) -> np.ndarray:
    """|ε| = (a − b)/(a + b) from e = 1 − b/a."""
    e = np.asarray(e, float)
    return e / (2.0 - e)


def prior_eps(shapes: Table, min_semimajor_px: float, min_snr: float) -> np.ndarray:
    """|ε| of the blank-field sources passing the `arcs` cuts other than the ellipticity cut."""
    sel = (
        np.asarray(shapes["is_extended"], bool)
        & (np.asarray(shapes["semimajor_px"], float) >= min_semimajor_px)
        & (np.asarray(shapes["snr"], float) >= min_snr)
        & np.isfinite(np.asarray(shapes["ellipticity"], float))
    )
    eps = eps_from_e(np.asarray(shapes["ellipticity"], float)[sel])
    return eps[(eps >= 0) & (eps < 0.95)]


def p_anti_emp(
    g: float,
    eps: np.ndarray,
    min_ellipticity: float = 0.5,
    anti_deg: float = 60.0,
    n: int = N_DRAW,
    seed: int = SEED,
) -> float:
    """``anti_probability`` with |ε| resampled from ``eps`` and uniform phases."""
    rng = np.random.default_rng(seed)
    mod = rng.choice(np.asarray(eps, float), n, replace=True)
    es = mod * np.exp(1j * rng.uniform(0.0, 2 * np.pi, n))
    g = float(g)
    eo = (es + g) / (1 + g * es) if abs(g) <= 1 else (1 + g * np.conj(es)) / (np.conj(es) + g)
    e_obs = 2 * np.abs(eo) / (1 + np.abs(eo))  # 1 - b/a
    angle = np.degrees(np.abs(np.angle(eo))) / 2
    sel = e_obs >= min_ellipticity
    return float(np.mean(angle[sel] >= anti_deg)) if sel.any() else float("nan")


def gate(weak_k: int, weak_lam: float, p_iso: float) -> dict:
    """Rule 3: weak-row two-sided P ≥ GATE_P and isotropy at g = 0 within ISO_TOL."""
    p_two = nc.two_sided_p(weak_k, weak_lam)
    return {
        "weak_anti": weak_k,
        "weak_expected": weak_lam,
        "weak_p_two_sided": p_two,
        "p_anti_g0": p_iso,
        "passes": bool(p_two >= GATE_P and abs(p_iso - 1.0 / 3.0) <= ISO_TOL),
    }


def outcome(gate_passes: bool, p_b: float, alpha: float = ALPHA) -> str:
    """Rule 4 (frozen)."""
    if not gate_passes:
        return "tool fails; undetermined"
    if p_b >= alpha:
        return "pooled excess explained by the empirical shape prior"
    return "excess persists under the empirical prior"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--keep", action="store_true", help="keep the downloaded catalogue")
    args = p.parse_args()
    meta = json.loads((ap.ROOT / ap.S1063).parent.joinpath("arcs.json").read_text())["assumptions"]
    paths = {**{f: ap.ROOT / q for f, q in ap.FIELDS.items()}, "S1063": ap.ROOT / ap.S1063}
    for f, q in paths.items():
        if json.loads(q.parent.joinpath("arcs.json").read_text())["assumptions"] != meta:
            raise SystemExit(f"error: {f} uses different `arcs` assumptions")
    e_min, anti_deg = float(meta["min_ellipticity"]), float(meta["anti_deg"])

    sha = ab.manifest_sha(MANIFEST, PRIOR_FILE)
    tmp = Path(tempfile.mkdtemp(prefix="shape-prior-"))
    try:
        path = fetch_catalog(ab.s3_url(PRIOR_FILE), sha, cache_dir=tmp)
        shapes = lc.load_shapes(path)
        if args.keep:
            shutil.copy(path, OUT / path.name)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    eps = prior_eps(shapes, float(meta["min_semimajor_px"]), float(meta["min_snr"]))

    @cache
    def pe(g: float) -> float:
        return p_anti_emp(g, eps, e_min, anti_deg)

    def lam(t: Table) -> float:
        return float(sum(pe(round(float(g), 2)) for g in np.asarray(t["g_min"], float)))

    weak = {f: nc.rows(q, strong=False) for f, q in paths.items()}
    strong = {f: nc.rows(q, strong=True) for f, q in paths.items()}

    def pool(tabs: dict[str, Table], fields: tuple[str, ...]) -> Table:
        return vstack(
            [tabs[f]["label", "g_min", "orientation_class"] for f in fields],
            metadata_conflicts="silent",
        )

    wb = pool(weak, ap.VARIANTS["B"])
    g_res = gate(int(ap.anti_mask(wb).sum()), lam(wb), pe(0.0))

    res = {}
    for v, fields in ap.VARIANTS.items():
        t = pool(strong, fields)
        k, lam_e = int(ap.anti_mask(t).sum()), lam(t)
        lam_g = float(nc.p_anti(t, nc.SIGMA0, e_min, anti_deg).sum())
        res[v] = {
            "n": len(t),
            "anti": k,
            "expected_emp": lam_e,
            "p_emp": float(poisson.sf(k - 1, lam_e)),
            "expected_gauss_0.26": lam_g,
        }
    s = strong["S1063"]
    k = int(ap.anti_mask(s).sum())
    res["S1063"] = {
        "n": len(s),
        "anti": k,
        "expected_emp": lam(s),
        "p_emp": float(poisson.sf(k - 1, lam(s))),
    }
    weak_a = pool(weak, ap.VARIANTS["A"])
    res["weak_A"] = {"anti": int(ap.anti_mask(weak_a).sum()), "expected_emp": lam(weak_a)}

    rows_out = {}
    for f, label in OPEN_ROWS.items():
        t = strong[f]
        hit = np.flatnonzero(np.asarray(t["label"]) == label)
        if len(hit) != 1:
            raise SystemExit(f"error: {f} {label} not once in the strong-shear rows")
        r = t[int(hit[0])]
        e_row = float(eps_from_e(float(r["ellipticity"])))
        g = round(float(r["g_min"]), 2)
        rows_out[f"{f} {label}"] = {
            "g_min": g,
            "ellipticity": float(r["ellipticity"]),
            "eps": e_row,
            "p_anti_emp": pe(g),
            "p_anti_gauss": ap._p_anti(g, nc.SIGMA0, e_min, anti_deg),
            "prior_frac_eps_ge_row": float(np.mean(eps >= e_row)),
        }

    gx = np.random.default_rng(SEED).normal(0.0, nc.SIGMA0, (2, N_DRAW))
    gauss_frac = float(np.mean(np.hypot(gx[0], gx[1]) >= 1 / 3))
    summary = {
        "plan": "results/arcs/shape-prior/PLAN.md",
        "inputs_sha256": {PRIOR_FILE: sha},
        "params": {
            "alpha": ALPHA,
            "gate_p": GATE_P,
            "iso_tol": ISO_TOL,
            "n_draw": N_DRAW,
            "seed": SEED,
            **{
                k2: meta[k2]
                for k2 in ("min_ellipticity", "anti_deg", "min_semimajor_px", "min_snr")
            },
        },
        "prior": {
            "n": int(len(eps)),
            "eps_median": float(np.median(eps)),
            "eps_rms_per_component": float(np.sqrt(np.mean(eps**2) / 2)),
            "frac_e_ge_0.5": float(np.mean(eps >= 1 / 3)),
            "gauss_0.26_frac_e_ge_0.5": gauss_frac,
            "eps_quantiles": {
                str(q): float(np.quantile(eps, q)) for q in (0.5, 0.75, 0.9, 0.95, 0.99)
            },
        },
        "gate": g_res,
        "strong": res,
        "open_rows": rows_out,
        "outcome": outcome(g_res["passes"], res["B"]["p_emp"]),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
