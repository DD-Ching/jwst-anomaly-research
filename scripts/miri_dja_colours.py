"""SMACS MIRI colours from DJA fixed apertures vs pipeline EE apertures.

Plan: results/miri_dja_colours/PLAN.md. Reads a run of ``configs/smacs_jwst3_miri_dja.yaml``
(sources + ``photometry.ecsv`` with ``dja07``) and, optionally, a baseline run of
``configs/smacs_jwst3_miri_confirm.yaml`` for the top-20 Jaccard. The frozen rules V, C and B are
the functions below; ``decide`` returns the verdict.

    python scripts/miri_dja_colours.py --run-dir <runs/ID> \
        [--baseline-run-dir <runs/ID>] --out summary.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table, join
from scipy.stats import spearmanr

from jwst_anomaly import schema
from jwst_anomaly.features import column_as_float, snr_from_mag_err

SAMPLE = "smacs0723_miri"
LABEL = "dja07"
V_MAX_MAD = 0.15  # ASSUMPTION (PLAN.md rule V)
C_MIN_COVERAGE = 0.80  # ASSUMPTION (rule C)
B_MIN_DROP = 0.10  # ASSUMPTION (rule B)
B_MIN_RHO = 0.10  # ASSUMPTION (rule B: below this there is no bias to remove)
MIN_SNR = 10.0


def mad(x: np.ndarray) -> float:
    x = np.asarray(x, float)
    return float(np.median(np.abs(x - np.median(x)))) if x.size else float("nan")


def rule_v(diff_compact: np.ndarray, max_mad: float = V_MAX_MAD) -> bool:
    """V: MAD of (dja07 − aper50) colour in the lowest-area quartile ≤ max_mad (needs ≥ 5 rows)."""
    d = np.asarray(diff_compact, float)
    d = d[np.isfinite(d)]
    return bool(d.size >= 5 and mad(d) <= max_mad)


def rule_c(n_ranked_with_colour: int, n_ranked: int, min_cov: float = C_MIN_COVERAGE) -> bool:
    """C: fraction of ranked sources with a finite dja07 colour ≥ min_cov."""
    return bool(n_ranked > 0 and n_ranked_with_colour / n_ranked >= min_cov)


def rule_b(rho_aper50: float, rho_dja: float) -> str:
    """B: 'no_bias' if |ρ_aper50| < B_MIN_RHO; 'pass' if |ρ_dja| ≤ |ρ_aper50| − B_MIN_DROP.

    Otherwise (or with a non-finite ρ) 'fail'.
    """
    if not np.isfinite(rho_aper50):
        return "fail"
    if abs(rho_aper50) < B_MIN_RHO:
        return "no_bias"
    if not np.isfinite(rho_dja):
        return "fail"
    return "pass" if abs(rho_dja) <= abs(rho_aper50) - B_MIN_DROP else "fail"


def decide(v: bool, c: bool, b: str) -> str:
    """Adopt dja07 only if V, C and B all pass; a failed V stops the evaluation."""
    if not v:
        return "stopped_v_failed"
    return "adopt" if (c and b == "pass") else "not_adopted"


def jaccard(a: list[str], b: list[str]) -> float:
    sa, sb = set(a), set(b)
    return len(sa & sb) / len(sa | sb) if sa | sb else float("nan")


def _col(t: Table, band: str, q: str) -> np.ndarray:
    return column_as_float(t, schema.band_column(band, q))


def colour_rows(t: Table, a: str, b: str) -> dict:
    """Common rows (both detected, ``a`` aper50 S/N > MIN_SNR, finite in all three colours)."""
    det = (_col(t, a, "detected") > 0) & (_col(t, b, "detected") > 0)
    snr = snr_from_mag_err(_col(t, a, "aper50_abmag_err"))
    area = _col(t, "F770W", "isophotal_area")
    cols = {
        n: _col(t, a, f"{n}_abmag") - _col(t, b, f"{n}_abmag")
        for n in ["aper50", "isophotal", LABEL]
    }
    ok = det & (snr > MIN_SNR) & np.isfinite(area) & (area > 0)
    for v in cols.values():
        ok &= np.isfinite(v)
    out: dict = {"n_common": int(ok.sum())}
    for n, v in cols.items():
        if ok.sum() >= 5:
            rho, p = spearmanr(np.log10(area[ok]), v[ok])
            out[n] = {"rho": float(rho), "p": float(p), "median": float(np.median(v[ok]))}
    q25 = np.quantile(area[ok], 0.25) if ok.any() else np.nan
    compact = ok & (area <= q25)
    out["compact_quartile_area_max_px"] = float(q25)
    out["compact_n"] = int(compact.sum())
    out["compact_diff_dja_minus_aper50"] = {
        "median": float(np.median((cols[LABEL] - cols["aper50"])[compact]))
        if compact.any()
        else None,
        "mad": mad((cols[LABEL] - cols["aper50"])[compact]),
    }
    out["_compact_diff"] = (cols[LABEL] - cols["aper50"])[compact]
    return out


def load(run_dir: Path) -> Table:
    t = Table.read(run_dir / SAMPLE / "sources.ecsv")
    phot_path = run_dir / SAMPLE / "photometry.ecsv"
    if phot_path.exists() and f"f770w_{LABEL}_abmag" not in t.colnames:
        t = join(t, Table.read(phot_path), keys="source_uid", join_type="left")
    return t


def ranked(run_dir: Path) -> Table:
    """Ranked (gate-passing, non-star) sources of the run, best first."""
    r = Table.read(run_dir / SAMPLE / "scores.ecsv")
    r.sort("rank")
    return r


def top20(run_dir: Path) -> list[str]:
    """The run's reported top 20 after the cutout screens (D-019 to D-021): ``targets.ecsv``."""
    t = Table.read(run_dir / SAMPLE / "targets.ecsv")
    t.sort("rank")
    return [str(u) for u in t["source_uid"][:20]]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--baseline-run-dir", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args(argv)
    t = load(args.run_dir)
    rk_tab = ranked(args.run_dir)
    rk = join(rk_tab[["source_uid"]], t, keys="source_uid", join_type="left")
    c770 = _col(rk, "F770W", f"{LABEL}_abmag") - _col(rk, "F1000W", f"{LABEL}_abmag")
    n_cov = int(np.isfinite(c770).sum())
    res: dict = {
        "run_id": args.run_dir.name,
        "provenance": "derived (observed DJA v7.4 and MAST L3 jwst 3.0.0 MIRI catalogues)",
        "thresholds_assumption": {
            "V_max_mad": V_MAX_MAD,
            "C_min_coverage": C_MIN_COVERAGE,
            "B_min_drop": B_MIN_DROP,
            "B_min_rho": B_MIN_RHO,
            "min_snr_f770w_aper50": MIN_SNR,
        },
        "n_sources": len(t),
        "n_ranked": len(rk_tab),
        "n_ranked_with_dja07_f770w_f1000w": n_cov,
        "coverage": n_cov / len(rk_tab) if len(rk_tab) else None,
        "colours": {},
    }
    for a, b in [("F770W", "F1000W"), ("F1000W", "F1500W")]:
        res["colours"][f"{a}-{b}"] = colour_rows(t, a, b)
    main_c = res["colours"]["F770W-F1000W"]
    v = rule_v(main_c.pop("_compact_diff"))
    res["colours"]["F1000W-F1500W"].pop("_compact_diff")
    c = rule_c(n_cov, len(rk_tab))
    b = rule_b(
        main_c.get("aper50", {}).get("rho", np.nan), main_c.get(LABEL, {}).get("rho", np.nan)
    )
    res["rules"] = {"V": v, "C": c, "B": b}
    res["verdict"] = decide(v, c, b)
    res["top20"] = top20(args.run_dir)
    if args.baseline_run_dir:
        base = top20(args.baseline_run_dir)
        res["baseline_run_id"] = args.baseline_run_dir.name
        res["baseline_top20"] = base
        res["top20_jaccard"] = jaccard(res["top20"], base)
        res["top20_new"] = [u for u in res["top20"] if u not in base]
        res["top20_dropped"] = [u for u in base if u not in res["top20"]]
    args.out.write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps({k: res[k] for k in ["coverage", "rules", "verdict"]}, indent=1))
    print(json.dumps(main_c, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
