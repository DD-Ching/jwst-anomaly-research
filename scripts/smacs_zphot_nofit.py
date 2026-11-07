"""What the DJA eazy ``z_phot = -1`` (no-fit) flag marks among SMACS ranked galaxy-stratum sources.

Plan: ``results/smacs_jwst3/zphot_nofit/PLAN.md`` (rule and definitions frozen there before the run).
Reads one pipeline run of ``configs/smacs_jwst3_gaia_only.yaml``, matches ranked galaxy-stratum and star-stratum
sources to the DJA v7.4 zout (0.2", same ``id`` order as the DJA phot catalog), and compares no-fit with fitted
rows: usable filters (``nusefilt``) and size (DJA ``flux_radius``, detection image, pixels). Writes ``summary.json``.

Provenance: catalog columns are ``observed`` / ``derived`` (DJA, MAST); every statistic is ``derived``; the
thresholds are ASSUMPTIONs from the plan.

Usage: python scripts/smacs_zphot_nofit.py <run_dir> --zout <zout.fits> --phot <phot.fits> --out <dir>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
from science_control_features import ZOUT_MATCH_ARCSEC  # noqa: E402
from smacs_zphot_stratify import TOP_K, ZOUT_SHA256  # noqa: E402

from jwst_anomaly.provenance import file_sha256  # noqa: E402

SAMPLE = "smacs0723_nircam"
PHOT_SHA256 = "361a4eb944fefa3afbbbfe08192d2d5c5764ae7137e41d77c2ce63c5b4a47166"  # configs (DJA v7.4 phot)
COV_FRACTION = 0.8  # ASSUMPTION (plan): share of no-fit rows below the fitted 5th percentile of nusefilt
COV_PCT = 5.0  # ASSUMPTION (plan)
SIZE_PCT = 10.0  # ASSUMPTION (plan): no-fit median flux_radius below this fitted percentile
STAR_PCT = 90.0  # ASSUMPTION (plan): "compact" = flux_radius <= this percentile of matched stars
MIN_STARS = 5  # ASSUMPTION (plan)
MIN_COMPACT_TOP = 2  # ASSUMPTION (plan): of the top-20 no-fit rows


def decide(
    nofit_nuse: np.ndarray,
    fit_nuse: np.ndarray,
    nofit_r: np.ndarray,
    fit_r: np.ndarray,
    n_compact_top: int | None,
) -> str:
    """Frozen plan rule. ``n_compact_top`` is None when compactness is undetermined."""
    cut = np.percentile(fit_nuse, COV_PCT)
    if np.mean(nofit_nuse < cut) >= COV_FRACTION:
        return "H-cov: coverage flag"
    if (
        n_compact_top is not None
        and np.median(nofit_r) < np.percentile(fit_r, SIZE_PCT)
        and n_compact_top >= MIN_COMPACT_TOP
    ):
        return "H-star supported"
    return "H-mix / undetermined"


def match(src: Table, uids: list[str], cat: SkyCoord) -> tuple[np.ndarray, np.ndarray]:
    """DJA row index of each uid and whether it lies within ZOUT_MATCH_ARCSEC."""
    row = {str(x): i for i, x in enumerate(src["source_uid"])}
    pos = np.array([row[x] for x in uids])
    c = SkyCoord(np.asarray(src["ra"], float)[pos] * u.deg, np.asarray(src["dec"], float)[pos] * u.deg)
    idx, sep, _ = c.match_to_catalog_sky(cat)
    return idx, sep.arcsec <= ZOUT_MATCH_ARCSEC


def q(x: np.ndarray) -> dict[str, float | int]:
    x = x[np.isfinite(x)]
    if not x.size:
        return {"n": 0}
    p = np.percentile(x, [5, 10, 50, 90])
    return {"n": int(x.size)} | {f"p{k}": round(float(v), 3) for k, v in zip((5, 10, 50, 90), p, strict=True)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--zout", type=Path, required=True)
    ap.add_argument("--phot", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if file_sha256(args.zout) != ZOUT_SHA256 or file_sha256(args.phot) != PHOT_SHA256:
        raise SystemExit("zout or phot is not the pinned DJA v7.4 SMACS file")
    z, p = Table.read(args.zout), Table.read(args.phot)
    if not np.array_equal(np.asarray(z["id"]), np.asarray(p["id"])):
        raise SystemExit("zout and phot rows are not aligned")
    cat = SkyCoord(np.asarray(z["ra"], float) * u.deg, np.asarray(z["dec"], float) * u.deg)
    zp = np.asarray(z["z_phot"], float)
    nuse = np.asarray(z["nusefilt"], int)
    rad = np.asarray(p["flux_radius"], float)

    d = args.run_dir / SAMPLE
    src = Table.read(d / "sources.ecsv")
    scores = Table.read(d / "scores.ecsv")
    order = np.argsort(np.asarray(scores["rank"]))
    uids = [str(x) for x in np.asarray(scores["source_uid"])[order]]
    idx, ok = match(src, uids, cat)
    nofit = ok & (zp[idx] == -1)
    fit = ok & (zp[idx] >= 0)

    stars = [str(x) for x in Table.read(args.run_dir / f"{SAMPLE}-stars" / "scores.ecsv")["source_uid"]]
    s_idx, s_ok = match(src, stars, cat)
    star_r = rad[s_idx[s_ok]]
    compact_cut = float(np.percentile(star_r, STAR_PCT)) if star_r.size >= MIN_STARS else None

    top = []
    for k in range(TOP_K):
        if not nofit[k]:
            continue
        r = float(rad[idx[k]])
        top.append(
            {
                "rank": k + 1,
                "source_uid": uids[k],
                "dja_id": int(z["id"][idx[k]]),
                "nusefilt": int(nuse[idx[k]]),
                "flux_radius_pix": round(r, 3),
                "compact": None if compact_cut is None else bool(r <= compact_cut),
            }
        )
    n_compact = None if compact_cut is None else sum(t["compact"] for t in top)
    verdict = decide(nuse[idx[nofit]], nuse[idx[fit]], rad[idx[nofit]], rad[idx[fit]], n_compact)

    summary = {
        "run": args.run_dir.name,
        "zout_sha256": ZOUT_SHA256,
        "phot_sha256": PHOT_SHA256,
        "n_ranked": len(uids),
        "n_matched": int(ok.sum()),
        "n_nofit": int(nofit.sum()),
        "n_fit": int(fit.sum()),
        "nusefilt": {
            "nofit": q(nuse[idx[nofit]].astype(float)),
            "fit": q(nuse[idx[fit]].astype(float)),
            "nofit_below_fit_p5": round(float(np.mean(nuse[idx[nofit]] < np.percentile(nuse[idx[fit]], COV_PCT))), 4),
            "nofit_counts": {
                str(int(v)): int(c) for v, c in zip(*np.unique(nuse[idx[nofit]], return_counts=True), strict=True)
            },
        },
        "flux_radius_pix": {
            "nofit": q(rad[idx[nofit]]),
            "fit": q(rad[idx[fit]]),
            "stars": q(star_r),
            "n_stars_matched": int(s_ok.sum()),
            "compact_cut": compact_cut,
            "nofit_compact_fraction": None
            if compact_cut is None
            else round(float(np.mean(rad[idx[nofit]] <= compact_cut)), 4),
            "fit_compact_fraction": None
            if compact_cut is None
            else round(float(np.mean(rad[idx[fit]] <= compact_cut)), 4),
        },
        "top20_nofit": top,
        "n_compact_top20_nofit": n_compact,
        "verdict": verdict,
        "assumptions": {
            "match_arcsec": ZOUT_MATCH_ARCSEC,
            "cov_fraction": COV_FRACTION,
            "cov_pct": COV_PCT,
            "size_pct": SIZE_PCT,
            "star_pct": STAR_PCT,
            "min_stars": MIN_STARS,
            "min_compact_top": MIN_COMPACT_TOP,
            "top_k": TOP_K,
        },
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: summary[k] for k in ("n_nofit", "n_fit", "nusefilt", "flux_radius_pix", "verdict")}))
    print(json.dumps(top))


if __name__ == "__main__":
    main()
