"""SMACS galaxy ranking scored within DJA z_phot strata (members vs background) vs the pipeline.

Plan: ``results/smacs_jwst3/zphot_stratify/PLAN.md`` (decision rule frozen there before a re-score).
Reads one pipeline run of ``configs/smacs_jwst3_gaia_only.yaml``, re-scores the saved galaxy stratum
(``rank.score_anomalies``, random_state 0; the baseline must reproduce the pipeline ranks exactly),
then scores the member (z_phot 0.25-0.55) and background (z_phot >= 0.8) strata separately. Writes
``summary.json``.

Provenance: features and z_phot are ``derived`` (from ``observed`` MAST/DJA catalogs); ICLv2 image
positions are ``observed`` (Mahler+2022 catalogue); every statistic here is ``derived``; the
ranges and thresholds are ASSUMPTIONs from the plan.

Usage: python scripts/smacs_zphot_stratify.py <run_dir> --zout <zout.fits> --out <dir>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402
from science_control_features import Z_BACKGROUND_MIN, Z_MEMBER, zphot_for  # noqa: E402

from jwst_anomaly.rank import score_anomalies  # noqa: E402

SAMPLE = "smacs0723_nircam"
TOP_K = 20  # ASSUMPTION (plan)
IMAGE_MATCH_ARCSEC = 0.3  # ASSUMPTION (plan)
ENRICH_FACTOR = 2.0  # ASSUMPTION (plan): member share of the top k vs of the stratum
ENRICH_MIN = 4  # ASSUMPTION (plan): members in the top k


def stratum_of(z: np.ndarray) -> np.ndarray:
    """``member`` / ``background`` / ``other`` (NaN or outside both ranges)."""
    out = np.full(len(z), "other", dtype=object)
    with np.errstate(invalid="ignore"):
        out[(z >= Z_MEMBER[0]) & (z <= Z_MEMBER[1])] = "member"
        out[z >= Z_BACKGROUND_MIN] = "background"
    return out


def ranks_of(feats: Table) -> dict[str, int]:
    s = score_anomalies(feats, random_state=0)
    return {str(u): int(r) for u, r in zip(s["source_uid"], s["rank"], strict=True)}


def decide(n_top_member: int, frac_top: float, frac_stratum: float) -> str:
    """Frozen plan rule."""
    if n_top_member >= ENRICH_MIN and frac_top >= ENRICH_FACTOR * frac_stratum:
        return "population effect confirmed (H0): stratify SMACS top-k statements by z_phot"
    return "no rule change"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--zout", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    d = args.run_dir / SAMPLE
    feats = Table.read(d / "features.ecsv")
    scores = Table.read(d / "scores.ecsv")
    pipeline = {str(u): int(r) for u, r in zip(scores["source_uid"], scores["rank"], strict=True)}
    uids = np.array([str(u) for u in feats["source_uid"]])
    ranked = feats[np.isin(uids, list(pipeline))]
    r_uids = np.array([str(u) for u in ranked["source_uid"]])

    base = ranks_of(ranked)
    n_diff = sum(base[u] != pipeline[u] for u in pipeline)
    if n_diff:
        raise SystemExit(f"baseline re-score differs from the pipeline: {n_diff} ranks")

    z = zphot_for(args.run_dir, SAMPLE, r_uids, args.zout)
    strat = stratum_of(z)
    per = {s: ranks_of(ranked[strat == s]) for s in ("member", "background")}
    n = len(r_uids)
    counts = {s: int((strat == s).sum()) for s in ("member", "background", "other")}

    order = sorted(r_uids, key=lambda u: base[u])
    top = order[:TOP_K]
    s_of = dict(zip(r_uids, strat, strict=True))
    z_of = dict(zip(r_uids, z, strict=True))
    top_counts = {s: sum(s_of[u] == s for u in top) for s in counts}
    frac_stratum = counts["member"] / n
    frac_top = top_counts["member"] / TOP_K

    # Overlap: each stratum's own top k_s (its share of TOP_K, at least 1) vs the overall top 20.
    overlap = {}
    for s in ("member", "background"):
        k_s = max(1, round(TOP_K * counts[s] / n))
        own = {u for u, r in per[s].items() if r <= k_s}
        overlap[s] = {
            "k_s": k_s,
            "own_top_in_overall_top20": len(own & set(top)),
            "own_top_overall_ranks": sorted(base[u] for u in own),
        }
    bg_top20 = sorted(per["background"], key=lambda u: per["background"][u])[:TOP_K]

    # Descriptive: ICLv2 catalogued multiple images.
    _, files, par = lc.load_model("smacs0723-iclv2")
    images, _ = lc.image_list("smacs0723-iclv2", files, par)
    src = Table.read(d / "sources.ecsv")
    srow = {str(u): i for i, u in enumerate(src["source_uid"])}
    pos = np.array([srow[u] for u in r_uids])
    cs = SkyCoord(np.asarray(src["ra"], float)[pos], np.asarray(src["dec"], float)[pos], unit="deg")
    ci = SkyCoord(images["ra"], images["dec"], unit="deg")
    j, sep, _ = ci.match_to_catalog_sky(cs)
    img_rows = []
    for k in np.flatnonzero(sep.arcsec <= IMAGE_MATCH_ARCSEC):
        u = r_uids[j[k]]
        img_rows.append(
            {
                "image": str(images["image_id"][k]),
                "source_uid": u,
                "sep_arcsec": round(float(sep.arcsec[k]), 3),
                "stratum": s_of[u],
                "overall_rank": base[u],
                "background_rank": per["background"].get(u),
            }
        )
    bg_imgs = [r for r in img_rows if r["background_rank"] is not None]
    med_overall = float(np.median([r["overall_rank"] for r in bg_imgs])) if bg_imgs else None
    med_bg = float(np.median([r["background_rank"] for r in bg_imgs])) if bg_imgs else None
    # Background ranks live on a smaller list; compare percentiles too.
    nb = counts["background"]

    # Post hoc (deviation, added after the first pass): the "other" stratum split into eazy's
    # z_phot = -1 (no fit), no zout match within 0.2", and finite z_phot outside both ranges.
    with np.errstate(invalid="ignore"):
        other_kind = np.where(~np.isfinite(z), "unmatched", np.where(z < 0, "zphot_-1", "between"))
    other_kind = np.where(strat == "other", other_kind, "")
    kind_of = dict(zip(r_uids, other_kind, strict=True))
    post_hoc_other = {
        k: {"ranked": int((other_kind == k).sum()), "top20": sum(kind_of[u] == k for u in top)}
        for k in ("zphot_-1", "unmatched", "between")
    }

    screened = {}
    if (d / "screened.ecsv").exists():  # D-019-D-021 spike / weight screen of the top rows
        t = Table.read(d / "screened.ecsv")
        screened = {str(u): str(r) for u, r in zip(t["source_uid"], t["reason"], strict=True)}

    summary = {
        "run": args.run_dir.name,
        "post_hoc_other_split": post_hoc_other,
        "n_ranked": n,
        "strata": counts,
        "top20_strata": top_counts,
        "member_fraction": {"stratum": round(frac_stratum, 4), "top20": frac_top},
        "verdict": decide(top_counts["member"], frac_top, frac_stratum),
        "overlap": overlap,
        "top20": [
            {
                "rank": base[u],
                "source_uid": u,
                "stratum": s_of[u],
                "z_phot": None if not np.isfinite(z_of[u]) else round(float(z_of[u]), 3),
                "background_rank": per["background"].get(u),
                "member_rank": per["member"].get(u),
                "other_kind": kind_of[u] or None,
                "screened": screened.get(u),
            }
            for u in top
        ],
        "background_top20_overall_ranks": [base[u] for u in bg_top20],
        "iclv2_images": {
            "n_catalogued": len(images),
            "n_matched_ranked": len(img_rows),
            "n_matched_background": len(bg_imgs),
            "median_overall_rank_of_background_images": med_overall,
            "median_overall_percentile": None if med_overall is None else round(med_overall / n, 4),
            "median_background_rank": med_bg,
            "median_background_percentile": None if med_bg is None else round(med_bg / nb, 4),
            "rows": img_rows,
        },
        "assumptions": {
            "z_member": Z_MEMBER,
            "z_background_min": Z_BACKGROUND_MIN,
            "top_k": TOP_K,
            "image_match_arcsec": IMAGE_MATCH_ARCSEC,
            "enrich_factor": ENRICH_FACTOR,
            "enrich_min": ENRICH_MIN,
        },
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(
        json.dumps(
            {
                k: summary[k]
                for k in ("strata", "top20_strata", "member_fraction", "verdict", "overlap")
            }
        )
    )
    print(json.dumps({k: v for k, v in summary["iclv2_images"].items() if k != "rows"}))


if __name__ == "__main__":
    main()
