"""SMACS point-source screen: DJA ``flux_radius`` against the stellar locus (TASKS Next M2).

Plan: ``results/smacs_jwst3/point_source/PLAN.md`` (rule and validation frozen before the run).
Reads one run of ``configs/smacs_jwst3_gaia_only.yaml`` and the pinned DJA v7.4 SMACS phot catalog;
derives the locus from Gaia-split star-stratum sources, validates it on held-out stars (V1) and on
Mahler+2022 ICLv2 catalogued multiple images (V2), then labels the plan's two targets.

Provenance: catalog columns are ``observed`` / ``derived`` (DJA, MAST, ICLv2 ``model_prediction``
positions of catalogued images are observed identifications); every statistic is ``derived``; the
thresholds are ASSUMPTIONs from the plan.

``--plan v2`` applies ``results/smacs_jwst3/point_source_v2/PLAN.md``: locus window [20, 24], no
fallback, plus a spatial held-out check (V1b) and a faint-star extrapolation check (V3).

Usage: python scripts/smacs_point_source.py <run_dir> --phot <phot.fits> --out <dir> [--plan v2]
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
from smacs_zphot_nofit import PHOT_SHA256, SAMPLE  # noqa: E402
from smacs_zphot_nofit import match as uid_match  # noqa: E402

from jwst_anomaly.provenance import file_sha256  # noqa: E402

LOCUS_MAG = (19.0, 24.0)  # ASSUMPTION (plan): brighter stars saturate in F200W
LOCUS_MAG_FALLBACK = (18.0, 24.0)  # ASSUMPTION (plan)
MIN_LOCUS = 10  # ASSUMPTION (plan)
N_SIGMA = 3.0  # ASSUMPTION (plan)
V1_MIN_RECOVERY = 0.90  # ASSUMPTION (plan)
V2_MIN_IMAGES = 10  # ASSUMPTION (plan)
V2_MAX_POINT = 0.10  # ASSUMPTION (plan)
TARGETS = (
    "jw02736-o001_t001_nircam_f200w_1816",
    "jw02736-o001_t001_nircam_f200w_384",
)
MAG_EDGES = (18.0, 22.0, 24.0, 26.0, 27.0, 28.0, 31.0)  # descriptive bins
TOP_K = 20
# v2 plan (results/smacs_jwst3/point_source_v2/PLAN.md), frozen before the v2 run
LOCUS_MAG_V2 = (
    20.0,
    24.0,
)  # ASSUMPTION (v2 plan; chosen after v1's outliers, so V1a is not independent)
V3_MAG = (24.0, 27.0)  # ASSUMPTION (v2 plan): faint extrapolation range, (lo, hi]
V3_MIN_STARS = 10  # ASSUMPTION (v2 plan): fewer → V3 descriptive
V3_MIN_POINT = 0.80  # ASSUMPTION (v2 plan)


def locus_cut(radius: np.ndarray) -> tuple[float, float, float]:
    """``(r_psf, s, cut)``: median, 1.4826 MAD and the point-like upper bound."""
    r = radius[np.isfinite(radius)]
    med = float(np.median(r))
    s = float(1.4826 * np.median(np.abs(r - med)))
    return med, s, med + N_SIGMA * s


def locus_mask(mag: np.ndarray, radius: np.ndarray) -> tuple[np.ndarray, list[float] | None]:
    """Locus stars in the plan window (fallback window if too few); None = undetermined."""
    for lo, hi in (LOCUS_MAG, LOCUS_MAG_FALLBACK):
        m = (mag >= lo) & (mag <= hi) & np.isfinite(radius)
        if m.sum() >= MIN_LOCUS:
            return m, [lo, hi]
    return np.zeros_like(mag, bool), None


def held_out_recovery(ids: np.ndarray, radius: np.ndarray) -> list[float]:
    """V1: cut from one ``id``-parity half, recovery fraction on the other, both directions."""
    out = []
    for train in (0, 1):
        a, b = ids % 2 == train, ids % 2 != train
        if a.sum() < 2 or b.sum() < 1:
            return []
        cut = locus_cut(radius[a])[2]
        out.append(float(np.mean(radius[b] <= cut)))
    return out


def held_out_by_group(group: np.ndarray, radius: np.ndarray) -> list[float]:
    """Cut from one group (``group`` True, then False), recovery fraction on the other; [] if a
    group is too small to split."""
    out = []
    for train in (True, False):
        a, b = group == train, group != train
        if a.sum() < 2 or b.sum() < 1:
            return []
        cut = locus_cut(radius[a])[2]
        out.append(float(np.mean(radius[b] <= cut)))
    return out


def spatial_groups(ra: np.ndarray) -> np.ndarray:
    """V1b halves: True = east (``ra`` above the locus median ``ra``)."""
    return ra > np.median(ra)


def decide_v2(
    rec_parity: list[float],
    rec_spatial: list[float],
    n_images: int,
    image_point_frac: float | None,
    n_faint: int,
    faint_point_frac: float | None,
) -> str:
    """Frozen v2 plan rule: V1a, V1b, V2 combined with "and"; V3 sets the magnitude reach."""
    if len(rec_parity) != 2 or len(rec_spatial) != 2:
        return "undetermined (V1 could not be split)"
    fails = []
    if min(rec_parity) < V1_MIN_RECOVERY:
        fails.append("V1a")
    if min(rec_spatial) < V1_MIN_RECOVERY:
        fails.append("V1b")
    if not (
        n_images >= V2_MIN_IMAGES
        and image_point_frac is not None
        and image_point_frac <= V2_MAX_POINT
    ):
        fails.append("V2")
    if fails:
        return "not validated (" + ", ".join(fails) + " fail)"
    v3 = (
        n_faint >= V3_MIN_STARS
        and faint_point_frac is not None
        and faint_point_frac >= V3_MIN_POINT
    )
    return "validated to ref_mag 27" if v3 else "validated to ref_mag 24 only"


def decide(recovery: list[float], n_images: int, image_point_frac: float | None) -> str:
    """Frozen plan rule for the validation gates. An empty ``recovery`` (a parity half too small to
    split) means V1 could not be run: undetermined, not a failure."""
    if len(recovery) != 2:
        return "undetermined (V1 could not be split)"
    v1 = min(recovery) >= V1_MIN_RECOVERY
    v2 = (
        n_images >= V2_MIN_IMAGES
        and image_point_frac is not None
        and image_point_frac <= V2_MAX_POINT
    )
    if v1 and v2:
        return "validated"
    return "not validated" + ("" if v1 else " (V1 fails)") + ("" if v2 else " (V2 fails)")


def match(ra: np.ndarray, dec: np.ndarray, cat: SkyCoord) -> tuple[np.ndarray, np.ndarray]:
    idx, sep, _ = SkyCoord(ra * u.deg, dec * u.deg).match_to_catalog_sky(cat)
    return idx, sep.arcsec <= ZOUT_MATCH_ARCSEC


def _stratum(run_dir: Path, name: str, src: Table, cat: SkyCoord, fmag: dict) -> Table:
    """Ranked uids of a stratum (rank order) with DJA index, match flag and ref_mag."""
    sc = Table.read(run_dir / name / "scores.ecsv")
    sc = sc[np.argsort(np.asarray(sc["rank"]))]
    uids = [str(x) for x in sc["source_uid"]]
    idx, ok = uid_match(src, uids, cat)
    return Table(
        {
            "source_uid": uids,
            "dja": idx,
            "ok": ok,
            "ref_mag": [fmag.get(x, np.nan) for x in uids],
        }
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--phot", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--plan", choices=("v1", "v2"), default="v1", help="frozen plan to apply")
    args = ap.parse_args()
    v2 = args.plan == "v2"
    if file_sha256(args.phot) != PHOT_SHA256:
        raise SystemExit("phot is not the pinned DJA v7.4 SMACS file")
    p = Table.read(args.phot)
    cat = SkyCoord(np.asarray(p["ra"], float) * u.deg, np.asarray(p["dec"], float) * u.deg)
    rad = np.asarray(p["flux_radius"], float)
    ids = np.asarray(p["id"], int)

    src = Table.read(args.run_dir / SAMPLE / "sources.ecsv")
    feats = Table.read(args.run_dir / SAMPLE / "features.ecsv")  # shared by both strata
    fmag = {str(x): float(m) for x, m in zip(feats["source_uid"], feats["ref_mag"], strict=True)}
    stars = _stratum(args.run_dir, f"{SAMPLE}-stars", src, cat, fmag)
    gal = _stratum(args.run_dir, SAMPLE, src, cat, fmag)

    s = stars[stars["ok"]]
    s_rad, s_mag, s_ids = rad[s["dja"]], np.asarray(s["ref_mag"], float), ids[s["dja"]]
    s_ra = np.asarray(p["ra"], float)[s["dja"]]
    if v2:  # no fallback window
        lo, hi = LOCUS_MAG_V2
        lm = (s_mag >= lo) & (s_mag <= hi) & np.isfinite(s_rad)
        window = [lo, hi] if lm.sum() >= MIN_LOCUS else None
    else:
        lm, window = locus_mask(s_mag, s_rad)
    summary: dict = {
        "plan": args.plan,
        "run": args.run_dir.name,
        "phot_sha256": PHOT_SHA256,
        "n_stars_ranked": len(stars),
        "n_stars_matched": int(stars["ok"].sum()),
        "locus_window": window,
        "n_locus": int(lm.sum()),
    }
    if window is None:
        summary["verdict"] = "undetermined (too few locus stars)"
    else:
        r_psf, sig, cut = locus_cut(s_rad[lm])
        rec = held_out_recovery(s_ids[lm], s_rad[lm])

        # lazy: loads the ICLv2 files from the model cache (lens_consistency.MODELS)
        import lens_consistency as lc

        _, files, par = lc.load_model("smacs0723-iclv2")
        images, _ = lc.image_list("smacs0723-iclv2", files, par)
        i_idx, i_ok = match(np.asarray(images["ra"], float), np.asarray(images["dec"], float), cat)
        i_rad = rad[i_idx[i_ok]]
        i_fin = np.isfinite(i_rad)  # a missing radius is neither kept nor removed
        i_ok[np.flatnonzero(i_ok)[~i_fin]] = False
        i_rad = i_rad[i_fin]
        i_point = i_rad <= cut
        n_img = int(i_ok.sum())
        frac = float(np.mean(i_point)) if n_img else None
        if v2:
            rec_sp = held_out_by_group(spatial_groups(s_ra[lm]), s_rad[lm])
            fm = (s_mag > V3_MAG[0]) & (s_mag <= V3_MAG[1]) & np.isfinite(s_rad)
            n_faint = int(fm.sum())
            faint_frac = float(np.mean(s_rad[fm] <= cut)) if n_faint else None
            verdict = decide_v2(rec, rec_sp, n_img, frac, n_faint, faint_frac)
            summary |= {
                "v1b_recovery_spatial": [round(x, 4) for x in rec_sp],
                "v1b_split_ra_deg": round(float(np.median(s_ra[lm])), 6),
                "v3_n_faint_stars": n_faint,
                "v3_point_fraction": None if faint_frac is None else round(faint_frac, 4),
                "v3_is_gate": n_faint >= V3_MIN_STARS,
                # Post hoc (deviation, added after V1a and V1b failed): locus radius by ref_mag and
                # where the outliers sit in the V1b halves.
                "post_hoc_locus_by_ref_mag": [
                    {
                        "ref_mag": [lo, hi],
                        "n": int(b.sum()),
                        "median_radius_pix": round(float(np.median(s_rad[b])), 3)
                        if b.any()
                        else None,
                        "above_cut": int(np.sum(s_rad[b] > cut)),
                    }
                    for lo, hi in ((20.0, 20.7), (20.7, 22.0), (22.0, 24.0))
                    for b in [lm & (s_mag >= lo) & (s_mag < hi if hi < 24.0 else s_mag <= hi)]
                ],
                "post_hoc_outliers_east_west": [
                    int(np.sum((s_rad[lm] > cut) & spatial_groups(s_ra[lm]))),
                    int(np.sum((s_rad[lm] > cut) & ~spatial_groups(s_ra[lm]))),
                ],
                "v3_faint_stars": [
                    {"ref_mag": round(float(m), 3), "flux_radius_pix": round(float(r), 3)}
                    for m, r in sorted(zip(s_mag[fm], s_rad[fm], strict=True))
                ],
            }
        else:
            verdict = decide(rec, n_img, frac)

        g_ok = np.asarray(gal["ok"])
        g_rad = np.where(g_ok, rad[gal["dja"]], np.nan)
        g_mag = np.asarray(gal["ref_mag"], float)
        bins = []
        for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:], strict=True):
            b = g_ok & (g_mag > lo) & (g_mag <= hi)
            bins.append(
                {
                    "ref_mag": [lo, hi],
                    "n": int(b.sum()),
                    "point_like": int(np.sum(g_rad[b] <= cut)),
                }
            )

        def top_row(k: int) -> dict:
            fin = bool(g_ok[k] and np.isfinite(g_rad[k]))
            return {
                "rank": k + 1,
                "source_uid": str(gal["source_uid"][k]),
                "ref_mag": round(float(g_mag[k]), 3) if np.isfinite(g_mag[k]) else None,
                "flux_radius_pix": round(float(g_rad[k]), 3) if fin else None,
                "point_like": bool(g_rad[k] <= cut) if fin else None,
            }

        top = [top_row(k) for k in range(min(TOP_K, len(gal)))]
        pos = {x: k for k, x in enumerate(gal["source_uid"])}
        targets = {t: top_row(pos[t]) if t in pos else None for t in TARGETS}
        if v2:  # v2 plan: targets (ref_mag > 24) are labelled only with the 27 mag reach
            summary["targets_labelled"] = verdict == "validated to ref_mag 27"
        point_images = [
            {
                "id": str(images["id"][j]) if "id" in images.colnames else int(j),
                "flux_radius_pix": round(float(r), 3),
            }
            for j, r in zip(np.flatnonzero(i_ok)[i_point], i_rad[i_point], strict=True)
        ]
        # Post hoc (deviation, added after V1 failed): which locus stars fall above the cut.
        out_m = s_rad[lm] > cut
        locus_outliers = [
            {"ref_mag": round(float(m), 3), "flux_radius_pix": round(float(r), 3)}
            for m, r in sorted(zip(s_mag[lm][out_m], s_rad[lm][out_m], strict=True))
        ]
        # Review check: several ranked uids (per-band detections) can match one DJA object.
        l_ids = s_ids[lm]
        _, first = np.unique(l_ids, return_index=True)
        dedup = np.zeros(l_ids.size, bool)
        dedup[first] = True
        rec_dedup = held_out_recovery(l_ids[dedup], s_rad[lm][dedup])
        binned = sum(b["n"] for b in bins)
        summary |= {
            "post_hoc_locus_outliers": locus_outliers,
            "locus_duplicate_dja_ids": int(l_ids.size - dedup.sum()),
            "review_v1_recovery_dedup": [round(x, 4) for x in rec_dedup],
            "galaxy_matched": int(g_ok.sum()),
            "galaxy_matched_outside_bins_or_nan_mag": int(g_ok.sum() - binned),
            "v2_images_nan_radius": int((~i_fin).sum()),
            "r_psf": round(r_psf, 4),
            "mad_sigma": round(sig, 4),
            "cut": round(cut, 4),
            "locus_radius_quantiles": [
                round(float(x), 3) for x in np.percentile(s_rad[lm], [5, 50, 95])
            ],
            "stars_all_point_fraction": round(float(np.mean(s_rad <= cut)), 4),
            "v1_recovery": [round(x, 4) for x in rec],
            "v2_n_images": n_img,
            "v2_n_images_listed": len(images),
            "v2_point_fraction": None if frac is None else round(frac, 4),
            "v2_point_images": point_images,
            "image_radius_quantiles": [
                round(float(x), 3) for x in np.percentile(i_rad, [5, 50, 95])
            ]
            if n_img
            else None,
            "verdict": verdict,
            "galaxy_point_like_by_ref_mag": bins,
            "top20": top,
            "targets": targets,
        }
    summary["assumptions"] = {
        "match_arcsec": ZOUT_MATCH_ARCSEC,
        "locus_mag": LOCUS_MAG_V2 if v2 else LOCUS_MAG,
        "locus_mag_fallback": None if v2 else LOCUS_MAG_FALLBACK,
        "min_locus": MIN_LOCUS,
        "n_sigma": N_SIGMA,
        "v1_min_recovery": V1_MIN_RECOVERY,
        "v2_min_images": V2_MIN_IMAGES,
        "v2_max_point": V2_MAX_POINT,
    }
    if v2:
        summary["assumptions"] |= {
            "v3_mag": V3_MAG,
            "v3_min_stars": V3_MIN_STARS,
            "v3_min_point": V3_MIN_POINT,
        }
    args.out.mkdir(parents=True, exist_ok=True)
    text = json.dumps(summary, indent=1, allow_nan=False)  # strict JSON: no bare NaN
    (args.out / "summary.json").write_text(text + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k not in ("top20",)}, indent=1))


if __name__ == "__main__":
    main()
