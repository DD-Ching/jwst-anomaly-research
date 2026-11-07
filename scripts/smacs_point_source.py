"""SMACS point-source screen: DJA ``flux_radius`` against the stellar locus (TASKS Next M2).

Plan: ``results/smacs_jwst3/point_source/PLAN.md`` (rule and validation frozen before the run).
Reads one run of ``configs/smacs_jwst3_gaia_only.yaml`` and the pinned DJA v7.4 SMACS phot catalog;
derives the locus from Gaia-split star-stratum sources, validates it on held-out stars (V1) and on
Mahler+2022 ICLv2 catalogued multiple images (V2), then labels the plan's two targets.

Provenance: catalog columns are ``observed`` / ``derived`` (DJA, MAST, ICLv2 ``model_prediction``
positions of catalogued images are observed identifications); every statistic is ``derived``; the
thresholds are ASSUMPTIONs from the plan.

Usage: python scripts/smacs_point_source.py <run_dir> --phot <phot.fits> --out <dir>
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


def decide(recovery: list[float], n_images: int, image_point_frac: float | None) -> str:
    """Frozen plan rule for the validation gates."""
    v1 = len(recovery) == 2 and min(recovery) >= V1_MIN_RECOVERY
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


def _stratum(run_dir: Path, name: str, src: Table, cat: SkyCoord) -> Table:
    """Ranked uids of a stratum (rank order) with DJA index, match flag and ref_mag."""
    d = run_dir / name
    sc = Table.read(d / "scores.ecsv")
    sc = sc[np.argsort(np.asarray(sc["rank"]))]
    uids = [str(x) for x in sc["source_uid"]]
    row = {str(x): i for i, x in enumerate(src["source_uid"])}
    pos = np.array([row.get(x, -1) for x in uids])
    have = pos >= 0
    ra = np.where(have, np.asarray(src["ra"], float)[pos], 0.0)
    dec = np.where(have, np.asarray(src["dec"], float)[pos], 0.0)
    idx, ok = match(ra, dec, cat)
    feats = Table.read(run_dir / SAMPLE / "features.ecsv")  # shared by both strata
    fmag = {str(x): float(m) for x, m in zip(feats["source_uid"], feats["ref_mag"], strict=True)}
    return Table(
        {
            "source_uid": uids,
            "dja": idx,
            "ok": have & ok,
            "ref_mag": [fmag.get(x, np.nan) for x in uids],
        }
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--phot", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if file_sha256(args.phot) != PHOT_SHA256:
        raise SystemExit("phot is not the pinned DJA v7.4 SMACS file")
    p = Table.read(args.phot)
    cat = SkyCoord(np.asarray(p["ra"], float) * u.deg, np.asarray(p["dec"], float) * u.deg)
    rad = np.asarray(p["flux_radius"], float)
    ids = np.asarray(p["id"], int)

    src = Table.read(args.run_dir / SAMPLE / "sources.ecsv")
    stars = _stratum(args.run_dir, f"{SAMPLE}-stars", src, cat)
    gal = _stratum(args.run_dir, SAMPLE, src, cat)

    s = stars[stars["ok"]]
    s_rad, s_mag, s_ids = rad[s["dja"]], np.asarray(s["ref_mag"], float), ids[s["dja"]]
    lm, window = locus_mask(s_mag, s_rad)
    summary: dict = {
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
        i_point = i_rad <= cut
        n_img = int(i_ok.sum())
        frac = float(np.mean(i_point)) if n_img else None
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
        top = [
            {
                "rank": k + 1,
                "source_uid": gal["source_uid"][k],
                "ref_mag": round(float(g_mag[k]), 3),
                "flux_radius_pix": None if not g_ok[k] else round(float(g_rad[k]), 3),
                "point_like": None if not g_ok[k] else bool(g_rad[k] <= cut),
            }
            for k in range(min(TOP_K, len(gal)))
        ]
        tmap = {t["source_uid"]: t for t in top}
        targets = {}
        for t in TARGETS:
            row = tmap.get(t)
            if row is None:
                k = list(gal["source_uid"]).index(t) if t in list(gal["source_uid"]) else None
                row = (
                    None
                    if k is None
                    else {
                        "rank": k + 1,
                        "flux_radius_pix": round(float(g_rad[k]), 3) if g_ok[k] else None,
                        "point_like": bool(g_rad[k] <= cut) if g_ok[k] else None,
                    }
                )
            targets[t] = row
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
        summary |= {
            "post_hoc_locus_outliers": locus_outliers,
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
        "locus_mag": LOCUS_MAG,
        "locus_mag_fallback": LOCUS_MAG_FALLBACK,
        "min_locus": MIN_LOCUS,
        "n_sigma": N_SIGMA,
        "v1_min_recovery": V1_MIN_RECOVERY,
        "v2_min_images": V2_MIN_IMAGES,
        "v2_max_point": V2_MAX_POINT,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k not in ("top20",)}, indent=1))


if __name__ == "__main__":
    main()
