"""Two-epoch check of SMACS ``f200w_1816`` and ``f200w_384`` (2736 o001 -> VENUS 6882 o057).

Plan: ``results/smacs_jwst3/epochs_1816_384/PLAN.md`` (frozen in the claim commit). ``decide``
encodes that plan's outcome rules. The empirical null (``empirical_null``) is post hoc: it compares
the target's F150W residual with mutual matches of similar magnitude and size, because the D-017
error model (tie + ``semimajor_sigma`` / S/N) is not calibrated for faint compact sources.

    python scripts/smacs_epochs_1816_384.py --cat-dir $JWST_ANOMALY_DATA/cache/s3cat \\
        --out results/smacs_jwst3/epochs_1816_384/summary.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

from jwst_anomaly import catalog

sys.path.insert(0, str(Path(__file__).resolve().parent))
from epoch_compare import compare_epochs, mutual_matches, offsets_mas  # noqa: E402

TARGETS = {"f200w_1816": 1816, "f200w_384": 384}
BANDS = ("f150w", "f444w")
E1 = "jw02736-o001_t001_nircam_clear-{b}_cat.ecsv"
E2 = "jw06882-o057_t057_nircam_clear-{b}_cat.ecsv"
# Frozen plan thresholds (ASSUMPTIONs).
MOTION_SIG = 3.0
NO_MOTION_SIG = 2.0
NO_MOTION_BOUND_MAS = 15.0
ANGLE_DEG = 45.0
VAR_FACTOR = 5.0
# Post hoc null selection (ASSUMPTIONs).
NULL_DMAG = 0.5
NULL_SIZE_FRAC = 0.3
NULL_RADIUS_ARCSEC = 60.0
# Epoch baseline 2022-06-07 -> 2026-06-05 (D-017 evidence), in years.
BASELINE_YR = 4.0


def _angle_deg(a: list[float], b: list[float]) -> float:
    cos = np.dot(a, b) / (np.hypot(*a) * np.hypot(*b))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def decide(per_band: dict[str, dict[str, Any] | None]) -> dict[str, str]:
    """Plan outcome rules. ``per_band[b]`` is a ``compare_epochs`` result, or None if not measured.

    Motion: every band in BANDS measured, each with significance >= MOTION_SIG and residual
    directions within ANGLE_DEG -> "moves: star"; every band measured with significance <
    NO_MOTION_SIG and 3 x error <= NO_MOTION_BOUND_MAS -> "no motion detected"; else "undetermined".
    Variability: every band measured with |dmag| >= VAR_FACTOR x reference sigma -> "variable";
    every band measured and none so -> "not variable" (a deviation: the plan names only "variable");
    else "undetermined". A band with zero significance has no error estimate, so it cannot give
    "no motion detected".
    """
    res = [per_band.get(b) for b in BANDS]
    if any(r is None for r in res):
        return {"motion": "undetermined", "variability": "undetermined"}
    sig = [r["target_residual_significance"] for r in res]
    err = [
        r["target_residual_total_mas"] / s if s > 0 else np.inf
        for r, s in zip(res, sig, strict=True)
    ]
    if (
        all(s >= MOTION_SIG for s in sig)
        and _angle_deg(*(r["target_residual_mas"] for r in res)) <= ANGLE_DEG
    ):
        motion = "moves: star"
    elif all(s < NO_MOTION_SIG for s in sig) and all(3 * e <= NO_MOTION_BOUND_MAS for e in err):
        motion = "no motion detected"
    else:
        motion = "undetermined"
    dm = [r["target_dmag"] for r in res]
    if any(d is None for d in dm):
        var = "undetermined"
    else:
        hits = [
            abs(d) >= VAR_FACTOR * r["reference_dmag_sigma"] for d, r in zip(dm, res, strict=True)
        ]
        var = "variable" if all(hits) else ("not variable" if not any(hits) else "undetermined")
    return {"motion": motion, "variability": var}


def empirical_null(
    c1: SkyCoord, c2: SkyCoord, cat1: Table, k1: int, k2: int, frame_offset_mas: list[float]
) -> dict[str, Any]:
    """Post hoc: the target's residual against mutual matches of similar magnitude and size.

    Uses the same frame tie as ``compare_epochs`` (``frame_offset_mas``) and drops both target rows.
    """
    i1, i2 = mutual_matches(c1, c2, 0.3)
    mag = np.asarray(cat1["aper50_abmag"], float)
    size = np.ma.filled(np.ma.asarray(cat1["semimajor_sigma"], float), np.nan)
    near = c1[k1].separation(c1[i1]).arcsec <= NULL_RADIUS_ARCSEC
    med = np.asarray(frame_offset_mas, float)
    d_all = offsets_mas(c1[i1], c2[i2])
    peer = (
        near
        & (np.abs(mag[i1] - mag[k1]) <= NULL_DMAG)
        & (np.abs(size[i1] / size[k1] - 1.0) <= NULL_SIZE_FRAC)
        & (i1 != k1)
        & (i2 != k2)
        & np.isfinite(d_all).all(axis=1)
    )
    r_peer = np.hypot(*(d_all[peer] - med).T)
    r_t = float(np.hypot(*(offsets_mas(c1[[k1]], c2[[k2]])[0] - med)))
    return {
        "n_peers": int(peer.sum()),
        "target_residual_mas": r_t,
        "peer_residual_median_mas": float(np.median(r_peer)) if peer.any() else None,
        "peer_residual_p95_mas": float(np.percentile(r_peer, 95)) if peer.any() else None,
        "fraction_of_peers_with_larger_residual": float(np.mean(r_peer > r_t))
        if peer.any()
        else None,
        "n_peers_with_larger_residual": int(np.sum(r_peer > r_t)),
        "label": "derived (post hoc)",
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cat-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    ref = Table.read(args.cat_dir / E1.format(b="f200w"))
    out: dict[str, Any] = {"plan": "results/smacs_jwst3/epochs_1816_384/PLAN.md", "targets": {}}
    for name, label in TARGETS.items():
        row = ref[ref["label"] == label][0]
        target = SkyCoord(row["sky_centroid"].ra, row["sky_centroid"].dec)
        rec: dict[str, Any] = {
            "ra": float(target.ra.deg),
            "dec": float(target.dec.deg),
            "f200w_aper50_abmag": float(row["aper50_abmag"]),
            "bands": {},
        }
        per_band: dict[str, dict[str, Any] | None] = {}
        for b in BANDS:
            cat1 = catalog.load_pipeline_catalog(args.cat_dir / E1.format(b=b))
            cat2 = catalog.load_pipeline_catalog(args.cat_dir / E2.format(b=b))
            entry: dict[str, Any] = {}
            c1 = SkyCoord(cat1["ra"], cat1["dec"], unit="deg")
            c2 = SkyCoord(cat2["ra"], cat2["dec"], unit="deg")
            ks = []
            for ep, cat, c in (("epoch1", cat1, c1), ("epoch2", cat2, c2)):
                sep = target.separation(c).arcsec
                k = int(np.argmin(sep))
                ks.append(k)
                entry[f"{ep}_nearest_arcsec"] = float(sep[k])
                entry[f"{ep}_nearest_aper50_abmag"] = float(cat["aper50_abmag"][k])
            k1, k2 = ks  # the same rows compare_epochs uses as the target (nearest in each epoch)
            try:
                r = compare_epochs(cat1, cat2, target)
                entry["compare"] = r
                entry["target_aper50_abmag"] = [
                    float(cat1["aper50_abmag"][k1]),
                    float(cat2["aper50_abmag"][k2]),
                ]
                sig = r["target_residual_significance"]
                entry["derived"] = {
                    "tool_error_mas": r["target_residual_total_mas"] / sig if sig > 0 else None,
                    "dmag_over_reference_sigma": abs(r["target_dmag"]) / r["reference_dmag_sigma"]
                    if r["target_dmag"] is not None
                    else None,
                }
                null = empirical_null(c1, c2, cat1, k1, k2, r["frame_offset_mas"])
                null["p95_motion_mas_per_yr"] = (
                    null["peer_residual_p95_mas"] / BASELINE_YR
                    if null["peer_residual_p95_mas"] is not None
                    else None
                )
                entry["empirical_null"] = null
            except ValueError as exc:
                r = None
                entry["not_measured"] = str(exc)
            per_band[b] = r
            rec["bands"][b] = entry
        rec["decision"] = decide(per_band)
        out["targets"][name] = rec
    out["assumptions"] = {
        "motion_sig": MOTION_SIG,
        "no_motion_sig": NO_MOTION_SIG,
        "no_motion_bound_mas": NO_MOTION_BOUND_MAS,
        "angle_deg": ANGLE_DEG,
        "var_factor": VAR_FACTOR,
        "null_dmag": NULL_DMAG,
        "null_size_frac": NULL_SIZE_FRAC,
        "null_radius_arcsec": NULL_RADIUS_ARCSEC,
        "baseline_yr": BASELINE_YR,
    }
    args.out.write_text(json.dumps(out, indent=1, default=float) + "\n")
    for name, rec in out["targets"].items():
        print(name, rec["decision"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
