"""Two-epoch astrometry and photometry of one target from level-3 pipeline catalogs (vetting).

Compares a target's position and magnitude in two catalogs of the same filter taken at different
epochs (e.g. program 2736 in 2022 and a later program covering the same sky). Extragalactic
sources do not move, so the median offset of matched sources around the target fixes the frame
tie between the two catalogs; its error is about 1.25 sigma / sqrt(N) per axis. The target's
residual from that median, compared with the frame-tie error (a bright point source centroids far
better than the galaxy references), tests for proper motion (stars, brown dwarfs). The same is
done for the magnitude difference (variability; zero points cancel in the median). Every result
is ``derived``. Optionally lists the target's magnitude in every catalog of a second observation
(e.g. medium bands).

    python scripts/epoch_compare.py --target 110.668714 -73.504767 \\
        --epoch1 $JWST_ANOMALY_DATA/cache/mast/<obs_id>/<obs_id>_cat.ecsv \\
        --epoch2-obs-id jw06882-o057_t057_nircam_clear-f444w \\
        --sed-obs-prefix jw06882-o057_t057_nircam_clear-

``--epoch2-obs-id`` and ``--sed-obs-prefix`` fetch level-3 ``CAT`` products from MAST into the
data cache, recording them in ``--manifest`` (default under the outputs directory, not in
``data/manifests/``).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time

from jwst_anomaly import acquire, catalog, paths, query
from jwst_anomaly.features import snr_from_mag_err

MAD_TO_SIGMA = 1.4826


def mutual_matches(a: SkyCoord, b: SkyCoord, radius_arcsec: float) -> tuple[np.ndarray, np.ndarray]:
    """Index pairs ``(i_a, i_b)`` that are each other's nearest neighbour within the radius."""
    j, sep, _ = a.match_to_catalog_sky(b)
    i_back, _, _ = b.match_to_catalog_sky(a)
    ia = np.arange(len(a))
    keep = (sep.arcsec <= radius_arcsec) & (i_back[j] == ia)
    return ia[keep], j[keep]


def offsets_mas(a: SkyCoord, b: SkyCoord) -> np.ndarray:
    """``(n, 2)`` offsets b - a: (dRA cos Dec, dDec) in mas."""
    dra, ddec = a.spherical_offsets_to(b)
    return np.column_stack([dra.to_value(u.mas), ddec.to_value(u.mas)])


def robust_center(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Median and MAD-based sigma along axis 0."""
    med = np.median(values, axis=0)
    return med, MAD_TO_SIGMA * np.median(np.abs(values - med), axis=0)


def compare_epochs(
    cat1: Table,
    cat2: Table,
    target: SkyCoord,
    *,
    match_radius_arcsec: float = 0.3,
    target_radius_arcsec: float = 0.3,
    ref_radius_arcsec: float = 60.0,
    min_snr: float = 30.0,
    n_ref_min: int = 20,
) -> dict[str, Any]:
    """The target's astrometric and photometric residuals relative to nearby matched sources.

    References are mutual matches within ``ref_radius_arcsec`` of the target with aper50 S/N >=
    ``min_snr`` in both catalogs, excluding the target itself. The target's two rows must be each
    other's nearest neighbours, so a neighbour cannot stand in for it in one epoch. Its centroid
    error per epoch is ``semimajor_sigma * pixel scale / S/N`` (when the catalog has them) and is
    added in quadrature to the frame-tie error for ``target_residual_significance``. Raises
    ``ValueError`` when the target is missing or ambiguous, or fewer than ``n_ref_min``
    references remain.
    """
    c1 = SkyCoord(cat1["ra"], cat1["dec"], unit="deg")
    c2 = SkyCoord(cat2["ra"], cat2["dec"], unit="deg")
    t1 = int(np.argmin(target.separation(c1).arcsec))
    t2 = int(np.argmin(target.separation(c2).arcsec))
    for name, c, k in (("epoch1", c1, t1), ("epoch2", c2, t2)):
        if target.separation(c[k]).arcsec > target_radius_arcsec:
            raise ValueError(f"target not within {target_radius_arcsec} arcsec in {name}")
    back1 = int(np.argmin(c2[t2].separation(c1).arcsec))
    back2 = int(np.argmin(c1[t1].separation(c2).arcsec))
    if back1 != t1 or back2 != t2:
        raise ValueError("ambiguous target: its epoch rows are not each other's nearest neighbours")
    i1, i2 = mutual_matches(c1, c2, match_radius_arcsec)
    snr1 = snr_from_mag_err(np.asarray(cat1["aper50_abmag_err"], float))[i1]
    snr2 = snr_from_mag_err(np.asarray(cat2["aper50_abmag_err"], float))[i2]
    near = target.separation(c1[i1]).arcsec <= ref_radius_arcsec
    mag1 = np.asarray(cat1["aper50_abmag"], float)
    mag2 = np.asarray(cat2["aper50_abmag"], float)
    finite = np.isfinite(mag1[i1]) & np.isfinite(mag2[i2])
    # The target is excluded by both of its rows: a moved target may pair with another source.
    ref = near & finite & (snr1 >= min_snr) & (snr2 >= min_snr) & (i1 != t1) & (i2 != t2)
    if ref.sum() < n_ref_min:
        raise ValueError(f"only {int(ref.sum())} reference sources (need {n_ref_min})")
    d_ref = offsets_mas(c1[i1[ref]], c2[i2[ref]])
    med, sig = robust_center(d_ref)
    d_t = offsets_mas(c1[[t1]], c2[[t2]])[0] - med
    dm_ref = mag2[i2[ref]] - mag1[i1[ref]]
    dm_med, dm_sig = robust_center(dm_ref)
    dm_t = float(mag2[t2] - mag1[t1]) - float(dm_med)
    total = float(np.hypot(*d_t))
    tie = 1.2533 * sig / np.sqrt(ref.sum())  # standard error of a median
    snr_t = [float(snr_from_mag_err(c["aper50_abmag_err"][k])) for c, k in ((cat1, t1), (cat2, t2))]
    cen = [_centroid_error_mas(c, k, s) for c, k, s in ((cat1, t1, snr_t[0]), (cat2, t2, snr_t[1]))]
    cen_total = float(np.hypot(*cen)) if all(x is not None for x in cen) else None
    err = np.hypot(tie, cen_total / np.sqrt(2)) if cen_total is not None else tie  # per axis
    resid_ref = np.hypot(*(d_ref - med).T)
    return {
        "n_references": int(ref.sum()),
        "frame_offset_mas": [float(x) for x in med],
        "reference_sigma_mas": [float(x) for x in sig],
        "target_residual_mas": [float(x) for x in d_t],
        "target_residual_total_mas": total,
        "frame_tie_error_mas": [float(x) for x in tie],
        "target_centroid_error_mas": cen,
        "target_residual_significance": float(np.hypot(*(d_t / err))),
        "target_residual_over_reference_sigma": float(np.hypot(*(d_t / sig))),
        "fraction_of_references_with_larger_residual": float(np.mean(resid_ref > total)),
        "target_dmag": dm_t,
        "reference_dmag_sigma": float(dm_sig),
        "target_snr": snr_t,
        "label": "derived",
    }


def _centroid_error_mas(cat: Table, k: int, snr: float) -> float | None:
    """Approximate centroid error, ``semimajor_sigma * pixel scale / S/N``, in mas (or None)."""
    scale = cat.meta.get("pixel_scale_arcsec")
    if "semimajor_sigma" not in cat.colnames or not scale or not np.isfinite(snr) or snr <= 0:
        return None
    return float(np.asarray(cat["semimajor_sigma"], float)[k] * float(scale) * 1000.0 / snr)


def _fetch_cats(obs_filter: str, exact: bool, manifest: Path) -> list[Path]:
    obs = query.query_observations(obs_id=obs_filter if exact else f"{obs_filter}*")
    prods = query.list_products(obs, subgroups=("CAT",))
    if len(prods) == 0:
        raise SystemExit(f"error: no public level-3 CAT products for {obs_filter!r}")
    rows = acquire.fetch_products(prods, manifest_path=manifest)
    return [paths.data_root() / str(r["local_path"]) for r in rows]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--target", nargs=2, type=float, required=True, metavar=("RA", "DEC"))
    ap.add_argument("--epoch1", type=Path, required=True, help="level-3 _cat.ecsv, first epoch")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--epoch2", type=Path, help="level-3 _cat.ecsv, second epoch (same filter)")
    g.add_argument("--epoch2-obs-id", help="fetch the second epoch's CAT product from MAST")
    ap.add_argument("--sed-obs-prefix", help="also list the target in every CAT of these obs_ids")
    ap.add_argument("--manifest", type=Path, default=None)
    ap.add_argument("--ref-radius", type=float, default=60.0)
    ap.add_argument("--min-snr", type=float, default=30.0)
    ap.add_argument("--out", type=Path, default=None, help="write the result as JSON")
    args = ap.parse_args(argv)

    manifest = args.manifest or paths.outputs_dir() / "vet" / "epoch_manifest.ecsv"
    target = SkyCoord(args.target[0], args.target[1], unit="deg")
    cat1 = catalog.load_pipeline_catalog(args.epoch1)
    path2 = args.epoch2 or _fetch_cats(args.epoch2_obs_id, True, manifest)[0]
    cat2 = catalog.load_pipeline_catalog(path2)
    if cat1.meta["band"] != cat2.meta["band"]:
        ap.error(f"filters differ: {cat1.meta['band']} vs {cat2.meta['band']}")
    result: dict[str, Any] = {
        "target": list(args.target),
        "epoch1": {"file": Path(args.epoch1).name, "sha256": cat1.meta.get("input_sha256")},
        "epoch2": {"file": Path(path2).name, "sha256": cat2.meta.get("input_sha256")},
        "band": cat1.meta["band"],
        **compare_epochs(
            cat1, cat2, target, ref_radius_arcsec=args.ref_radius, min_snr=args.min_snr
        ),
    }
    if args.sed_obs_prefix:
        sed = {}
        for p in _fetch_cats(args.sed_obs_prefix, False, manifest):
            c = catalog.load_pipeline_catalog(p)
            sky = SkyCoord(c["ra"], c["dec"], unit="deg")
            k = int(np.argmin(target.separation(sky).arcsec))
            sep = float(target.separation(sky[k]).arcsec)
            sed[c.meta["band"]] = {
                "aper50_abmag": float(c["aper50_abmag"][k]) if sep <= 0.3 else None,
                "aper50_abmag_err": float(c["aper50_abmag_err"][k]) if sep <= 0.3 else None,
                "sep_arcsec": sep,
                "file": p.name,
                "sha256": c.meta.get("input_sha256"),
            }
        result["sed"] = sed
    result["script_utc"] = Time.now().isot
    text = json.dumps(result, indent=2)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
