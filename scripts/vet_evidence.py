"""Gather vetting evidence for ranked candidates of one pipeline run (used by /vet-candidate).

For each ``--uid`` this writes ``<out>/<uid>/evidence.json`` (strict JSON; NaN -> null) and a
multi-band contact sheet, and prints a summary. Stages that fail are recorded under ``errors``
and the rest still run. Evidence collected:

* the run's own score (model_prediction) and D-011 quality-gate row;
* catalog rows per band (observed) from the run's ``sources`` table;
* geometry relative to a cluster centre (derived): radius, position angle, and tangential
  alignment of the major axis from the catalog ``sky_orientation`` and, independently, from
  second moments of the reference-band cutout (grown until the segment is closed);
* the nearest image of a published multiple-image catalog (Lenstool ``arcs.dat``, absolute
  coordinates). Absence from such a catalog only means "not a constraint image of that model";
* the nearest star of the run's star/galaxy classification;
* SIMBAD/NED/Gaia matches (observed): one query at the largest radius, filtered for smaller ones,
  with redshifts, lens flags, query times and cache status;
* cutouts in every band of the sample (S3 byte ranges; the reference band uses the gate's
  reference weight).

Example (program 2736 with the Mahler et al. 2022 model, pinned commit; see SOURCES.md):

    python scripts/vet_evidence.py --run-dir $JWST_ANOMALY_OUTPUTS/runs/<run_id> \\
        --sample smacs0723_nircam --config configs/reference_sample.yaml \\
        --center 110.826750 -73.454628 --lens-images "$MAHLER_ARCS" \\
        --uid jw02736-o001_t001_nircam_f200w_2925

Use ``--sample <sample>-stars`` for candidates of a star stratum (D-012).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.wcs import WCS
from scipy import ndimage

from jwst_anomaly import crossmatch, cutouts, paths, pipeline, schema, viz
from jwst_anomaly.features import column_as_float, discover_bands

QUANTITIES = (
    "detected",
    "sep_arcsec",
    "aper50_abmag",
    "aper50_abmag_err",
    "isophotal_abmag",
    "CI_50_30",
    "ellipticity",
    "sky_orientation",
    "isophotal_area",
    "sharpness",
    "roundness",
    "nn_dist",
    "is_extended",
)
MAX_GROWTH = 3  # cutout doublings while the moments segment touches the border


def axis_offset_deg(pa_a: float, pa_b: float) -> float:
    """Smallest angle between two axes (orientation mod 180), in [0, 90]."""
    d = (pa_a - pa_b) % 180.0
    return float(min(d, 180.0 - d))


def moment_orientation(path: str | Path, nsigma: float = 3.0) -> dict[str, Any] | None:
    """Major-axis position angle (deg E of N) and axis ratio of the central segment.

    Pixels above median + nsigma * 1.4826 * MAD are segmented (no deblending); the segment
    containing the cutout centre is measured with flux-weighted second moments (derived).
    ``touches_border`` marks a segment clipped by the cutout, whose moments are biased.
    """
    with fits.open(path) as hdul:
        hdu = hdul[0] if hdul[0].data is not None else hdul["SCI"]
        data = np.asarray(hdu.data, float)
        wcs = WCS(hdu.header).celestial
    finite = np.isfinite(data)
    if finite.sum() < 10:
        return None
    med = np.median(data[finite])
    mad = 1.4826 * np.median(np.abs(data[finite] - med))
    labels, _ = ndimage.label(finite & (data > med + nsigma * mad))
    cy, cx = (np.array(data.shape) - 1) / 2
    lab = labels[int(round(cy)), int(round(cx))]
    if lab == 0:
        return None
    seg = labels == lab
    yy, xx = np.nonzero(seg)
    if len(xx) < 5:
        return None
    border = np.zeros_like(seg)
    border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
    w = data[yy, xx] - med
    xm, ym = np.average(xx, weights=w), np.average(yy, weights=w)
    cov = np.cov(np.vstack([xx - xm, yy - ym]), aweights=w)
    vals, vecs = np.linalg.eigh(cov)
    if not (vals[1] > 0 and np.all(np.isfinite(vals))):
        return None
    major = vecs[:, 1]
    p0 = wcs.pixel_to_world(xm, ym)
    p1 = wcs.pixel_to_world(xm + major[0], ym + major[1])
    return {
        "pa_deg": float(p0.position_angle(p1).to_value(u.deg)) % 180.0,
        "axis_ratio": float(np.sqrt(max(vals[0], 0.0) / vals[1])),
        "n_pix": int(len(xx)),
        "touches_border": bool((seg & border).any()),
        "n_border": int((seg & border).sum()),
    }


def load_lens_images(source: str) -> Table:
    """Lenstool ``arcs.dat`` with absolute coordinates: ``id ra dec [a b theta z mag]``.

    Lines that are blank or start with ``#`` (after whitespace) are skipped; ``#REFERENCE``
    with a non-zero mode (relative arcsec offsets) is rejected rather than misread.
    """
    if source.startswith(("http://", "https://")):
        with urllib.request.urlopen(source, timeout=60) as resp:  # noqa: S310 (user-given URL)
            raw = resp.read()
    else:
        raw = Path(source).read_bytes()
    text = raw.decode()
    ids, ras, decs = [], [], []
    for n, line in enumerate(text.splitlines(), 1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            parts = stripped.lstrip("#").split()
            if parts and parts[0].upper() == "REFERENCE" and len(parts) > 1 and parts[1] != "0":
                raise ValueError(f"{source}: relative Lenstool coordinates (REFERENCE {parts[1]})")
            continue
        parts = stripped.split()
        if len(parts) < 3:
            raise ValueError(f"{source}:{n}: expected 'id ra dec ...', got {stripped!r}")
        ids.append(parts[0])
        ras.append(float(parts[1]))
        decs.append(float(parts[2]))
    out = Table({"image_id": ids, "ra": ras, "dec": decs})
    out.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source=source,
        sha256=hashlib.sha256(raw).hexdigest(),
    )
    return out


def json_safe(obj: Any) -> Any:
    """Recursively replace NaN/inf by None and numpy scalars by Python ones."""
    if isinstance(obj, dict):
        return {str(k): json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list | tuple):
        return [json_safe(v) for v in obj]
    if isinstance(obj, np.generic):
        obj = obj.item()
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj


def _read_run_table(sample_dir: Path, name: str, required: bool = True) -> Table | None:
    for ext in ("ecsv", "parquet"):
        path = sample_dir / f"{name}.{ext}"
        if path.exists():
            return Table.read(path)
    if required:
        raise SystemExit(f"error: {sample_dir}/{name}.(ecsv|parquet) not found")
    return None


def _moments_with_growth(
    uri: str, target: Table, size: float, out_root: Path
) -> tuple[dict[str, Any] | None, str | None]:
    """Moments of the reference-band cutout, doubling the box while the segment is clipped."""
    moments = None
    for attempt in range(MAX_GROWTH + 1):
        used = size * 2**attempt
        try:
            t = cutouts.make_cutouts(
                uri, target, size_arcsec=used, out_dir=out_root / "moments" / f"{used:g}arcsec"
            )
            moments = moment_orientation(t["path"][0]) if t["path"][0] else None
        except Exception as exc:
            return None, f"{type(exc).__name__}: {exc}"
        if moments is not None:
            moments["cutout_arcsec"] = used  # the size actually measured
        if moments is None or not moments["touches_border"]:
            break
    return moments, None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--sample", required=True, help="sample id, or '<sample>-stars'")
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--uid", action="append", required=True)
    ap.add_argument("--center", nargs=2, type=float, metavar=("RA", "DEC"))
    ap.add_argument("--lens-images", help="Lenstool arcs.dat path or URL (absolute coordinates)")
    ap.add_argument("--radii", nargs="+", type=float, default=[1.0, 3.0])
    ap.add_argument("--size-arcsec", type=float, default=4.0)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)

    # --- validate everything before any network work --------------------------------------
    config = pipeline.load_config(args.config)
    base = args.sample.removesuffix("-stars")
    sample = next((s for s in config["samples"] if str(s["id"]) == base), None)
    if sample is None:
        ap.error(f"sample {base!r} not in {args.config}")
    if not sample.get("obs_ids"):
        ap.error(f"sample {base!r} has no obs_ids; cutouts need them")
    sources = _read_run_table(args.run_dir / base, "sources")
    scores = _read_run_table(args.run_dir / args.sample, "scores")
    gate = _read_run_table(args.run_dir / base, "quality", required=False)
    pops = _read_run_table(args.run_dir / base, "populations", required=False)
    phot = _read_run_table(args.run_dir / base, "photometry", required=False)
    wanted = list(dict.fromkeys(args.uid))
    uids = [str(x) for x in sources["source_uid"]]
    scored = {str(x) for x in scores["source_uid"]}
    unknown = [x for x in wanted if x not in uids]
    unscored = [x for x in wanted if x in uids and x not in scored]
    if unknown or unscored:
        ap.error(
            f"not in sources: {unknown}; not ranked in {args.sample!r}: {unscored} "
            "(gated out by D-011, or in another stratum: try --sample <sample>-stars)"
        )
    radii = sorted(set(args.radii))
    ref_band = str(sample["ref_band"]).upper()
    bands = discover_bands(sources)  # lower-case, wavelength order
    out_root = args.out or paths.outputs_dir() / "vetting" / args.run_dir.name
    errors: dict[str, str] = {}
    script_utc = datetime.now(UTC).isoformat(timespec="seconds")

    lens = None
    if args.lens_images:
        try:
            lens = load_lens_images(args.lens_images)
        except Exception as exc:  # recorded, not fatal
            errors["lens_images"] = f"{type(exc).__name__}: {exc}"
    center = SkyCoord(args.center[0], args.center[1], unit="deg") if args.center else None
    stars = None
    if pops is not None:
        star_uids = {
            str(p)
            for p, k in zip(pops["source_uid"], pops["population"], strict=True)
            if str(k) == "star"
        }
        idx = [i for i, x in enumerate(uids) if x in star_uids]
        stars = SkyCoord(sources["ra"][idx], sources["dec"][idx], unit="deg") if idx else None

    targets = Table(
        {
            "source_uid": wanted,
            "ra": [float(sources["ra"][uids.index(x)]) for x in wanted],
            "dec": [float(sources["dec"][uids.index(x)]) for x in wanted],
        }
    )
    targets["ra"].unit = targets["dec"].unit = "deg"
    targets.meta.update(
        provenance=schema.Provenance.DERIVED.value, source=f"vetting {args.run_dir.name}"
    )

    matches = None
    try:
        matches = crossmatch.query_matches(targets, radius_arcsec=radii[-1])
    except Exception as exc:
        errors["crossmatch"] = f"{type(exc).__name__}: {exc}"

    qcfg = (config.get("stages") or {}).get("quality") or {}
    uri_of: dict[str, str] = {}
    for obs in sample["obs_ids"]:
        band = pipeline.band_from_name(obs)
        if band:
            uri_of[band.upper()] = pipeline.l3_image_uri(config.get("cloud"), obs)
    weight_ref = None
    if ref_band in uri_of:
        try:
            weight_ref = cutouts.sample_weight_map(
                uri_of[ref_band], grid_arcsec=float(qcfg.get("grid_arcsec", 1.0))
            ).reference_weight
        except Exception as exc:
            errors["weight_map"] = f"{type(exc).__name__}: {exc}"
    cut_tables = []
    for band, uri in uri_of.items():
        try:
            kw = {"weight_ref": weight_ref} if band == ref_band and weight_ref else {}
            cut_tables.append(
                cutouts.make_cutouts(
                    uri, targets, size_arcsec=args.size_arcsec, out_dir=out_root / "cutouts", **kw
                )
            )
        except Exception as exc:
            errors[f"cutouts[{band}]"] = f"{type(exc).__name__}: {exc}"
    all_cuts = vstack(cut_tables, metadata_conflicts="silent") if cut_tables else None

    for uid in wanted:
        i = uids.index(uid)
        row = sources[i]
        pos = SkyCoord(row["ra"], row["dec"], unit="deg")
        ev: dict[str, Any] = {
            "source_uid": uid,
            "run_id": args.run_dir.name,
            "sample": args.sample,
            "script_utc": script_utc,
            "ra": float(row["ra"]),
            "dec": float(row["dec"]),
            "errors": dict(errors),
        }
        k = [str(x) for x in scores["source_uid"]].index(uid)
        ev["score"] = {
            "rank": int(scores["rank"][k]),
            "score": float(scores["score"][k]),
            "top_features": str(scores["top_features"][k]),
            "label": "model_prediction",
        }
        if gate is not None:
            g = gate[[str(x) for x in gate["source_uid"]].index(uid)]
            ev["run_quality_gate"] = {
                "quality_ok": bool(g["quality_ok"]),
                "quality_reason": pipeline._text(g["quality_reason"]),  # "" (masked in ECSV) = ok
                "rel_weight": float(g["rel_weight"]),
                "edge_dist_arcsec": float(g["edge_dist_arcsec"]),
                "label": "derived",
            }
        one = sources[i : i + 1]
        ev["catalog"] = {
            b.upper(): {
                q: (
                    float(column_as_float(one, schema.band_column(b, q))[0])
                    if schema.band_column(b, q) in sources.colnames
                    else None
                )
                for q in QUANTITIES
            }
            for b in bands
        }
        if phot is not None:
            k_phot = [str(x) for x in phot["source_uid"]].index(uid)
            sep_cols = [c for c in phot.colnames if c.endswith("_match_sep_arcsec")]
            label = sep_cols[0].removesuffix("_match_sep_arcsec") if sep_cols else ""
            ev["matched_photometry"] = {
                "label": label,
                "match_sep_arcsec": float(phot[sep_cols[0]][k_phot]) if sep_cols else None,
                "abmag": {
                    b.upper(): float(column_as_float(phot[k_phot : k_phot + 1], c)[0])
                    for b in bands
                    if (c := schema.band_column(b, f"{label}_abmag")) in phot.colnames
                },
                "note": "these colours ranked the source (D-013); aper50 colours are size-biased",
                "label_provenance": "derived",
            }
        ev["cutouts"] = {}
        if all_cuts is not None:
            for r in all_cuts[np.asarray(all_cuts["source_uid"]).astype(str) == uid]:
                band = str(r["band"])
                ev["cutouts"][band] = {
                    "path": str(r["path"]),
                    "quality_flag": str(r["quality_flag"]),
                    "weight_ref": (
                        "D-011 gate map" if band == ref_band and weight_ref else "32 sampled rows"
                    ),
                }
        moments = None
        if ref_band in uri_of:
            moments, err = _moments_with_growth(
                uri_of[ref_band], targets[[wanted.index(uid)]], args.size_arcsec, out_root
            )
            if err:
                ev["errors"]["moments"] = err
        if center is not None:
            pa_radial = float(center.position_angle(pos).to_value(u.deg))
            tangential = (pa_radial + 90.0) % 180.0
            cat_pa = ev["catalog"].get(ref_band, {}).get("sky_orientation")
            ev["geometry"] = {
                "center": [float(center.ra.deg), float(center.dec.deg)],
                "radius_arcsec": float(center.separation(pos).arcsec),
                "pa_from_center_deg": pa_radial,
                "tangential_pa_deg": tangential,
                "catalog_sky_orientation_deg": cat_pa,
                "offset_from_tangential_catalog_deg": (
                    axis_offset_deg(cat_pa, tangential)
                    if cat_pa is not None and math.isfinite(cat_pa)
                    else None
                ),
                "moments": moments,
                "offset_from_tangential_moments_deg": (
                    axis_offset_deg(moments["pa_deg"], tangential) if moments else None
                ),
                "label": "derived",
            }
        if lens is not None:
            lc = SkyCoord(lens["ra"], lens["dec"], unit="deg")
            j = int(np.argmin(pos.separation(lc)))
            ev["lens_images"] = {
                "catalog": args.lens_images,
                "sha256": lens.meta["sha256"],
                "n_images": len(lens),
                "nearest_image": str(lens["image_id"][j]),
                "separation_arcsec": float(pos.separation(lc[j]).arcsec),
                "meaning": "absent = not a constraint image of this model (not: singly imaged)",
                "label": "observed",
            }
        if stars is not None:
            sep = pos.separation(stars).arcsec
            sep = sep[sep > 0.05]
            ev["nearest_star_arcsec"] = float(sep.min()) if sep.size else None
        if matches is not None:
            mine = matches[np.asarray(matches["source_uid"]).astype(str) == uid]
            rows = [
                {
                    "service": str(m["service"]),
                    "id": str(m["match_id"]),
                    "type": str(m["match_type"]),
                    "sep_arcsec": float(m["sep_arcsec"]),
                    "redshift": float(np.ma.filled(m["redshift"], np.nan)),
                    "is_star": bool(m["is_star"]),
                    "is_lens_related": bool(m["is_lens_related"]),
                    "lens_types": pipeline._text(m["lens_types"]),
                }
                for m in mine
            ]
            ev["crossmatch"] = {f"{r:g}": [m for m in rows if m["sep_arcsec"] <= r] for r in radii}
            ev["crossmatch_meta"] = {
                key: matches.meta.get(key)
                for key in (
                    "services_requested",
                    "services_failed",
                    "service_errors",
                    "catalog",
                    "query_utc",
                    "from_cache",
                    "radius_arcsec",
                )
            }
        out_dir = out_root / uid
        out_dir.mkdir(parents=True, exist_ok=True)
        try:
            if all_cuts is not None:
                mine_t = all_cuts[np.asarray(all_cuts["source_uid"]).astype(str) == uid]
                ev["contact_sheet"] = str(
                    viz.contact_sheet(mine_t, out_dir / "bands.png", title=f"{uid} (unvetted)")
                )
        except Exception as exc:
            ev["errors"]["contact_sheet"] = f"{type(exc).__name__}: {exc}"
        (out_dir / "evidence.json").write_text(
            json.dumps(json_safe(ev), indent=2, allow_nan=False), encoding="utf-8"
        )
        geo = ev.get("geometry", {})
        li = ev.get("lens_images", {})
        mo = geo.get("moments") or {}
        print(
            f'{uid}: rank {ev["score"]["rank"]}, r={geo.get("radius_arcsec", float("nan")):.1f}", '
            f"tangential offset cat={geo.get('offset_from_tangential_catalog_deg')} "
            f"mom={geo.get('offset_from_tangential_moments_deg')} "
            f'(cutout {mo.get("cutout_arcsec")}", border={mo.get("touches_border")}), '
            f"nearest constraint image {li.get('nearest_image')} at "
            f'{li.get("separation_arcsec", float("nan")):.1f}", '
            f'nearest star {ev.get("nearest_star_arcsec")}", errors={sorted(ev["errors"])}'
        )
        print(f"  -> {out_dir / 'evidence.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
