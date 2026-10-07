"""Gather vetting evidence for ranked candidates of one pipeline run (used by /vet-candidate).

For each ``--uid`` this writes ``<out>/<uid>/evidence.json`` and a multi-band contact sheet and
prints a summary. It collects:

* catalog rows per band (observed) from the run's ``sources`` table;
* geometry relative to a cluster centre (derived): radius, position angle, and tangential
  alignment of the major axis from the catalog ``sky_orientation`` and, independently, from
  second moments of the reference-band cutout;
* the nearest image of a published multiple-image catalog (Lenstool ``arcs.dat`` format);
* the nearest star of the run's star/galaxy classification (bright-star neighbours);
* SIMBAD/NED/Gaia cross-matches at two radii (observed);
* cutouts in every band of the sample (S3 byte ranges, never full downloads).

Example (program 2736, Mahler et al. 2022 model, CC0):

    python scripts/vet_evidence.py --run-dir $JWST_ANOMALY_OUTPUTS/runs/<run_id> \\
        --sample smacs0723_nircam --config configs/reference_sample.yaml \\
        --center 110.826750 -73.454628 --lens-images "$MAHLER_ARCS" \\
        --uid jw02736-o001_t001_nircam_f200w_2925

where ``MAHLER_ARCS`` is the ``ICLv0/arcs.dat`` raw URL listed in SOURCES.md ("Vetting").
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import yaml
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from scipy import ndimage

from jwst_anomaly import crossmatch, cutouts, paths, pipeline, schema, viz
from jwst_anomaly.features import column_as_float

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


def axis_offset_deg(pa_a: float, pa_b: float) -> float:
    """Smallest angle between two axes (orientation mod 180), in [0, 90]."""
    d = (pa_a - pa_b) % 180.0
    return float(min(d, 180.0 - d))


def moment_orientation(path: str | Path, nsigma: float = 3.0) -> dict[str, float] | None:
    """Major-axis position angle (deg E of N) and axis ratio from the central segment.

    Pixels above median + nsigma * 1.4826 * MAD are segmented; the segment containing the
    cutout centre is measured with flux-weighted second moments (derived).
    """
    with fits.open(path) as hdul:
        data = np.asarray(hdul[0].data if hdul[0].data is not None else hdul["SCI"].data, float)
        header = hdul[0].header if hdul[0].data is not None else hdul["SCI"].header
    wcs = WCS(header).celestial
    finite = np.isfinite(data)
    if finite.sum() < 10:
        return None
    med = np.median(data[finite])
    mad = 1.4826 * np.median(np.abs(data[finite] - med))
    mask = finite & (data > med + nsigma * mad)
    labels, _ = ndimage.label(mask)
    cy, cx = (np.array(data.shape) - 1) / 2
    lab = labels[int(round(cy)), int(round(cx))]
    if lab == 0:
        return None
    yy, xx = np.nonzero(labels == lab)
    w = data[yy, xx] - med
    xm, ym = np.average(xx, weights=w), np.average(yy, weights=w)
    cov = np.cov(np.vstack([xx - xm, yy - ym]), aweights=w)
    vals, vecs = np.linalg.eigh(cov)
    major = vecs[:, 1]
    p0 = wcs.pixel_to_world(xm, ym)
    p1 = wcs.pixel_to_world(xm + major[0], ym + major[1])
    pa = float(p0.position_angle(p1).to_value(u.deg)) % 180.0
    return {"pa_deg": pa, "axis_ratio": float(np.sqrt(vals[0] / vals[1])), "n_pix": int(len(xx))}


def load_lens_images(source: str) -> Table:
    """Lenstool arcs.dat: id ra dec a b theta z mag (whitespace separated)."""
    if source.startswith(("http://", "https://")):
        with urllib.request.urlopen(source, timeout=60) as resp:  # noqa: S310 (fixed URL)
            text = resp.read().decode()
    else:
        text = Path(source).read_text(encoding="utf-8")
    rows = [line.split()[:3] for line in text.splitlines() if line.strip() and line[0] != "#"]
    return Table(
        {
            "image_id": [r[0] for r in rows],
            "ra": [float(r[1]) for r in rows],
            "dec": [float(r[2]) for r in rows],
        }
    )


def band_columns(sources: Table) -> list[str]:
    return sorted({c.split("_", 1)[0].upper() for c in sources.colnames if c.endswith("_detected")})


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--uid", action="append", required=True)
    ap.add_argument("--center", nargs=2, type=float, metavar=("RA", "DEC"))
    ap.add_argument("--lens-images", help="Lenstool arcs.dat path or URL")
    ap.add_argument("--size-arcsec", type=float, default=4.0)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)

    run_dir = args.run_dir
    sources = Table.read(run_dir / args.sample / "sources.ecsv")
    scores = Table.read(run_dir / args.sample / "scores.ecsv")
    pops_path = run_dir / args.sample / "populations.ecsv"
    pops = Table.read(pops_path) if pops_path.exists() else None
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    sample = next(s for s in config["samples"] if str(s["id"]) == args.sample)
    out_root = args.out or paths.outputs_dir() / "vetting" / run_dir.name
    uids = [str(x) for x in sources["source_uid"]]
    bands = band_columns(sources)
    lens = load_lens_images(args.lens_images) if args.lens_images else None
    center = SkyCoord(args.center[0], args.center[1], unit="deg") if args.center else None
    stars = None
    if pops is not None:
        is_star = np.asarray(pops["population"]).astype(str) == "star"
        idx = [uids.index(str(u_)) for u_ in pops["source_uid"][is_star]]
        stars = SkyCoord(sources["ra"][idx], sources["dec"][idx], unit="deg") if idx else None

    targets = Table(
        {
            "source_uid": args.uid,
            "ra": [float(sources["ra"][uids.index(x)]) for x in args.uid],
            "dec": [float(sources["dec"][uids.index(x)]) for x in args.uid],
        }
    )
    targets["ra"].unit = targets["dec"].unit = "deg"
    targets.meta.update(
        provenance=schema.Provenance.DERIVED.value, source=f"vetting {run_dir.name}"
    )
    xm = {r: crossmatch.query_matches(targets, radius_arcsec=r) for r in (1.0, 3.0)}

    cut_tables = []
    for obs_id in sample["obs_ids"]:
        uri = pipeline.l3_image_uri(config.get("cloud"), obs_id)
        t = cutouts.make_cutouts(
            uri, targets, size_arcsec=args.size_arcsec, out_dir=out_root / "cutouts"
        )
        cut_tables.append(t)
    all_cuts = Table(np.hstack([np.asarray(t.as_array()) for t in cut_tables]))
    ref_band = str(sample["ref_band"]).upper()

    for uid in args.uid:
        i = uids.index(uid)
        row = sources[i]
        pos = SkyCoord(row["ra"], row["dec"], unit="deg")
        ev: dict[str, Any] = {
            "source_uid": uid,
            "run_id": run_dir.name,
            "sample": args.sample,
            "queried_utc": datetime.now(UTC).isoformat(timespec="seconds"),
            "ra": float(row["ra"]),
            "dec": float(row["dec"]),
        }
        k = list(scores["source_uid"]).index(uid)
        ev["score"] = {
            "rank": int(scores["rank"][k]),
            "score": float(scores["score"][k]),
            "top_features": str(scores["top_features"][k]),
            "label": "model_prediction",
        }
        ev["catalog"] = {
            b: {
                q: (
                    float(column_as_float(sources[i : i + 1], f"{b.lower()}_{q}")[0])
                    if f"{b.lower()}_{q}" in sources.colnames
                    else None
                )
                for q in QUANTITIES
            }
            for b in bands
        }
        ref_cut = [
            r for r in all_cuts if str(r["source_uid"]) == uid and str(r["band"]) == ref_band
        ]
        moments = moment_orientation(ref_cut[0]["path"]) if ref_cut and ref_cut[0]["path"] else None
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
                    None
                    if cat_pa is None or not np.isfinite(cat_pa)
                    else axis_offset_deg(cat_pa, tangential)
                ),
                "moments": moments,
                "offset_from_tangential_moments_deg": (
                    None if moments is None else axis_offset_deg(moments["pa_deg"], tangential)
                ),
                "label": "derived",
            }
        if lens is not None:
            lc = SkyCoord(lens["ra"], lens["dec"], unit="deg")
            j = int(np.argmin(pos.separation(lc)))
            ev["lens_images"] = {
                "catalog": args.lens_images,
                "nearest_image": str(lens["image_id"][j]),
                "separation_arcsec": float(pos.separation(lc[j]).arcsec),
                "label": "observed",
            }
        if stars is not None:
            sep = pos.separation(stars).arcsec
            sep = sep[sep > 0.05]
            ev["nearest_star_arcsec"] = float(sep.min()) if sep.size else None
        ev["crossmatch"] = {
            str(r): [
                {
                    "service": str(m["service"]),
                    "id": str(m["match_id"]),
                    "type": str(m["match_type"]),
                    "sep_arcsec": round(float(m["sep_arcsec"]), 3),
                }
                for m in tab
                if str(m["source_uid"]) == uid
            ]
            for r, tab in xm.items()
        }
        ev["crossmatch_meta"] = {
            "services_failed": list(xm[3.0].meta.get("services_failed") or []),
            "query_utc": str(xm[3.0].meta.get("query_utc", "")),
        }
        mine = all_cuts[np.asarray(all_cuts["source_uid"]).astype(str) == uid]
        ev["cutouts"] = {
            str(r["band"]): {"path": str(r["path"]), "quality_flag": str(r["quality_flag"])}
            for r in mine
        }
        out_dir = out_root / uid
        out_dir.mkdir(parents=True, exist_ok=True)
        mine_t = Table(mine)
        mine_t.meta.update(provenance="observed", source=f"cutouts of {uid}")
        ev["contact_sheet"] = str(
            viz.contact_sheet(mine_t, out_dir / "bands.png", title=f"{uid} (unvetted)")
        )
        (out_dir / "evidence.json").write_text(json.dumps(ev, indent=2), encoding="utf-8")
        g = ev.get("geometry", {})
        li = ev.get("lens_images", {})
        print(
            f'{uid}: rank {ev["score"]["rank"]}, r={g.get("radius_arcsec", float("nan")):.1f}", '
            f"tangential offset cat={g.get('offset_from_tangential_catalog_deg')} "
            f"mom={g.get('offset_from_tangential_moments_deg')}, nearest multiple image "
            f'{li.get("nearest_image")} at {li.get("separation_arcsec", float("nan")):.1f}", '
            f'nearest star {ev.get("nearest_star_arcsec")}", xm3={len(ev["crossmatch"]["3.0"])}'
        )
        print(f"  -> {out_dir / 'evidence.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
