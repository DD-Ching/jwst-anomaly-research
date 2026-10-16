"""Arc axes re-measured after removing the local background (ICL) gradient.

An ordinary-explanation test for ``arcs`` anti rows.

A faint source on a steep intracluster-light gradient can get an isophotal segment that extends
toward the brighter side, i.e. radially. This script re-measures the major axis of the
strong-shear background rows of a stored ``arcs`` table on S3 byte-range cutouts of the same
band: a ``--filter-arcsec`` median filter is subtracted (the halo-filter idea of D-029), then
Gaussian-weighted second moments (``--weight-px`` σ, centred on the brightest pixel within
``--peak-px`` of the catalogue position) give the axis. ASSUMPTIONs: the filter size, weight and
peak radius.

Known-case check first: rows whose pipeline axis is ``aligned`` with the z_s = 2 prediction must
stay aligned with the re-measured axis against the same prediction (the method's own calibration;
thresholds from ``arcs.json``). Each row also gets its offset from the radial
direction to the model centre. Output: ``results/arcs/<model>/bgsub_shapes.ecsv`` and
``bgsub_shapes.json`` (derived; pixels are observed, predictions model_prediction).

    python scripts/arcs_bgsub_shapes.py --model abells1063-cats \\
        --image s3://stpubdata/jwst/public/jw03293/L3/t/o001/jw03293-o001_t001_nircam_clear-f200w_i2d.fits
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from scipy.ndimage import median_filter

from jwst_anomaly import cutouts, lensmodel, paths

_spec = importlib.util.spec_from_file_location(
    "lens_consistency", Path(__file__).resolve().parent / "lens_consistency.py"
)
lc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(lc)


axis_offset = lensmodel.axis_offset_deg


def axis_class(off, aligned_deg: float, anti_deg: float) -> np.ndarray:
    """``aligned`` (≤ aligned_deg), ``anti`` (≥ anti_deg) or ``mixed`` (also for NaN offsets)."""
    off = np.asarray(off, float)
    return np.where(off <= aligned_deg, "aligned", np.where(off >= anti_deg, "anti", "mixed"))


def weighted_axis(img: np.ndarray, x0: float, y0: float, sigma_px: float) -> tuple[float, float]:
    """Major-axis angle (deg, counter-clockwise from +x, mod 180) and ellipticity 1 - b/a.

    Gaussian-weighted second moments about (x0, y0); negative pixels are kept (unbiased
    under noise after background subtraction).
    """
    yy, xx = np.indices(img.shape, dtype=float)
    dx, dy = xx - x0, yy - y0
    w = np.exp(-0.5 * (dx**2 + dy**2) / sigma_px**2) * np.nan_to_num(img)
    s = w.sum()
    if s <= 0:
        return float("nan"), float("nan")
    qxx, qyy, qxy = (w * dx * dx).sum() / s, (w * dy * dy).sum() / s, (w * dx * dy).sum() / s
    theta = 0.5 * np.degrees(np.arctan2(2 * qxy, qxx - qyy)) % 180.0
    tr, det = qxx + qyy, qxx * qyy - qxy**2
    disc = np.sqrt(max(tr * tr / 4 - det, 0.0))
    l1, l2 = tr / 2 + disc, tr / 2 - disc
    e = 1.0 - np.sqrt(l2 / l1) if l1 > 0 and l2 > 0 else float("nan")
    return float(theta), float(e)


def stamp_pa(path: str, ra: float, dec: float, args) -> tuple[float, float]:
    """Sky PA (deg E of N, mod 180) and ellipticity on one background-subtracted cutout."""
    with fits.open(path) as h:
        sci = np.asarray(h["SCI"].data, float)
        wcs = WCS(h["SCI"].header)
    px = float(np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))) * 3600.0
    sub = sci - median_filter(np.nan_to_num(sci), size=int(round(args.filter_arcsec / px)) | 1)
    x, y = wcs.world_to_pixel_values(ra, dec)
    r = args.peak_px
    xi, yi = int(round(float(x))), int(round(float(y)))
    x_lo, y_lo = max(xi - r, 0), max(yi - r, 0)
    box = sub[y_lo : yi + r + 1, x_lo : xi + r + 1]
    if box.size == 0 or not np.isfinite(box).any():
        return float("nan"), float("nan")
    j, i = np.unravel_index(np.nanargmax(box), box.shape)
    x0, y0 = x_lo + i, y_lo + j
    theta, e = weighted_axis(sub, x0, y0, args.weight_px)
    # pixel-frame angle -> sky PA through the WCS (one step along the axis)
    t = np.radians(theta)
    c1 = wcs.pixel_to_world(x0, y0)
    c2 = wcs.pixel_to_world(x0 + np.cos(t), y0 + np.sin(t))
    return float(c1.position_angle(c2).deg % 180.0), e


def run(args) -> dict:
    d = args.arcs_dir / args.model
    t = Table.read(d / "arcs.ecsv")
    t = t[np.asarray(t["strong_shear"], bool) & (t["population"] == "background")]
    stored = json.loads((d / "arcs.json").read_text())["assumptions"]
    model, _, _ = lc.load_model(args.model)
    lc.apply_frame_offset(args.model, model)
    ra0, dec0 = model.to_sky(np.array([0.0]), np.array([0.0]))
    centre = SkyCoord(ra0[0], dec0[0], unit="deg")
    tg = Table({"source_uid": [str(v) for v in t["label"]], "ra": t["ra"], "dec": t["dec"]})
    with tempfile.TemporaryDirectory() as tmp:
        cut = cutouts.make_cutouts(args.image, tg, size_arcsec=args.size_arcsec, out_dir=Path(tmp))
        res = [
            stamp_pa(p, r["ra"], r["dec"], args) if p else (float("nan"), float("nan"))
            for p, r in zip(cut["path"], t, strict=True)
        ]
    t["pa_bgsub"] = [v[0] for v in res]
    t["e_bgsub"] = [v[1] for v in res]
    pa_rad = centre.position_angle(SkyCoord(t["ra"], t["dec"], unit="deg")).deg % 180.0
    t["dpa_radial_obs"] = axis_offset(t["pa_obs"], pa_rad)
    t["dpa_radial_bgsub"] = axis_offset(t["pa_bgsub"], pa_rad)
    t["dpa_pipeline_bgsub"] = axis_offset(t["pa_obs"], t["pa_bgsub"])
    # both axes against the same z_s = 2 prediction, so the control isolates the measurement
    lim = (stored["aligned_deg"], stored["anti_deg"])
    t["class_obs_z2"] = axis_class(axis_offset(t["pa_obs"], t["pa_pred_z2"]), *lim)
    t["offset_pred_bgsub"] = axis_offset(t["pa_bgsub"], t["pa_pred_z2"])
    t["class_bgsub"] = axis_class(t["offset_pred_bgsub"], *lim)
    keep = [
        "label",
        "ra",
        "dec",
        "mag",
        "orientation_class",
        "class_obs_z2",
        "pa_obs",
        "pa_bgsub",
        "e_bgsub",
        "pa_pred_z2",
        "offset_pred_bgsub",
        "class_bgsub",
        "dpa_pipeline_bgsub",
        "dpa_radial_obs",
        "dpa_radial_bgsub",
    ]
    out = t[keep]
    out.meta.update(
        provenance="derived",
        source=f"{args.model} arcs.ecsv strong-shear background + {args.image}",
        assumptions={k: getattr(args, k) for k in ("filter_arcsec", "weight_px", "peak_px")},
    )
    out.write(d / "bgsub_shapes.ecsv", overwrite=True)
    ctrl = out[out["class_obs_z2"] == "aligned"]
    anti = out[out["orientation_class"] == "anti"]

    def num(v):
        v = float(v)
        return round(v, 2) if np.isfinite(v) else None

    summary = {
        "model": args.model,
        "image": args.image,
        "label": "ordinary-explanation test (derived); predictions model_prediction",
        "assumptions": out.meta["assumptions"]
        | {"aligned_deg": lim[0], "anti_deg": lim[1], "pred_z": 2.0},
        "centre_radec": [round(float(ra0[0]), 6), round(float(dec0[0]), 6)],
        "control_aligned_z2": {
            "n": len(ctrl),
            "stay_aligned": int(np.sum(ctrl["class_bgsub"] == "aligned")),
            "median_dpa_pipeline_bgsub": num(np.median(ctrl["dpa_pipeline_bgsub"]))
            if len(ctrl)
            else None,
        },
        "anti": [
            {
                "label": int(r["label"]),
                "class_obs_z2": str(r["class_obs_z2"]),
                "class_bgsub": str(r["class_bgsub"]),
            }
            | {
                k: num(r[k])
                for k in (
                    "pa_obs",
                    "pa_bgsub",
                    "e_bgsub",
                    "offset_pred_bgsub",
                    "dpa_radial_obs",
                    "dpa_radial_bgsub",
                )
            }
            for r in anti
        ],
    }
    (d / "bgsub_shapes.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--model", required=True)
    p.add_argument("--image", required=True, help="_i2d URI (S3 byte ranges) of the shape band")
    p.add_argument("--arcs-dir", type=Path, default=paths.repo_root() / "results" / "arcs")
    p.add_argument("--size-arcsec", type=float, default=4.0)
    p.add_argument("--filter-arcsec", type=float, default=1.0, help="ASSUMPTION")
    p.add_argument("--weight-px", type=float, default=3.0, help="ASSUMPTION")
    p.add_argument("--peak-px", type=int, default=3, help="ASSUMPTION")
    args = p.parse_args(argv)
    print(json.dumps(run(args), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
