"""`arcs` on RXJ2129 with the DJA grizli v7.0 photometric catalogue (a re-test of #283's field).

Plan (frozen before the run): ``results/arcs/rxj2129-dja/PLAN.md``. The DJA SEP catalogue (``ir``
detection image, 0.04" px, north-up) is converted to the JWST pipeline ``_cat.ecsv`` columns that
``lens_consistency.load_shapes`` reads; the PA convention is checked against the program 2767 F200W
pipeline catalogue (the frozen gate ``convention_passes``) before ``arcs_rxj2129.main`` runs the
same maps, rules and outcome as #283. Shapes are ``observed``; the conversion is ``derived``.

    python scripts/arcs_rxj2129_dja.py --dja <rxj2129-grizli-v7.0-fix_phot_apcorr.fits> \
        --pipeline <jw02767-o002_t001_nircam_clear-f200w_cat.ecsv> \
        --photoz <rxj2129-grizli-v7.0-fix.eazypy.zout.fits>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

from jwst_anomaly import lensmodel
from jwst_anomaly.photometry import sha256_file

sys.path.insert(0, str(Path(__file__).resolve().parent))
import arcs_rxj2129  # noqa: E402

DJA_PX_ARCSEC = 0.04  # DJA v7 mosaic pixel (catalogue meta BACK_PIXEL_SCALE)
MATCH_ARCSEC = 0.2  # PLAN.md: validation match radius
MAX_MEDIAN_DPA = 10.0  # PLAN.md: convention gate, degrees (ASSUMPTION)
OUT = Path(__file__).resolve().parents[1] / "results" / "arcs" / "rxj2129-dja"


def dja_pa(theta_image_rad) -> np.ndarray:
    """Sky PA (deg east of north, mod 180) of the major axis for a north-up, east-left image,
    from SEP ``theta_image`` (radians counter-clockwise from +x)."""
    return np.mod(np.degrees(np.asarray(theta_image_rad, float)) + 90.0, 180.0)


def pipeline_px_arcsec(cat: Table) -> float:
    """Pixel scale of the pipeline mosaic from its centroids (median over consecutive rows)."""
    sky = cat["sky_centroid"]
    dx = np.hypot(np.diff(np.asarray(cat["xcentroid"], float)), np.diff(cat["ycentroid"]))
    ds = sky[:-1].separation(sky[1:]).arcsec
    ok = dx > 50
    return float(np.median(ds[ok] / dx[ok]))


def to_pipeline(dja: Table, pipe_px_arcsec: float) -> Table:
    """DJA rows as the pipeline columns ``load_shapes`` reads (PLAN.md 'Conversion')."""
    flux = np.asarray(dja["flux_iso"], float)
    # flux_iso is in detection-image units, not the ZP 28.9 unit, so the AB magnitude is the
    # catalogue's mag_auto (deviation from PLAN.md; mag is not used by the selection or classes)
    mag = np.asarray(dja["mag_auto"], float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ell = 1.0 - np.asarray(dja["b_image"], float) / np.asarray(dja["a_image"], float)
    return Table(
        {
            "label": np.asarray(dja["id"], int),
            "sky_centroid": SkyCoord(dja["ra"], dja["dec"], unit="deg"),
            "isophotal_flux": flux,
            "isophotal_flux_err": np.asarray(dja["fluxerr_iso"], float),
            "sky_orientation": dja_pa(dja["theta_image"]),
            "ellipticity": ell,
            "semimajor_sigma": np.asarray(dja["a_image"], float) * DJA_PX_ARCSEC / pipe_px_arcsec,
            "isophotal_area": np.asarray(dja["area_iso"], float),
            "is_extended": np.ones(len(dja), bool),
            "isophotal_abmag": mag,
        }
    )


def convention_passes(median_dpa: float, n: int) -> bool:
    """PLAN.md gate: at least one matched pair and median |dPA| <= 10 deg."""
    return n > 0 and np.isfinite(median_dpa) and median_dpa <= MAX_MEDIAN_DPA


def convention_check(conv: Table, pipe: Table) -> dict:
    """Median |dPA| of rows matched within 0.2" with e >= 0.5 in both and S/N >= 20 in both."""

    def sel(t):
        with np.errstate(divide="ignore", invalid="ignore"):
            snr = np.asarray(t["isophotal_flux"], float) / np.asarray(
                t["isophotal_flux_err"], float
            )
        return (np.asarray(t["ellipticity"], float) >= 0.5) & (snr >= 20)

    a, b = conv[sel(conv)], pipe[sel(pipe)]
    idx, sep, _ = a["sky_centroid"].match_to_catalog_sky(b["sky_centroid"])
    m = sep < MATCH_ARCSEC * u.arcsec
    dpa = lensmodel.axis_offset_deg(
        np.asarray(a["sky_orientation"], float)[m],
        np.mod(np.asarray(b["sky_orientation"], float)[idx[m]], 180.0),
    )
    med = float(np.median(dpa)) if m.any() else float("nan")
    # the alternative convention (PA = 90 - theta) as a diagnostic only
    alt = lensmodel.axis_offset_deg(
        np.mod(180.0 - np.asarray(a["sky_orientation"], float)[m], 180.0),
        np.mod(np.asarray(b["sky_orientation"], float)[idx[m]], 180.0),
    )
    return {
        "n_matched": int(m.sum()),
        "median_abs_dpa_deg": round(med, 2),
        "fraction_within_20deg": round(float(np.mean(dpa <= 20)), 3) if m.any() else None,
        "median_abs_dpa_alternative_90_minus_theta_deg": (
            round(float(np.median(alt)), 2) if m.any() else None
        ),
        "passes": bool(convention_passes(med, int(m.sum()))),
    }


def posthoc(out: Path) -> dict:
    """Post hoc (not in PLAN.md; never changes the frozen outcome): strong-shear background
    classes split at z_phot > 15 (eazy fits at the grid edge)."""
    res: dict = {"note": "post hoc diagnostic (not in the plan); does not change the outcome"}
    for model in (arcs_rxj2129.PRIMARY, arcs_rxj2129.SECOND):
        t = Table.read(out / model / "arcs.ecsv")
        b = t[t["strong_shear"] & (t["population"] == "background")]
        bad = np.asarray(b["z_phot"], float) > 15
        res[model] = {
            name: {
                "n": int(s.sum()),
                "aligned": int((b["orientation_class"][s] == "aligned").sum()),
                "anti": int((b["orientation_class"][s] == "anti").sum()),
            }
            for name, s in (("z_phot>15", bad), ("z_phot<=15", ~bad))
        }
    return res


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--posthoc-only", action="store_true", help="re-write posthoc.json only")
    ap.add_argument("--dja", type=Path)
    ap.add_argument("--pipeline", type=Path)
    ap.add_argument("--photoz", type=Path)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--work", type=Path, help="scratch dir for the converted catalogue")
    args = ap.parse_args(argv)
    if args.posthoc_only:
        (args.out / "posthoc.json").write_text(json.dumps(posthoc(args.out), indent=1) + "\n")
        return 0
    if None in (args.dja, args.pipeline, args.photoz, args.work):
        ap.error("--dja, --pipeline, --photoz and --work are required for a run")
    pipe = Table.read(args.pipeline)
    px = pipeline_px_arcsec(pipe)
    conv = to_pipeline(Table.read(args.dja), px)
    check = convention_check(conv, pipe)
    check["pipeline_px_arcsec"] = round(px, 5)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "convention.json").write_text(json.dumps(check, indent=1) + "\n")
    print(json.dumps(check))
    if not check["passes"]:
        print("PA convention check fails: no arcs run (PLAN.md)")
        return 1
    args.work.mkdir(parents=True, exist_ok=True)
    cat = args.work / "rxj2129-dja-v7.0-as-pipeline_cat.ecsv"
    conv.write(cat, overwrite=True)
    arcs_rxj2129.main(["--catalog", str(cat), "--photoz", str(args.photoz), "--out", str(args.out)])
    summary = json.loads((args.out / "summary.json").read_text())
    summary["plan"] = "results/arcs/rxj2129-dja/PLAN.md"
    summary["inputs"]["dja_catalog"] = args.dja.name
    summary["inputs"]["dja_catalog_sha256"] = sha256_file(args.dja)
    summary["inputs"]["pipeline_catalog_sha256"] = sha256_file(args.pipeline)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    (args.out / "posthoc.json").write_text(json.dumps(posthoc(args.out), indent=1) + "\n")
    return 0


if __name__ == "__main__":
    np.seterr(invalid="ignore")
    raise SystemExit(main())
