"""MACS0417: extra-image test with CANUCS DR1 (BCG-subtracted, deblended) reference fluxes.

Frozen plan: ``results/lens_consistency/macs0417_canucsref/PLAN.md``. Runs ``m0417_arcref.run``
with its frozen rules unchanged and the JWST pipeline catalogue replaced by the CANUCS DR1
MACS0417 cluster-field photometry catalogue (column mapping in the plan).

    python scripts/m0417_canucsref.py --cat <path to the CANUCS DR1 photometry-cat.fits.gz>
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lens_consistency as lc
import m0417_arcref as ar
import numpy as np
from astropy.table import Table

from jwst_anomaly.provenance import file_sha256

OUT = Path("results/lens_consistency/macs0417_canucsref")
URL = (
    "https://archive.stsci.edu/hlsps/canucs/dr1/macs0417/clu/"
    "hlsp_canucs_jwst-hst_multi_macs0417-clu_multi_v1_photometry-cat.fits.gz"
)
BAND = "F200W"
# sky PA (East of North) = PA (CCW from +x) + 90 on a north-up image (ASSUMPTION)
PA_OFFSET_DEG = 90.0


def pixel_scale(x, y, ra, dec) -> float:
    """Arcsec per pixel from a linear fit of the catalogue's pixel to sky positions."""
    x, y, ra, dec = (np.asarray(v, float) for v in (x, y, ra, dec))
    d0 = np.median(dec)
    a = np.c_[x - np.median(x), y - np.median(y), np.ones_like(x)]
    e = np.linalg.lstsq(a, (ra - np.median(ra)) * np.cos(np.radians(d0)) * 3600.0, rcond=None)[0]
    n = np.linalg.lstsq(a, (dec - d0) * 3600.0, rcond=None)[0]
    return float(np.sqrt(abs(e[0] * n[1] - e[1] * n[0])))


def sky_pa(pa_image) -> np.ndarray:
    """Major-axis PA East of North (deg, 0-180) from the image angle CCW from +x (north up)."""
    return np.mod(np.asarray(pa_image, float) + PA_OFFSET_DEG, 180.0)


def depth_mag(flux_jy, snr) -> float:
    """``lens_consistency.depth_mag`` (median AB of rows with S/N 4-6) on this catalogue."""
    f, s = np.asarray(flux_jy, float), np.asarray(snr, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        mag = np.where(f > 0, -2.5 * np.log10(f / 3631.0), np.nan)
    return lc.depth_mag(Table({"snr": s, "mag": mag}))


def to_cat(t: Table, pa_offset: float = PA_OFFSET_DEG, pa_sign: float = 1.0) -> dict:
    """The CANUCS table in ``m0417_arcref``'s catalogue form (plan's column mapping)."""
    flux = np.asarray(t[f"FLUX_KRON_TOTAL_{BAND}"], float) * 1e-9
    keep = np.isfinite(flux)
    t = t[keep]
    flux = flux[keep]
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = flux / (np.asarray(t[f"FLUXERR_KRON_TOTAL_{BAND}"], float) * 1e-9)
    scale = pixel_scale(t["X"], t["Y"], t["RA"], t["DEC"])
    a = np.asarray(t["A"], float)
    b = np.asarray(t["B"], float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ell = 1.0 - b / a
    return {
        "label": np.asarray(t["SOURCE"]),
        "ra": np.asarray(t["RA"], float),
        "dec": np.asarray(t["DEC"], float),
        "a": a * scale,
        "b": b * scale,
        "pa": np.mod(pa_sign * np.asarray(t["PA"], float) + pa_offset, 180.0),
        "flux": flux,
        "snr": snr,
        "ext": ~np.asarray(t["FLAG_POINTSRC"], bool),
        "ell": ell,
        "scale": scale,
    }


def pa_sensitivity(path: Path) -> int:
    """Post hoc (not a frozen rule): gate C2 with the PA read in other conventions."""
    t = Table.read(path)
    res = {"label": "post hoc sensitivity; frozen convention is ccw_from_x_plus90"}
    conv = {
        "ccw_from_x_plus90": (PA_OFFSET_DEG, 1.0),
        "as_east_of_north": (0.0, 1.0),
        "cw_from_x_plus90": (PA_OFFSET_DEG, -1.0),
        "cw_from_x": (0.0, -1.0),
    }
    for name, (off, sign) in conv.items():
        cat = to_cat(t, off, sign)
        depth = depth_mag(cat["flux"], cat["snr"])
        tmp = OUT / "_pa_tmp"
        s = ar.run(Path("."), cat=cat, depth=depth, cat_inputs={}, out=tmp)
        res[name] = {
            "gate_C2": s["gate_C2"],
            "flags": s["flags"],
            "target_classes": s["target_classes"],
        }
        for f in tmp.iterdir():
            f.unlink()
        tmp.rmdir()
    (OUT / "pa_sensitivity.json").write_text(json.dumps(res, indent=1, allow_nan=False) + "\n")
    print(json.dumps(res, indent=1))
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cat", type=Path, required=True)
    ap.add_argument(
        "--pa-sensitivity",
        action="store_true",
        help="post hoc: gate C2 under other PA conventions, to pa_sensitivity.json",
    )
    args = ap.parse_args(argv)
    if args.pa_sensitivity:
        return pa_sensitivity(args.cat)
    cat = to_cat(Table.read(args.cat))
    inputs = {
        "canucs_cat": args.cat.name,
        "canucs_url": URL,
        "canucs_cat_sha256": file_sha256(args.cat),
        "n_rows_finite_flux": int(len(cat["flux"])),
        "flux_column": f"FLUX_KRON_TOTAL_{BAND}",
    }
    depth = depth_mag(cat["flux"], cat["snr"])
    summary = ar.run(Path("."), cat=cat, depth=depth, cat_inputs=inputs, out=OUT)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
