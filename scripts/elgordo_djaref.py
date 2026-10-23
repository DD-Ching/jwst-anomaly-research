"""El Gordo: the reference-catalogue extra-image test (stages 1 and 2).

Frozen plans: ``results/lens_consistency/elgordo_djaref/PLAN_1.md`` (stage 1): rule 1 runs gate K
of ``m0417_images`` on the 56 spectroscopic Caminha+2023 images (``elgordo-caminha23``, both model
and images in the JWST frame, D-034 rule 1); rule 2 reads only the FITS headers of the DJA v7.0
photometry (HTTP range) to check the stage-2 reference inputs. ``PLAN_2.md`` (``--stage 2``):
predicted / control tables on GO 1176 F200W, the DJA reference with grizli total fluxes, local
noise and gate C2 shared, as Abell 2744 stage 2b. Writes ``summary.json``.

    python scripts/elgordo_djaref.py [--stage 2] [--cache data/cache/elgordox]
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

import a2744_uncoverref as uref
import lens_consistency as lc
import m0417_canucsref as cr
import m0417_extra as mx
import m0417_images as mi
import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from m0417_canucsref_posthoc import _clean

from jwst_anomaly import lensmodel

FIELD = "elgordo"
MODEL = "elgordo-caminha23"
OUT = Path("results/lens_consistency/elgordo_djaref")
IMAGES = "caminha23_images_jwstframe.ecsv"
DJA_PHOT = "https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/elgordo-grizli-v7.0-fix_phot.fits"
SHAPE_COLS = ("x_image", "y_image", "a_image", "b_image", "theta_image")
GRID_STEP = 0.25  # m0417_images default (PLAN_1.md rule 1)
APER_KEYS = ("APER_0", "APER_1", "APER_2", "APER_3")  # diameters in pixels
HEAD_BYTES = 2_000_000  # primary + table header blocks (cards only; no rows are read)
# stage 2 (PLAN_2.md; ASSUMPTIONs)
OBS, TILE = (
    "o241",
    "t012",
)  # GO 1176 F200W tile containing all 56 images (PLAN_2.md footprint check)
S3_1176 = "s3://stpubdata/jwst/public/jw01176/L3/t/{0}/jw01176-{0}_{1}_nircam_clear-f200w_{2}"
FRAME = OUT / "frame"
BAND, APER_IX = "f200w", 1  # 0.5-arcsec aperture (rule 2)
MASK_FRAC = 0.2  # grizli: aperture invalid when the masked area exceeds this fraction
UNITS_JY = {23.9: 1e-6, 28.9: 1e-8}  # AB zero point of the flux unit -> Jy (rule 2a)
UNIT_TOL_MAG = 0.05
ROT_MAX_DEG = 1.0  # rule 2b


def write_images(out: Path = OUT) -> Path:
    """The model's ``arcs.dat`` in the JWST frame, ``z_spec`` = true (all 56 spectroscopic)."""
    _, files, _ = lc.load_model(MODEL)
    t = lensmodel.load_lenstool_images(files["arcs.dat"])
    lc.shift_images(MODEL, t)
    t["z_spec"] = np.ones(len(t), bool)
    t.meta["provenance"] = "observed"
    t.meta["source"] = (
        f"{MODEL} arcs.dat (sha256 in lensmodel.ELGORDO_CAMINHA23), shifted by "
        f"frame_offset_arcsec {lc.MODELS[MODEL]['frame_offset_arcsec']} (D-034 rule 1)"
    )
    out.mkdir(parents=True, exist_ok=True)
    path = out / IMAGES
    t.write(path, format="ascii.ecsv", overwrite=True)
    return path


def header_cards(blob: bytes) -> list[fits.Header]:
    """Every complete FITS header in ``blob`` (2880-byte blocks), skipping the data units."""
    heads, pos = [], 0
    while pos < len(blob):
        end = blob.find(b"END" + b" " * 77, pos)
        if end < 0 or (end - pos) % 80:
            break
        stop = pos + ((end - pos) // 2880 + 1) * 2880
        h = fits.Header.fromstring(blob[pos:stop].decode("ascii"))
        heads.append(h)
        size = abs(h.get("BITPIX", 8)) // 8 * h.get("NAXIS1", 0) * h.get("NAXIS2", 0)
        size += h.get("PCOUNT", 0)
        pos = stop + -(-size // 2880) * 2880
    return heads


def reference_inputs(blob: bytes, total_bytes: int) -> dict:
    """PLAN_1.md rule 2 on the header bytes of the DJA photometry table."""
    heads = header_cards(blob)
    tab = next((h for h in heads if h.get("XTENSION") == "BINTABLE"), None)
    if tab is None:
        return {"reachable": True, "bytes": total_bytes, "table_header": False, "passes": False}
    cols = [str(tab[f"TTYPE{i}"]).lower() for i in range(1, tab["TFIELDS"] + 1)]
    tot = sorted(c for c in cols if c.endswith("_tot_1") or c.endswith("_tot_0"))
    has_shape = all(c in cols for c in SHAPE_COLS)
    return {
        "reachable": True,
        "bytes": total_bytes,
        "rows": int(tab["NAXIS2"]),
        "n_columns": len(cols),
        "shape_columns_present": has_shape,
        "total_flux_columns": tot,
        "passes": bool(has_shape and tot),
        # post hoc (read after rule 2 failed): the table carries aperture fluxes with a per-band
        # aperture-to-total correction instead of total-flux columns
        "post_hoc_tot_corr_columns": sorted(c for c in cols if c.endswith("_tot_corr")),
        "post_hoc_aperture_flux_columns_f200w": sorted(
            c for c in cols if c.startswith("f200w_flux_aper_")
        ),
        "post_hoc_detection_aperture_pix": {k: float(tab[k]) for k in APER_KEYS if k in tab},
    }


def fetch_head(url: str = DJA_PHOT) -> tuple[bytes, int]:
    """The first ``HEAD_BYTES`` of ``url`` and its total size (-1 when not stated); a server that
    ignores the range is refused, so the whole table is never read."""
    with requests.get(
        url, headers={"Range": f"bytes=0-{HEAD_BYTES - 1}"}, timeout=60, stream=True
    ) as r:
        r.raise_for_status()
        if r.status_code != 206:
            raise RuntimeError(f"{url}: range request not honoured (HTTP {r.status_code})")
        blob = r.raw.read(HEAD_BYTES, decode_content=True)
        total = r.headers.get("Content-Range", "").rsplit("/", 1)[-1]
    return blob, int(total) if total.isdigit() else -1


def flux_unit_jy(t) -> dict:
    """PLAN_2.md rule 2a: the flux unit from ``mag_auto + 2.5 log10(flux_auto)``."""
    f = np.ma.filled(np.ma.asarray(t["flux_auto"], float), np.nan)
    m = np.ma.filled(np.ma.asarray(t["mag_auto"], float), np.nan)
    ok = np.isfinite(f) & (f > 0) & np.isfinite(m)
    u = float(np.median(m[ok] + 2.5 * np.log10(f[ok]))) if ok.any() else float("nan")
    zp = next((z for z in UNITS_JY if abs(u - z) <= UNIT_TOL_MAG), None)
    return {
        "zero_point_mag": u,
        "n_rows": int(ok.sum()),
        "unit_jy": UNITS_JY[zp] if zp else None,
        "passes": zp is not None,
    }


def grid_rotation_deg(t) -> dict:
    """PLAN_2.md rule 2b: affine (x_image, y_image) -> (xi, eta) fit; north-up, east-left check."""
    x, y = np.asarray(t["x_image"], float), np.asarray(t["y_image"], float)
    ra, dec = np.asarray(t["ra"], float), np.asarray(t["dec"], float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(ra) & np.isfinite(dec)
    x, y, ra, dec = x[ok], y[ok], ra[ok], dec[ok]
    d0, r0 = np.median(dec), np.median(ra)
    xi = (ra - r0 + 180.0) % 360.0 - 180.0
    xi = xi * np.cos(np.radians(d0)) * 3600.0
    eta = (dec - d0) * 3600.0
    a = np.column_stack([x, y, np.ones_like(x)])
    (dxi_dx, dxi_dy, _), *_ = np.linalg.lstsq(a, xi, rcond=None)
    (deta_dx, deta_dy, _), *_ = np.linalg.lstsq(a, eta, rcond=None)
    det = dxi_dx * deta_dy - dxi_dy * deta_dx
    rot = float(np.degrees(np.arctan2(dxi_dy, deta_dy)))
    return {
        "rotation_deg": rot,
        "determinant": float(det),
        "pixel_scale_arcsec": float(np.sqrt(abs(det))),
        "n_rows": int(ok.sum()),
        "passes": bool(abs(rot) <= ROT_MAX_DEG and det < 0),
    }


def to_cat_dja(t, unit_jy: float, pix: float, band: str = BAND, k: int = APER_IX) -> dict:
    """PLAN_2.md rule 2: the DJA table in ``m0417_arcref``'s catalogue form, grizli total flux
    ``{band}_flux_aper_k x flux_auto / flux_aper_k x {band}_tot_corr``."""

    def col(name):
        return np.ma.filled(np.ma.asarray(t[name], float), np.nan)

    fa, fe = col(f"{band}_flux_aper_{k}"), col(f"{band}_fluxerr_aper_{k}")
    det_aper = col(f"flux_aper_{k}")
    area = np.pi * (float(t.meta[f"APER_{k}"]) / 2.0) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        flux = fa * col("flux_auto") / det_aper * col(f"{band}_tot_corr") * unit_jy
        snr = fa / fe
    keep = (
        np.isfinite(flux)
        & (flux > 0)
        & (det_aper > 0)
        & ~(col(f"{band}_mask_aper_{k}") > MASK_FRAC * area)
    )
    t, flux, snr = t[keep], flux[keep], snr[keep]
    a, b = np.asarray(t["a_image"], float), np.asarray(t["b_image"], float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ell = 1.0 - b / a
    return {
        "label": np.asarray(t["id"]),
        "ra": np.asarray(t["ra"], float),
        "dec": np.asarray(t["dec"], float),
        "a": a * pix,
        "b": b * pix,
        "pa": cr.sky_pa(np.degrees(np.asarray(t["theta_image"], float))),
        "flux": flux,
        "snr": snr,
        "ext": np.ones(len(t), bool),
        "ell": ell,
        "scale": pix,
    }


def register_2() -> None:
    """Point ``m0417_extra``'s field table at this field, tile and directory."""
    mx.FIELDS[FIELD] = (IMAGES, MODEL, OBS, TILE)
    mx.IMAGES = OUT
    mx.OUT = FRAME
    mx.S3 = S3_1176


def fetch_table(url: str = DJA_PHOT) -> tuple[Table, str, int]:
    """The whole DJA photometry table in memory, with its sha256 and size."""
    r = requests.get(url, timeout=300)
    r.raise_for_status()
    blob = r.content
    with fits.open(io.BytesIO(blob)) as hl:
        t = Table(hl[1].data)
        t.meta[f"APER_{APER_IX}"] = float(hl[1].header[f"APER_{APER_IX}"])
    return t, hashlib.sha256(blob).hexdigest(), len(blob)


def stage_2(cache: Path) -> dict:
    """PLAN_2.md: rules 2a / 2b, frame tables, DJA reference, local noise, gate C2 shared."""
    import s3fs

    t, sha, nbytes = fetch_table()
    unit, grid = flux_unit_jy(t), grid_rotation_deg(t)
    res: dict = {"dja": {"sha256": sha, "bytes": nbytes, "rows": len(t)}, "rule_2a": unit}
    res["rule_2b"] = grid
    if not (unit["passes"] and grid["passes"]):
        res["outcome"] = "stopped at rule 2a / 2b; nothing read"
        return res
    register_2()
    fs = s3fs.S3FileSystem(anon=True, default_block_size=64 * 1024)
    cache.mkdir(parents=True, exist_ok=True)
    shapes = cache / Path(S3_1176.format(OBS, TILE, "cat.ecsv")).name
    if not shapes.exists():
        fs.get(S3_1176.format(OBS, TILE, "cat.ecsv")[5:], str(shapes))
    FRAME.mkdir(parents=True, exist_ok=True)
    frame = mx.run_field(FIELD, cache)  # rule 1
    (FRAME / "summary.json").write_text(
        json.dumps(_clean({FIELD: frame}), indent=1, allow_nan=False) + "\n"
    )
    w, d, pix = uref._read_sci(fs, S3_1176.format(OBS, TILE, "i2d.fits")[5:])
    rc = to_cat_dja(t, unit["unit_jy"], grid["pixel_scale_arcsec"])
    inputs = {
        "reference": Path(DJA_PHOT).name,
        "reference_sha256": sha,
        "n_rows_kept": int(len(rc["flux"])),
        "flux": f"{BAND}_flux_aper_{APER_IX} x flux_auto / flux_aper_{APER_IX} x {BAND}_tot_corr",
        "noise_image": S3_1176.format(OBS, TILE, "i2d.fits"),
        "noise_pixel_scale_arcsec": pix,
    }
    depth = cr.depth_mag(rc["flux"], rc["snr"])
    src_x, src_y = (
        np.asarray(v, float) for v in w.world_to_pixel(SkyCoord(t["ra"], t["dec"], unit="deg"))
    )
    res["frame"] = frame
    for contain, out in ((False, OUT / "without_contain"), (True, OUT)):  # rule 3 / rule 7
        out.mkdir(parents=True, exist_ok=True)
        res["frozen" if contain else "without_contain"] = uref._reference_run(
            rc, depth, inputs, out, contain, d, w, src_x, src_y, pix, FIELD, FRAME
        )
    return res


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--stage", choices=("1", "2"), default="1")
    ap.add_argument("--cache", type=Path, default=Path("data/cache/elgordox"))
    args = ap.parse_args(argv)
    if args.stage == "2":
        summary = json.loads((OUT / "summary.json").read_text())
        summary["stage2"] = _clean(stage_2(args.cache))
        (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
        s2 = summary["stage2"]
        print(json.dumps({k: s2.get(k) for k in ("rule_2a", "rule_2b", "outcome")}, indent=1))
        if "frozen" in s2:
            print(json.dumps({k: s2["frozen"][k] for k in ("gate_C2_shared", "flags")}, indent=1))
        return 0
    write_images()
    mi.FIELDS[FIELD] = {"images": IMAGES, "models": (MODEL,)}
    summary = {"stage1_gate_K": mi.run_field(FIELD, OUT, GRID_STEP)}
    blob, total = fetch_head()
    summary["stage1_reference_inputs"] = reference_inputs(blob, total) | {"url": DJA_PHOT}
    summary["rules"] = {
        "K_MAX_RMS_ARCSEC": mi.K_MAX_RMS_ARCSEC,
        "K_MIN_MATCHED": mi.K_MIN_MATCHED,
        "MATCH_ARCSEC": mi.MATCH_ARCSEC,
        "grid_step_arcsec": GRID_STEP,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    g = summary["stage1_gate_K"][MODEL]["gate_K"]
    print(json.dumps({"gate_K": g, "ref": summary["stage1_reference_inputs"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
