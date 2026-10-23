"""El Gordo: gate K of the reference-catalogue extra-image test (stage 1).

Frozen plan: ``results/lens_consistency/elgordo_djaref/PLAN_1.md``. Rule 1 runs gate K of
``m0417_images`` on the 56 spectroscopic Caminha+2023 images (``elgordo-caminha23``, both model and
images in the JWST frame, D-034 rule 1). Rule 2 reads only the FITS headers of the DJA v7.0
photometry (HTTP range) to check the stage-2 reference inputs. Writes ``summary.json``.

    python scripts/elgordo_djaref.py
"""

from __future__ import annotations

import json
from pathlib import Path

import lens_consistency as lc
import m0417_images as mi
import numpy as np
import requests
from astropy.io import fits

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
    r = requests.get(url, headers={"Range": f"bytes=0-{HEAD_BYTES - 1}"}, timeout=60)
    r.raise_for_status()
    total = int(r.headers.get("Content-Range", "/0").rsplit("/", 1)[-1] or 0)
    return r.content, total


def main() -> int:
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
