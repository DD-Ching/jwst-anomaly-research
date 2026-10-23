"""SMACS 0723: the reference-catalogue extra-image test (stage 1).

Frozen plan: ``results/lens_consistency/smacs_djaref/PLAN_1.md``. Rule 0 takes the systems whose
redshift ``input.par`` keeps fixed (no ``z_m_limit``) as spectroscopic; rule 1 runs gate K of
``m0417_images`` on them (``smacs0723-iclv2``, no frame offset); rule 2 reads only the FITS headers
of the DJA v7.4 photometry (HTTP range) and checks the columns El Gordo stage 2b used. Writes
``summary.json``.

    python scripts/smacs_djaref.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import elgordo_djaref as eg
import lens_consistency as lc
import m0417_images as mi
import numpy as np

from jwst_anomaly import lensmodel

FIELD = "smacs0723"
MODEL = "smacs0723-iclv2"
OUT = Path("results/lens_consistency/smacs_djaref")
IMAGES = "mahler22_images.ecsv"
DJA_PHOT = "https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/smacs0723-grizli-v7.4-fix_phot.fits"
GRID_STEP = 0.25  # m0417_images default (PLAN_1.md rule 1)
# PLAN_1.md rule 2: the columns El Gordo stage 2b read (D-034 addenda p-q)
STAGE2_COLS = (
    *eg.SHAPE_COLS,
    "flux_auto",
    "mag_auto",
    "tot_corr",
    "flux_aper_1",
    "f200w_flux_aper_1",
    "f200w_fluxerr_aper_1",
    "f200w_mask_aper_1",
    "f200w_tot_corr",
)
STAGE2_CARDS = ("APER_1",)


def fixed_z_systems(input_par: str, systems) -> list[str]:
    """PLAN_1.md rule 0: systems without a ``z_m_limit`` line in ``input.par``."""
    free = set()
    for line in input_par.splitlines():
        f = line.split()
        if f and f[0] == "z_m_limit" and len(f) > 2:
            free.add(str(int(float(f[2]))))
    return sorted({str(s) for s in systems} - free, key=lambda s: (len(s), s))


def stage2_inputs(blob: bytes, total_bytes: int) -> dict:
    """PLAN_1.md rule 2 on the header bytes of the DJA photometry table."""
    heads = eg.header_cards(blob)
    tab = next((h for h in heads if h.get("XTENSION") == "BINTABLE"), None)
    if tab is None:
        return {"reachable": True, "bytes": total_bytes, "table_header": False, "passes": False}
    cols = {str(tab[f"TTYPE{i}"]).lower() for i in range(1, tab["TFIELDS"] + 1)}
    missing = [c for c in STAGE2_COLS if c not in cols]
    missing += [k for k in STAGE2_CARDS if k not in tab]
    return {
        "reachable": True,
        "bytes": total_bytes,
        "rows": int(tab["NAXIS2"]),
        "n_columns": len(cols),
        "missing": missing,
        "aper_1_pix": float(tab["APER_1"]) if "APER_1" in tab else None,
        "passes": not missing,
    }


def write_images(out: Path = OUT) -> tuple[Path, list[str]]:
    """The model's ``arcs.dat`` with ``z_spec`` from rule 0."""
    _, files, _ = lc.load_model(MODEL)
    t = lensmodel.load_lenstool_images(files["arcs.dat"])
    fixed = fixed_z_systems(Path(files["input.par"]).read_text(), t["system"])
    t["z_spec"] = np.isin(np.asarray(t["system"], str), fixed)
    t.meta["provenance"] = "observed"
    t.meta["source"] = (
        f"{MODEL} arcs.dat (sha256 in lensmodel.SMACS0723_MAHLER22_ICLV2); z_spec = systems "
        "without z_m_limit in input.par (PLAN_1.md rule 0, ASSUMPTION)"
    )
    out.mkdir(parents=True, exist_ok=True)
    path = out / IMAGES
    t.write(path, format="ascii.ecsv", overwrite=True)
    return path, fixed


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.parse_args(argv)
    _, fixed = write_images()
    # fitted_z's grid starts at z_lens; take it from best.par (0.39), as lc.MODELS has none
    lc.MODELS[MODEL].setdefault("z_lens", lc.load_model(MODEL)[0].z_lens)
    mi.FIELDS[FIELD] = {"images": IMAGES, "models": (MODEL,)}
    summary = {"stage1_fixed_z_systems": fixed}
    summary["stage1_gate_K"] = mi.run_field(FIELD, OUT, GRID_STEP)
    blob, total = eg.fetch_head(DJA_PHOT)
    summary["stage1_reference_inputs"] = stage2_inputs(blob, total) | {"url": DJA_PHOT}
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
