"""Abell 2744: Bergamini+2023 and MARS vs UNCOVER v2.0 images absent from their catalogues.

Frozen plan: ``results/lens_consistency/a2744_heldout/PLAN.md``. Rules 1-5 of the
``mars_hff_heldout`` stage, unchanged, per model, with that model's own catalogued images as
the "seen" list. MARS is reported only (its training gate K failed). No flag, flux or
extra-image reading. ``posthoc.json`` (not in the plan) splits the Bergamini residuals
by region.

    python scripts/a2744_heldout.py --zip DIR [--out DIR]

``DIR`` receives the MARS zip (md5 checked); the zip is deleted after the run unless ``--keep``.
"""

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

import m1149_hst_heldout as hh
import mars_hff_heldout as mo
import mars_hff_inputs as mh
import numpy as np
import requests
from astropy.table import Table
from m0417_canucsref_posthoc import _clean

OUT = Path("results/lens_consistency/a2744_heldout")
IMAGES = Path("results/lens_consistency/abell2744_uncoverref/uncover2_images_spec.ecsv")
MODEL = "abell2744-bergamini23"
MARS_FIELD = "Abell 2744"
MARS_GATE_K_RMS = 1.112  # D-034 addendum ao (fails; MARS is reported only)
NORTH_DEC = -30.385  # post hoc: images north of this Dec are near the northern subcluster


def shared_over_3(blocks: dict) -> list[str]:
    """PLAN.md readings: held-out images > 3 arcsec in every model."""
    sets = [set(b["over_3arcsec"]) for b in blocks.values()]
    return sorted(set.intersection(*sets)) if sets else []


def _stats(d: np.ndarray) -> dict:
    """rms over the matched images and the fraction within 3 arcsec; no prediction (inf) is a
    miss, as in ``gate_h``."""
    fin = d[np.isfinite(d)]
    return {
        "n": len(d),
        "rms_arcsec": round(float(np.sqrt(np.mean(fin**2))), 3) if len(fin) else None,
        "within_3arcsec": round(float(np.mean(np.isfinite(d) & (d <= 3.0))), 3) if len(d) else None,
    }


def _min_sep_arcsec(t: Table, a: np.ndarray, b: np.ndarray) -> float | None:
    """Smallest separation (arcsec) between the images in mask ``a`` and those in mask ``b``."""
    if not (a.any() and b.any()):
        return None
    _, _, d = mo._offsets_arcsec(t["ra"][a], t["dec"][a], t["ra"][b], t["dec"][b])
    return round(float(d.min()), 2)


def mars_vs_bergamini(out: Path) -> list[dict]:
    """PLAN.md reading (from the stored tables): residuals of the images both models tested."""
    rb = Table.read(out / f"heldout_{MODEL}.ecsv")
    rm = Table.read(out / "heldout_A2744.ecsv")
    b = dict(
        zip(np.asarray(rb["image_id"], str), np.asarray(rb["dtheta_arcsec"], float), strict=True)
    )
    rows = []
    for i, d, inside in zip(rm["image_id"], rm["dtheta_arcsec"], rm["inside"], strict=True):
        if bool(inside) and str(i) in b:
            rows.append(
                {
                    "image_id": str(i),
                    "bergamini_arcsec": round(b[str(i)], 2),
                    "mars_arcsec": round(float(d), 2),
                }
            )
    return rows


def posthoc(out: Path = OUT, images: Path = IMAGES) -> dict:
    """Post hoc (not in the plan): Bergamini held-out residuals split by region, by z within
    the northern group, and without family 80 (18 images of one source)."""
    r = Table.read(out / f"heldout_{MODEL}.ecsv")
    r = r[np.asarray(r["inside"], bool)]
    t = Table.read(images)
    row = {str(i): k for k, i in enumerate(t["image_id"])}
    t = t[[row[str(i)] for i in r["image_id"]]]  # aligned with r
    d = np.asarray(r["dtheta_arcsec"], float)
    north = np.asarray(t["dec"], float) > NORTH_DEC
    top = np.array([str(s).split(".")[0] for s in r["system"]])
    hiz = np.asarray(r["z"], float) >= 3.0
    return {
        "label": "post hoc; not in the frozen plan",
        "north_dec_cut": NORTH_DEC,
        "north": _stats(d[north]) | {"systems": sorted({str(x) for x in top[north]})},
        "north_z_lt_3": _stats(d[north & ~hiz])
        | {"systems": sorted({str(x) for x in top[north & ~hiz]})},
        "north_z_ge_3": _stats(d[north & hiz])
        | {"systems": sorted({str(x) for x in top[north & hiz]})},
        "core": _stats(d[~north]) | {"systems": sorted({str(x) for x in top[~north]})},
        "without_system_80": _stats(d[top != "80"]),
        "min_sep_north_z_ge_3_to_z_lt_3_arcsec": _min_sep_arcsec(t, north & hiz, north & ~hiz),
    }


def run(zdir: Path, out: Path = OUT) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    images = Table.read(IMAGES)
    blocks = {MODEL: hh.model_block(MODEL, images, out)}
    meta = mh.fetch_zip(requests.Session(), zdir)
    with zipfile.ZipFile(zdir / mh.ZIP) as zf:
        cat = mh.parse_catalogue(zf.read(mh.CATALOGUE).decode())
        blocks["mars-A2744"] = mo.field_block(
            zf, MARS_FIELD, cat, out, images=(IMAGES, MARS_GATE_K_RMS)
        )
    blocks["mars-A2744"]["role"] = "reported only (training gate K failed)"
    # field_block names its input list after CANUCS; here it is UNCOVER v2.0
    blocks["mars-A2744"]["reference_images"] = blocks["mars-A2744"].pop("canucs_images")
    summary = {
        "plan": str(out / "PLAN.md"),
        "images": str(IMAGES),
        "n_images": len(images),
        "zip": meta,
        "provenance": "positions and z: observed (UNCOVER v2.0; model image lists); "
        "residuals: model_prediction; thresholds: ASSUMPTION (mars_hff_heldout/PLAN.md)",
        "models": blocks,
        "over_3arcsec_in_both": shared_over_3(blocks),
    }
    (out / "summary.json").write_text(json.dumps(_clean(summary), indent=1, allow_nan=False) + "\n")
    p = posthoc(out) | {"mars_vs_bergamini (PLAN.md reading)": mars_vs_bergamini(out)}
    (out / "posthoc.json").write_text(json.dumps(p, indent=1) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--zip", type=Path, required=True, help="directory for the MARS zip")
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--keep", action="store_true", help="keep the downloaded zip")
    a = ap.parse_args(argv)
    try:
        s = run(a.zip, a.out)
    finally:
        if not a.keep:
            (a.zip / mh.ZIP).unlink(missing_ok=True)
    keys = ("frame_offset", "seen", "held_out", "types", "outside_map", "gate_H", "outcome")
    brief = {m: {k: b.get(k) for k in keys} for m, b in s["models"].items()}
    print(json.dumps(_clean(brief | {"shared": s["over_3arcsec_in_both"]}), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
