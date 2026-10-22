"""Abell 370: gate K of the BUFFALO Niemiec+2023 model A deflection maps (stage 1 only).

Frozen plan: ``results/lens_consistency/abell370_buffalo/PLAN.md``. Gate L checks the map
loading (κ from the deflection Hessian vs the published κ map); gate K runs ``m0417_images`` on
the Lagattuta et al. 2019 gold spectroscopic images. Writes ``summary.json`` to the plan's
directory. ``archive.stsci.edu`` needs a CA bundle with Mozilla's roots in cloud sessions
(docs/operations.md).

    python scripts/a370_buffalo.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lens_consistency as lc
import m0417_images as mi
import m1149_sharon as ms
import numpy as np
from astropy.table import Table
from m0417_canucsref_posthoc import _clean

from jwst_anomaly import lensmodel
from jwst_anomaly.photometry import fetch_catalog

FIELD = "abell370"
MODEL = "abell370-buffalo-a-maps"
OUT = Path("results/lens_consistency/abell370_buffalo")
IMAGES = "lagattuta19_gold_zspec.ecsv"
GRID_STEP_K = 0.25  # m0417_images default (PLAN.md rule 2)
L_MAX_MEDIAN_REL = ms.G1_MAX_MEDIAN_REL  # gate L (PLAN.md rule 1; ASSUMPTION)
L_MAP_STEP = 37  # κ sub-grid step in map pixels, as m0417_images.map_checks
FIELD_K2 = "abell370_k2"
IMAGES_K2 = "buffalo_sl_gold.ecsv"
BUFFALO_GOLD = (  # the model's own constraints (PLAN.md addendum K2)
    "https://archive.stsci.edu/hlsps/buffalo/abell370/catalogs/niemiec-lensing-dr1/"
    "hlsp_buffalo_hst_multi_abell370_multi_v1.0_sl-gold.dat",
    "cad4c02f52f509e1d6df77d0dff2e155c791a4decb5b0d761cc63dca7c5e711d",
)
LAGATTUTA19 = (
    "https://cdsarc.cds.unistra.fr/ftp/J/MNRAS/485/3738/table1.dat",
    "6c7e41483a179b07d43f062790b4e868bc29162019c59abe8ebf95290960cdc5",
)


def parse_table1(path: Path) -> Table:
    """CDS ``table1.dat`` of Lagattuta et al. 2019 by its byte-by-byte ReadMe (1-based columns)."""
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        line = line.ljust(55)
        zspec = line[26:32].strip()
        rows.append(
            (
                line[0:4].strip(),
                float(line[5:14]),
                float(line[15:25]),
                float(zspec) if zspec else np.nan,
                line[45:51].strip(),
                line[52:55].strip(),
            )
        )
    return Table(rows=rows, names=("image_id", "ra", "dec", "zspec", "class", "note"))


def select(t: Table) -> np.ndarray:
    """PLAN.md selection: gold, zspec present, Note not containing 1."""
    cls = np.asarray(t["class"], str)
    notes = np.asarray(t["note"], str)
    return (
        (cls == "gold")
        & np.isfinite(np.asarray(t["zspec"], float))
        & np.array(["1" not in n.split(",") for n in notes])
    )


def write_images(out: Path = OUT) -> Path:
    """The selected images in the ``m0417_images`` list format (z = zspec, fixed: ASSUMPTION)."""
    src = parse_table1(fetch_catalog(*LAGATTUTA19))
    sel = src[select(src)]
    t = Table({"image_id": sel["image_id"], "ra": sel["ra"], "dec": sel["dec"], "z": sel["zspec"]})
    t["ra"].unit = t["dec"].unit = "deg"
    t["z_spec"] = np.ones(len(t), bool)
    t.meta["provenance"] = "observed"
    t.meta["source"] = (
        f"{LAGATTUTA19[0]} (Lagattuta et al. 2019; gold, zspec present, Note without 1; "
        "z = zspec fixed, PLAN.md ASSUMPTION)"
    )
    out.mkdir(parents=True, exist_ok=True)
    path = out / IMAGES
    t.write(path, format="ascii.ecsv", overwrite=True)
    return path


def write_buffalo_images(out: Path = OUT) -> Path:
    """All rows of the BUFFALO gold list as released (z fixed: PLAN.md addendum K2 ASSUMPTION)."""
    src = lensmodel.load_lenstool_images(fetch_catalog(*BUFFALO_GOLD))
    z = np.asarray(src["z"], float)
    if not np.all(
        np.isfinite(z) & (z > 0)
    ):  # free / fitted rows would not be fixed spectroscopic z
        raise SystemExit(f"error: {BUFFALO_GOLD[0]}: rows without a fixed z")
    t = Table({k: src[k] for k in ("image_id", "ra", "dec", "z")})
    t["ra"].unit = t["dec"].unit = "deg"
    t["z_spec"] = np.ones(len(t), bool)
    t.meta["provenance"] = "observed"
    t.meta["source"] = f"{BUFFALO_GOLD[0]} (BUFFALO gold list; z fixed, PLAN.md addendum K2)"
    out.mkdir(parents=True, exist_ok=True)
    path = out / IMAGES_K2
    t.write(path, format="ascii.ecsv", overwrite=True)
    return path


def k1_overlap(out: Path = OUT, radius: float = 0.5) -> dict:
    """Post hoc (not in PLAN.md): K1 residuals split by whether the Lagattuta+2019 image has a
    BUFFALO gold image within ``radius`` arcsec (i.e. is one of the model's own constraints)."""
    k1 = Table.read(out / f"imageplane_{FIELD}_{MODEL}.ecsv")
    lag = Table.read(out / IMAGES)
    gold = Table.read(out / IMAGES_K2)
    cos = np.cos(np.deg2rad(float(np.mean(gold["dec"]))))
    dx = (np.asarray(lag["ra"])[:, None] - np.asarray(gold["ra"])[None, :]) * cos
    dy = np.asarray(lag["dec"])[:, None] - np.asarray(gold["dec"])[None, :]
    sep = np.hypot(dx, dy).min(axis=1) * 3600.0
    lag_ids = np.asarray(lag["image_id"], str)
    shared = dict(zip(lag_ids, sep <= radius, strict=True))
    ids = np.asarray(k1["image_id"], str)
    d = np.asarray(k1["dtheta_arcsec"], float)
    m = np.array([shared[i] for i in ids])
    out_d = {}
    for name, sel in (("shared", m), ("not_shared", ~m)):
        out_d[name] = {
            "n": int(sel.sum()),
            "n_finite": int(np.isfinite(d[sel]).sum()),
            "rms_arcsec": mi.rms_of(d[sel]),
            "within_3arcsec": int(np.sum(np.isfinite(d[sel]) & (d[sel] <= mi.MATCH_ARCSEC))),
            "image_ids": sorted(ids[sel].tolist()) if name == "not_shared" else None,
        }
    nearest = {i: round(float(v), 3) for i, v in zip(lag_ids, sep, strict=True)}
    return {"radius_arcsec": radius, **out_d, "nearest_gold_arcsec": nearest}


def k2_decision(k2: dict) -> str:
    """PLAN.md addendum K2 decision from ``run_field``'s summary."""
    if MODEL in k2["passing_models"]:
        return "passes: validated on its own constraints only; extra-image stages may follow"
    return "fails: the Abell 370 line stays stopped"


def gate_l(check: dict) -> dict:
    rel = check["median_rel_diff"]
    return {**check, "passes": bool(np.isfinite(rel) and rel <= L_MAX_MEDIAN_REL)}


def main(argv: list[str] | None = None) -> int:
    argparse.ArgumentParser(description=__doc__.split("\n\n")[0]).parse_args(argv)
    write_images()
    write_buffalo_images()
    model, _, _ = lc.load_model(MODEL)
    check = lc.map_check(model, lc.model_files(MODEL)["kappa_map"], "kappa", None, L_MAP_STEP)
    check["map"] = Path(check["map"]).name  # cache file name, not a local path
    summary: dict = {"gate_L": gate_l(check)}
    if not summary["gate_L"]["passes"]:
        summary["stopped"] = "gate L failed: map loading not verified (PLAN.md rule 1)"
    else:
        mi.FIELDS[FIELD] = {"images": IMAGES, "models": (MODEL,)}
        k = mi.run_field(FIELD, OUT, GRID_STEP_K)
        summary["stage1_gate_K"] = k
        if MODEL not in k["passing_models"]:
            # PLAN.md rule 2 as frozen; addendum K2 (a deviation) decides the model's further use
            summary["K1_stopped"] = "gate K1 failed: original route stopped (PLAN.md rule 2)"
        mi.FIELDS[FIELD_K2] = {"images": IMAGES_K2, "models": (MODEL,)}
        k2 = mi.run_field(FIELD_K2, OUT, GRID_STEP_K)
        summary["addendum_gate_K2"] = k2  # deviation (PLAN.md addendum K2)
        summary["post_hoc_k1_overlap"] = k1_overlap()
        summary["K2_decision"] = k2_decision(k2)
    summary["rules"] = {
        "L_MAX_MEDIAN_REL": L_MAX_MEDIAN_REL,
        "K_MAX_RMS_ARCSEC": mi.K_MAX_RMS_ARCSEC,
        "K_MIN_MATCHED": mi.K_MIN_MATCHED,
        "MATCH_ARCSEC": mi.MATCH_ARCSEC,
        "GRID_STEP_K": GRID_STEP_K,
    }
    summary = _clean(summary)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    gates = {
        g: summary.get(key, {}).get(MODEL, {}).get("gate_K")
        for g, key in (("K1", "stage1_gate_K"), ("K2", "addendum_gate_K2"))
    }
    print(json.dumps({"gate_L": summary["gate_L"]["median_rel_diff"], **gates}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
