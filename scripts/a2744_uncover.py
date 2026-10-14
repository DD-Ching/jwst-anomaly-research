"""Abell 2744: the UNCOVER v2.0 lens model (Furtak+2023, Price+2025) as an independent second model.

Step 1 of the UNCOVER cross-check (TASKS.md "Supporting" item 1): validate the published v2.0
best-fit deflection maps (0.1"/pix, D_LS/D_S = 1) before any counter-image test uses them.

* **Known cases (own list):** UNCOVER's own multiple-image table (``A2744_multiple-images_v2.0``,
  ids ``S.KI``: source ``S.K``, image ``I``), every image at its catalogued redshift, traced
  through ``lensmodel.MapLensModel`` and re-imaged (``imageplane_residuals``). The release quotes
  a lens-plane rms of 0.60"; ours must come within ``lens_consistency.IMAGE_LIST_RMS_FACTOR`` of
  it (ASSUMPTION, as for the HFF maps, D-035) before the model is used.
* **Cross-model (Bergamini+23 list):** the Bergamini+23 ``arcs.dat`` images with their
  ``best.par`` redshifts through the UNCOVER maps: per family the image-plane residuals and the
  number of predicted images against the number catalogued. Families with predicted images that
  are not catalogued (``n_extra``) are where the two models disagree on multiplicity; those are
  the inputs of the D-041 cross-model rule in the next step.

Everything computed is ``model_prediction``; the image positions are ``observed``.

Data: the maps are fetched once from the UNCOVER Google Drive release (sha256-pinned; 216 MB
compressed in three files, ~250 MB on disk while the run lasts; delete them afterwards in cloud
sessions).

Usage: python scripts/a2744_uncover.py [--out results/lens_consistency/abell2744-uncover2]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402
from jwst_anomaly.photometry import fetch_catalog  # noqa: E402

_DRIVE = "https://drive.usercontent.google.com/download?export=download&confirm=t&id="
# UNCOVER DR4 lensing maps v2.0 (2024-08-07), Best-model_low-resolution_100mas; accessed 2026-10-14
UNCOVER_V2 = {
    "alpha_x": (
        _DRIVE + "1BHZ63x9nYAd5O0s1K_mO1UndyaekeJjz",
        "7ba107377e1585327d9ad95f8d847941dbcb7864611c442af6b712b455df4803",
    ),
    "alpha_y": (
        _DRIVE + "1Ht5KmqeXifqI8Q0tSRcVzPvEz6iAaH3w",
        "23469d560df560f6c7136cf273695aebe0fedd09853f256599b55c30335d156a",
    ),
    "images": (
        _DRIVE + "1uQQQBA1W_zd-y8SILZfBYKpVGAGoTxMG",
        "0d2649b3991129df0b775965caadf10109714d32cd804dcc76afa17c9c789119",
    ),
}
Z_LENS = 0.308
COSMOLOGY = (70.0, 0.3)  # ASSUMPTION: the Zitrin-code default; maps are at D_LS/D_S = 1 anyway
QUOTED_RMS = 0.60  # README v2.0, lens-plane rms (arcsec)
MATCH_ARCSEC = 1.5  # a predicted image within this of a catalogued one is that image (ASSUMPTION)


def uncover_images(path: Path) -> Table:
    """UNCOVER's multiple-image table as a ``load_lenstool_images``-like table."""
    t = Table.read(path)
    ids = [f"{float(i):.2f}" for i in t["ID"]]
    out = Table(
        {
            "image_id": ids,
            "system": [i[:-1].rstrip(".") for i in ids],  # "1.11" -> source "1.1"
            "ra": np.asarray(t["RA"], float),
            "dec": np.asarray(t["Dec"], float),
            "a": np.full(len(t), np.nan),
            "z": np.asarray(t["z"], float),
        }
    )
    out.meta.update(provenance="observed", source=str(path))
    return out


def multiplicity(model, grid, bt: Table) -> Table:
    """Per family: catalogued and predicted image counts, and predicted images more than
    ``MATCH_ARCSEC`` from every catalogued image (``n_extra``; their |mu| listed)."""
    rows = []
    for fam in np.unique(np.asarray(bt["system"]).astype(str)):
        sel = bt[(np.asarray(bt["system"]).astype(str) == fam) & np.isfinite(bt["beta_x"])]
        if len(sel) < 2:
            continue
        z_s = float(sel["z_used"][0])
        pred = lensmodel.find_images(
            model, grid, float(np.mean(sel["beta_x"])), float(np.mean(sel["beta_y"])), z_s
        )
        ox, oy = model.to_frame(np.asarray(sel["ra"], float), np.asarray(sel["dec"], float))
        if len(pred):
            d = np.hypot(pred["x"][:, None] - ox[None, :], pred["y"][:, None] - oy[None, :])
            extra = np.flatnonzero(d.min(axis=1) > MATCH_ARCSEC)
        else:
            extra = np.array([], int)
        rows.append(
            (
                fam,
                z_s,
                len(sel),
                len(pred),
                len(extra),
                ",".join(f"{abs(float(pred['magnification'][k])):.1f}" for k in extra),
                ",".join(f"{float(pred['ra'][k]):.5f}:{float(pred['dec'][k]):.5f}" for k in extra),
            )
        )
    return Table(
        rows=rows or None,
        names=(
            "system",
            "z_used",
            "n_catalogued",
            "n_predicted",
            "n_extra",
            "extra_mu",
            "extra_radec",
        ),
        dtype=(str, float, int, int, int, str, str),
    )


def run(model, images: Table, zml: dict, label: str, out: Path) -> dict:
    bt = lensmodel.backtrace_images(model, images, zml)
    traced = images[np.isfinite(bt["beta_x"])]
    half = lc.grid_half_width(model, traced)
    grid = lensmodel.DeflectionGrid.compute(model, half, 0.25)
    ip = lensmodel.imageplane_residuals(model, grid, bt)
    d = np.asarray(ip["dtheta_arcsec"], float)
    ok = np.isfinite(d)
    mult = multiplicity(model, grid, bt)
    for t in (ip, mult):
        t.meta.update(provenance="model_prediction", source=f"abell2744-uncover2 vs {label}")
    ip.write(out / f"imageplane_{label}.ecsv", overwrite=True)
    mult.write(out / f"multiplicity_{label}.ecsv", overwrite=True)
    return {
        "image_list": label,
        "n_images": len(images),
        "n_traced": int(np.isfinite(bt["beta_x"]).sum()),
        "n_reimaged": int(ok.sum()),
        "rms_dtheta_arcsec": float(np.sqrt(np.mean(d[ok] ** 2))) if ok.any() else None,
        "median_dtheta_arcsec": float(np.median(d[ok])) if ok.any() else None,
        "max_dtheta_arcsec": float(d[ok].max()) if ok.any() else None,
        "n_over_2arcsec": int((d[ok] > 2.0).sum()),
        "families": len(mult),
        "families_with_extra": int((mult["n_extra"] > 0).sum()) if len(mult) else 0,
        "families_short": int((mult["n_predicted"] < mult["n_catalogued"]).sum())
        if len(mult)
        else 0,
    }


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=Path("results/lens_consistency/abell2744-uncover2"))
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    files = {k: fetch_catalog(url, sha) for k, (url, sha) in UNCOVER_V2.items()}
    h0, om = COSMOLOGY
    model = lensmodel.MapLensModel.from_fits(
        files["alpha_x"],
        files["alpha_y"],
        Z_LENS,
        FlatLambdaCDM(H0=h0, Om0=om),
        source="abell2744-uncover2",
    )
    summary = {
        "model": "abell2744-uncover2",
        "model_sha256": model.sha256,
        "quoted_rms_arcsec": QUOTED_RMS,
    }
    own = run(model, uncover_images(files["images"]), {}, "uncover2", args.out)
    own["gate_passes"] = bool(
        own["rms_dtheta_arcsec"] is not None
        and own["rms_dtheta_arcsec"] <= lc.IMAGE_LIST_RMS_FACTOR * QUOTED_RMS
    )
    summary["own_list"] = own
    _, bfiles, par = lc.load_model("abell2744-bergamini23")
    berg = lensmodel.load_lenstool_images(bfiles["arcs.dat"])
    summary["bergamini23_list"] = run(model, berg, par["z_m_limit"], "bergamini23", args.out)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return summary


if __name__ == "__main__":
    main()
