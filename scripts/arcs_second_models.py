"""Predicted stretch directions of stored ``arcs`` rows in other published HFF models.

A model-error test for ``arcs`` anti rows.

For the strong-shear background rows of ``results/arcs/<model>/arcs.ecsv`` each HFF team model
given with ``--team`` (``<team>/<version>`` under the HFF HLSP tree, e.g. ``glafic/v4``) supplies
κ, γ1, γ2 maps at D_ls/D_s = 1, read by HTTP byte ranges (nearest pixel). They are scaled by
D_ls/D_s at the row's photo-z (and at z_s = 2); the stretching axis is the shear axis where
κ < 1 and the perpendicular where κ > 1
(as ``lensmodel``), and the pixel-frame angle is turned into a sky PA through the map WCS.
Classes use the thresholds stored in ``arcs.json``; rows without a background photo-z use
z_s = 2 only. Known-case check: rows the stored model
calls ``aligned`` must be aligned in each team model (that model's calibration, including its
γ sign convention).

Output: ``results/arcs/<model>/second_models.ecsv`` and ``second_models.json``
(``model_prediction``; observed axes from the stored table).

    python scripts/arcs_second_models.py --model abells1063-cats --field abells1063 \\
        --team glafic/v4 --team sharon/v4 --team keeton/v4
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from jwst_anomaly import lensmodel, paths

_spec = importlib.util.spec_from_file_location(
    "lens_consistency", Path(__file__).resolve().parent / "lens_consistency.py"
)
lc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(lc)

HFF = "https://archive.stsci.edu/pub/hlsp/frontier"


axis_offset = lensmodel.axis_offset_deg


def axis_class(off, aligned_deg: float, anti_deg: float) -> np.ndarray:
    """``aligned`` (≤ aligned_deg), ``anti`` (≥ anti_deg) or ``mixed`` (also for NaN offsets)."""
    off = np.asarray(off, float)
    return np.where(off <= aligned_deg, "aligned", np.where(off >= anti_deg, "anti", "mixed"))


def combine_classes(c_z2, c_zphot, use_zphot, inside) -> np.ndarray:
    """Team class: z2 and photo-z classes must agree; photo-z ignored where ``use_zphot`` is
    False (no usable background photo-z); ``off_map`` outside the team map."""
    c_zphot = np.where(use_zphot, c_zphot, c_z2)
    both = np.where(c_z2 == c_zphot, c_z2, "mixed")
    return np.where(inside, both, "off_map")


def stretch_angle_px(kappa, g1, g2) -> np.ndarray:
    """Pixel-frame stretching axis (deg from +x): shear axis if κ < 1, perpendicular if κ > 1."""
    phi = 0.5 * np.degrees(np.arctan2(g2, g1))
    return (phi + np.where(np.asarray(kappa) > 1.0, 90.0, 0.0)) % 180.0


def read_pixels(url: str, ra: np.ndarray, dec: np.ndarray) -> tuple[np.ndarray, WCS, np.ndarray]:
    """Nearest-pixel map values at (ra, dec) by HTTP byte ranges, the map WCS and an inside mask."""
    with fits.open(url, use_fsspec=True, fsspec_kwargs={"block_size": 2**16}) as h:
        hdu = next(x for x in h if x.header.get("NAXIS", 0) == 2)
        wcs = WCS(hdu.header).celestial
        ny, nx = hdu.header["NAXIS2"], hdu.header["NAXIS1"]
        x, y = wcs.world_to_pixel_values(ra, dec)
        xi, yi = np.round(x).astype(int), np.round(y).astype(int)
        inside = (xi >= 0) & (yi >= 0) & (xi < nx) & (yi < ny)
        vals = np.full(len(ra), np.nan)
        for k in np.flatnonzero(inside):
            vals[k] = float(hdu.section[yi[k], xi[k]])
    return vals, wcs, inside


def team_predictions(t: Table, field: str, team: str, z_cols: dict, scale) -> dict:
    name, ver = team.split("/")
    stem = f"{HFF}/{field}/models/{name}/{ver}/hlsp_frontier_model_{field}_{name}_{ver}"
    ra, dec = np.asarray(t["ra"], float), np.asarray(t["dec"], float)
    v, urls, grids = {}, {}, []
    for q in ("kappa", "gamma1", "gamma2"):
        urls[q] = f"{stem}_{q}.fits"
        v[q], wcs, inside = read_pixels(urls[q], ra, dec)
        grids.append((wcs.to_header_string(), tuple(inside)))
    if len(set(grids)) != 1:
        raise SystemExit(f"error: {team} kappa/gamma maps differ in WCS or footprint")
    x, y = wcs.world_to_pixel_values(ra, dec)
    out = {"urls": urls, "inside": inside}
    for key, z in z_cols.items():
        s = scale(z)
        k, g1, g2 = s * v["kappa"], s * v["gamma1"], s * v["gamma2"]
        ang = np.radians(stretch_angle_px(k, g1, g2))
        c1 = wcs.pixel_to_world(x, y)
        c2 = wcs.pixel_to_world(x + np.cos(ang), y + np.sin(ang))
        out[key] = {
            "pa": np.where(inside, c1.position_angle(c2).deg % 180.0, np.nan),
            "kappa": np.where(inside, k, np.nan),
            "g": np.where(inside, np.hypot(g1, g2) / np.abs(1.0 - k), np.nan),
        }
    return out


def run(args) -> dict:
    d = args.arcs_dir / args.model
    t = Table.read(d / "arcs.ecsv")
    t = t[np.asarray(t["strong_shear"], bool) & (t["population"] == "background")]
    stored = json.loads((d / "arcs.json").read_text())["assumptions"]
    if any(lc.MODELS[args.model].get("frame_offset_arcsec", (0.0, 0.0))):
        # team maps are in the HFF frame; positions here are JWST-frame (D-034, D-040)
        raise SystemExit("error: models with a frame_offset_arcsec are not supported yet")
    model, _, _ = lc.load_model(args.model)
    zp = np.asarray(t["z_phot"], float)
    use_zphot = np.isfinite(zp) & (zp > model.z_lens + stored["z_margin"])
    z_cols = {"z2": np.full(len(t), 2.0), "zphot": np.where(use_zphot, zp, 2.0)}
    summary = {
        "model": args.model,
        "label": "model_prediction (other teams' published maps); observed axes from arcs.ecsv",
        "assumptions": {k: stored[k] for k in ("aligned_deg", "anti_deg", "z_margin")},
        "teams": {},
    }
    out = t["label", "mag", "orientation_class", "pa_obs"]
    ctrl = np.asarray(t["orientation_class"] == "aligned")
    anti = np.asarray(t["orientation_class"] == "anti")
    for team in args.team:
        pred = team_predictions(t, args.field, team, z_cols, model.dls_ds)
        col = team.replace("/", "_")
        cls = {}
        for key in z_cols:
            off = axis_offset(t["pa_obs"], pred[key]["pa"])
            out[f"{col}_off_{key}"] = np.round(off, 1)
            out[f"{col}_kappa_{key}"] = np.round(pred[key]["kappa"], 3)
            out[f"{col}_g_{key}"] = np.round(pred[key]["g"], 2)
            cls[key] = axis_class(off, stored["aligned_deg"], stored["anti_deg"])
        both = combine_classes(cls["z2"], cls["zphot"], use_zphot, pred["inside"])
        out[f"{col}_class"] = both
        summary["teams"][team] = {
            "urls": pred["urls"],
            "control_aligned": {
                "n": int(ctrl.sum()),
                "aligned": int(np.sum(both[ctrl] == "aligned")),
                "anti": int(np.sum(both[ctrl] == "anti")),
            },
            "stored_anti": {
                str(int(lab)): str(c) for lab, c in zip(t["label"][anti], both[anti], strict=True)
            },
        }
    out.meta.update(provenance="model_prediction", source=f"{args.model} arcs.ecsv + HFF team maps")
    out.write(d / "second_models.ecsv", overwrite=True)
    (d / "second_models.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--model", required=True)
    p.add_argument("--field", required=True, help="HFF field name, e.g. abells1063")
    p.add_argument("--team", action="append", required=True, help="<team>/<version>")
    p.add_argument("--arcs-dir", type=Path, default=paths.repo_root() / "results" / "arcs")
    args = p.parse_args(argv)
    print(json.dumps(run(args), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
