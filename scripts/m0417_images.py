"""MACS0417 and MACS1423: published multiple images against the pinned CANUCS v1 maps.

Gate K per model on the spectroscopic images (image-plane rms and match fraction), then a per-image
reading: spectroscopic images at z_spec, the other systems at a per-model fitted z (grid scan).
Frozen plan: ``results/lens_consistency/macs0417_images/PLAN.md``.

    python scripts/m0417_images.py [--out results/lens_consistency/macs0417_images]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lens_consistency as lc
import m1149_sharon as ms
import numpy as np
from astropy.table import Table, vstack

from jwst_anomaly import lensmodel

OUT = Path("results/lens_consistency/macs0417_images")
# Frozen rules (PLAN.md; ASSUMPTIONs)
K_MAX_RMS_ARCSEC = 1.0  # gate K: image-plane rms of the spectroscopic images
K_MIN_MATCHED = 0.8  # gate K: fraction of spectroscopic images matched within MATCH_ARCSEC
MATCH_ARCSEC = 3.0  # an image is reproduced when a predicted image lies within this
Z_STEPS = 60  # log-spaced fitted-z grid from z_lens + Z_MIN_OFFSET to Z_MAX
Z_MIN_OFFSET = 0.2
Z_MAX = 10.0

FIELDS = {
    "macs0417": {
        "images": "jauzac19_images.ecsv",
        "models": ("macs0417-canucs-maps", "macs0417-relics-lenstool"),
    },
    "macs1423": {
        "images": "kaleidoscope_m1423_images.ecsv",
        "models": ("macs1423-canucs-maps",),
    },
}


def load_images(path: Path) -> Table:
    t = ms.gold_images(path)
    t["z_spec"] = np.asarray(t["z_spec"], bool)
    return t


def rms_of(d: np.ndarray) -> float:
    ok = np.isfinite(d)
    return float(np.sqrt(np.mean(d[ok] ** 2))) if ok.any() else float("nan")


def gate_k(ip: Table) -> dict:
    d = np.asarray(ip["dtheta_arcsec"], float)
    rms = rms_of(d)
    frac = float(np.mean(np.isfinite(d) & (d <= MATCH_ARCSEC))) if len(d) else 0.0
    return {
        "n_images": len(d),
        "rms_arcsec": rms,
        "matched_fraction": frac,
        "passes": bool(np.isfinite(rms) and rms <= K_MAX_RMS_ARCSEC and frac >= K_MIN_MATCHED),
    }


def best_z(d_by_z: list[np.ndarray]) -> int:
    """Index of the fitted z: smallest rms among the z matching every image, else smallest rms."""
    rms = np.array([rms_of(d) if np.isfinite(d).any() else np.inf for d in d_by_z])
    full = np.array([bool(np.all(np.isfinite(d) & (d <= MATCH_ARCSEC))) for d in d_by_z])
    pool = np.where(full)[0] if full.any() else np.arange(len(d_by_z))
    return int(pool[np.argmin(rms[pool])])


class Solver:
    """One model and its deflection grid (built once), as in ``m1149_sharon.image_plane``."""

    def __init__(self, name: str, images: Table, grid_step: float):
        self.name = name
        self.model, _, _ = lc.load_model(name)
        lc.apply_frame_offset(name, self.model)
        half = lc.grid_half_width(self.model, images)
        self.grid = lensmodel.DeflectionGrid.cached(
            self.model, lc.grid_cache_path(self.model, half, grid_step), half, grid_step
        )

    def image_plane(self, images: Table) -> Table:
        bt = lensmodel.backtrace_images(self.model, images, {})
        ip = lensmodel.imageplane_residuals(self.model, self.grid, bt)
        ip.meta["model"] = self.name
        return ip


def fitted_z(solver: Solver, images: Table) -> Table:
    """Image-plane residuals of the non-spectroscopic families at each family's fitted z."""
    z_lens = lc.MODELS[solver.name]["z_lens"]
    grid_z = np.geomspace(z_lens + Z_MIN_OFFSET, Z_MAX, Z_STEPS)
    families = list(dict.fromkeys(np.asarray(images["system"], str)))
    per_z = []
    for z in grid_z:
        t = images.copy()
        t["z"] = z
        per_z.append(solver.image_plane(t))
    out = []
    for fam in families:
        sel = [np.asarray(ip["system"], str) == fam for ip in per_z]
        k = best_z(
            [np.asarray(ip["dtheta_arcsec"], float)[s] for ip, s in zip(per_z, sel, strict=True)]
        )
        rows = per_z[k][sel[k]]
        rows["z_fit"] = np.full(len(rows), grid_z[k])
        out.append(rows)
    return vstack(out)


def drop_unreproduced(solver: Solver, sp: Table, ip_sp: Table) -> dict:
    """Post hoc (deviation from PLAN.md): gate K and each affected family re-solved without the
    spectroscopic images that are not reproduced (their own back-traces moved the family source)."""
    d = np.asarray(ip_sp["dtheta_arcsec"], float)
    bad = ~(np.isfinite(d) & (d <= MATCH_ARCSEC))
    dropped = [str(i) for i in ip_sp["image_id"][bad]]
    keep = ~np.isin(np.asarray(sp["image_id"], str), dropped)
    ip = solver.image_plane(sp[keep])
    fams = sorted({str(f) for f in ip_sp["system"][bad]})
    m = np.isin(np.asarray(ip["system"], str), fams)
    # leave-one-out (over all its images) on spectroscopic families whose re-solved rms, without
    # the dropped images, still exceeds K_MAX_RMS_ARCSEC
    loo = {}
    for fam in dict.fromkeys(np.asarray(ip["system"], str)):
        f = np.asarray(ip["system"], str) == fam
        if rms_of(np.asarray(ip["dtheta_arcsec"], float)[f]) <= K_MAX_RMS_ARCSEC or f.sum() < 3:
            continue
        sub = sp[np.asarray(sp["system"], str) == fam]  # all its images, dropped ones included
        loo[fam] = {
            str(out): rms_of(
                np.asarray(solver.image_plane(sub[sub["image_id"] != out])["dtheta_arcsec"], float)
            )
            for out in sub["image_id"]
        }
    # the image whose removal helps most: its offset from the nearest image predicted for the
    # source of the remaining ones (and that prediction's magnification)
    held_out = {}
    for fam, by_out in loo.items():
        out_id = min(by_out, key=by_out.get)
        sub = sp[np.asarray(sp["system"], str) == fam]
        rest = sub[sub["image_id"] != out_id]
        bt = lensmodel.backtrace_images(solver.model, rest, {})
        pred = lensmodel.find_images(
            solver.model,
            solver.grid,
            float(np.mean(bt["beta_x"])),
            float(np.mean(bt["beta_y"])),
            float(bt["z_used"][0]),
        )
        row = sub[sub["image_id"] == out_id]
        ox, oy = solver.model.to_frame(np.asarray(row["ra"], float), np.asarray(row["dec"], float))
        dist = np.hypot(np.asarray(pred["x"]) - ox[0], np.asarray(pred["y"]) - oy[0])
        k = int(np.argmin(dist))
        held_out[out_id] = {
            "n_predicted": len(pred),
            "nearest_predicted_arcsec": float(dist[k]),
            "nearest_predicted_magnification": float(pred["magnification"][k]),
        }
    return {
        "dropped": dropped,
        "leave_one_out_rms": loo,
        "held_out_offset": held_out,
        "gate_K_without_dropped": gate_k(ip),
        "families_refit": {str(r["image_id"]): float(r["dtheta_arcsec"]) for r in ip[m]},
    }


def map_checks(model, name: str) -> dict:
    """Post hoc: κ and |μ(z = 2)| from the deflection maps vs the published maps (loading check)."""
    files = lc.model_files(name)
    out = {
        "kappa": lc.map_check(model, files["kappa_map"], "kappa", None, 37),
        "mu_z2": lc.map_check(model, files["mag_z2"], "mu", 2.0, 37),
    }
    for v in out.values():
        v["map"] = Path(v["map"]).name  # cache file name, not a local path
    return out


def readings(tables: dict[str, Table], passing: list[str], images: Table) -> Table:
    ids = np.asarray(images["image_id"], str)
    rows = Table({"image_id": ids, "system": np.asarray(images["system"], str)})
    rows["z_spec"] = np.asarray(images["z_spec"], bool)
    not_rep = np.ones(len(ids), bool)
    for name, ip in tables.items():
        lookup = dict(
            zip(
                np.asarray(ip["image_id"], str), np.asarray(ip["dtheta_arcsec"], float), strict=True
            )
        )
        d = np.array([lookup.get(i, np.nan) for i in ids])
        rows[f"dtheta_{name}"] = d
        rep = np.isfinite(d) & (d <= MATCH_ARCSEC)
        rows[f"reproduced_{name}"] = rep
        if name in passing:
            not_rep &= ~rep
    flag = not_rep if passing else np.zeros(len(ids), bool)
    rows["reading"] = np.where(flag, "lensing_residual_input", "no_flag")
    rows.meta["provenance"] = "model_prediction"
    rows.meta["source"] = f"m0417_images readings; gate-passing models {passing}"
    return rows


def run_field(field: str, out: Path, grid_step: float) -> dict:
    spec = FIELDS[field]
    images = load_images(out / spec["images"])
    sp, nsp = images[images["z_spec"]], images[~images["z_spec"]]
    summary: dict = {"images": len(images), "spectroscopic": len(sp)}
    tables, passing = {}, []
    for name in spec["models"]:
        solver = Solver(name, images, grid_step)
        ip_sp = solver.image_plane(sp)
        k = gate_k(ip_sp)
        entry = {"model_sha256": solver.model.sha256, "gate_K": k}
        ip_nsp = fitted_z(solver, nsp) if len(nsp) else ip_sp[:0]
        if len(nsp):
            entry["z_fit"] = {
                str(s): float(z)
                for s, z in dict(zip(ip_nsp["system"], ip_nsp["z_fit"], strict=True)).items()
            }
        ip = vstack([ip_sp, ip_nsp])
        lc._write(ip, out / f"imageplane_{field}_{name}.ecsv")
        counts = {}
        for fam in dict.fromkeys(np.asarray(ip["system"], str)):
            m = np.asarray(ip["system"], str) == fam
            counts[fam] = (int(m.sum()), int(np.max(ip["n_predicted"][m])))
        entry["extra_predicted"] = {f: c for f, c in counts.items() if c[1] > c[0]}
        entry["fewer_predicted"] = {f: c for f, c in counts.items() if c[1] < c[0]}
        entry["post_hoc"] = drop_unreproduced(solver, sp, ip_sp)
        if "kappa_map" in lc.MODELS[name]["files"] and "mag_z2" in lc.MODELS[name]["files"]:
            entry["post_hoc"]["map_check"] = map_checks(solver.model, name)
        tables[name] = ip
        if k["passes"]:
            passing.append(name)
        summary[name] = entry
    r = readings(tables, passing, images)
    lc._write(r, out / f"readings_{field}.ecsv")
    summary["passing_models"] = passing
    summary["flags"] = [str(i) for i in r["image_id"][r["reading"] == "lensing_residual_input"]]
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--grid-step", type=float, default=0.25)
    ap.add_argument("--fields", nargs="+", default=list(FIELDS))
    args = ap.parse_args(argv)
    summary = {f: run_field(f, args.out, args.grid_step) for f in args.fields}
    summary["rules"] = {
        "K_MAX_RMS_ARCSEC": K_MAX_RMS_ARCSEC,
        "K_MIN_MATCHED": K_MIN_MATCHED,
        "MATCH_ARCSEC": MATCH_ARCSEC,
        "Z_STEPS": Z_STEPS,
        "Z_MIN_OFFSET": Z_MIN_OFFSET,
        "Z_MAX": Z_MAX,
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({f: (summary[f]["passing_models"], summary[f]["flags"]) for f in args.fields}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
