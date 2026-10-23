"""Abell 2744: gate K of the UNCOVER-reference extra-image test (stages 1 and 1b).

Frozen plans: ``results/lens_consistency/abell2744_uncoverref/PLAN.md`` (stage 1), ``PLAN_1b.md``
(``--stage 1b``: gate K with the WLS source of ``a2744_uncover``; a deviation after stage 1 failed)
and ``PLAN_2a.md`` (``--stage 2a``: the position-angle convention of the UNCOVER DR3 SUPER
catalogue, measured on the release's segmentation map; needs ``--cat`` and ``--segmap``).
Runs gate K of ``m0417_images`` on the spectroscopic images of the UNCOVER v2.0 multiple-image
table with the v2.0 deflection maps (``abell2744-uncover2-maps``). The extra-image stages follow
in a later cycle only if the gate passes. Writes ``summary.json`` to the plan's directory.

    python scripts/a2744_uncoverref.py [--stage 1b]
    python scripts/a2744_uncoverref.py --stage 2a --cat <SUPER.fits> --segmap <SEGMAP.fits.gz>
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import a2744_uncover as au
import lens_consistency as lc
import m0417_canucsref as cr
import m0417_images as mi
import numpy as np
from astropy.io import fits
from astropy.table import Table
from m0417_canucsref_posthoc import _clean

from jwst_anomaly import lensmodel
from jwst_anomaly.photometry import fetch_catalog
from jwst_anomaly.provenance import file_sha256

FIELD = "abell2744"
MODEL = "abell2744-uncover2-maps"
OUT = Path("results/lens_consistency/abell2744_uncoverref")
IMAGES = "uncover2_images.ecsv"
GRID_STEP_K = 0.25  # m0417_images default (PLAN.md rule 1)
# stage 2a inputs (PLAN_2a.md), pinned from the file bytes at first download
SUPER_CAT = "UNCOVER_DR3_LW_SUPER_catalog.fits"
SUPER_SHA = "b8f265311d041b815e1511d74eadec2d88bf8c80ef2fee2394d7f78dbe4adb4c"
SEGMAP = "UNCOVER_DR3_SEGMAP.fits.gz"
SEGMAP_SHA = "5a64d01fcd75ecbea38efcc4ec673d2cb5fb3fb19b76ced21417fed6bc00087a"
# PLAN_2a.md rules 1 and 3 (ASSUMPTIONs)
MIN_A_PIX, MIN_ELL, MIN_SEG_PIX, N_MAX = 4.0, 0.4, 50, 500
MEDIAN_MAX_DEG, WITHIN_DEG, WITHIN_FRAC = 10.0, 20.0, 0.70

lc.MODELS[MODEL] = {
    "files": {k: au.UNCOVER_V2[k] for k in ("alpha_x", "alpha_y")},
    "kind": "maps",
    "z_lens": au.Z_LENS,
    "cosmology": au.COSMOLOGY,
}


def lenstool_id(image_id: str) -> str:
    """UNCOVER id -> an id whose ``lensmodel.image_family`` is ``au.uncover_system(image_id)``:
    "S.KI" -> "S.K" + letter I, "S.I" -> "S" + letter I (PLAN.md Inputs)."""
    idx = int(image_id.split(".")[1][-1])
    if idx == 0:  # "S.K0" / "S.10" cannot be told apart once stored as floats
        raise ValueError(f"UNCOVER id {image_id!r}: image number 0 is ambiguous")
    return f"{au.uncover_system(image_id)}{chr(ord('a') + idx - 1)}"


def image_list(raw: Table) -> Table:
    """The release's image table with ``z_spec`` from ``Redshift_type``; ``None`` rows dropped."""
    kind = np.asarray(raw["Redshift_type"], str)
    keep = kind != "None"
    ids = [repr(float(i)) for i in raw["ID"][keep]]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate UNCOVER image ids")
    t = Table(
        {
            "image_id": [lenstool_id(i) for i in ids],
            "uncover_id": ids,
            "ra": np.asarray(raw["RA"], float)[keep],
            "dec": np.asarray(raw["Dec"], float)[keep],
            "z": np.asarray(raw["z"], float)[keep],
            "z_spec": kind[keep] == "spec",
        }
    )
    t["ra"].unit = t["dec"].unit = "deg"
    t.meta["provenance"] = "observed"
    t.meta["source"] = (
        f"{au.UNCOVER_V2['images'][0]} (UNCOVER v2.0 multiple images; z_spec = Redshift_type spec)"
    )
    return t


def read_summary() -> dict:
    p = OUT / "summary.json"
    return json.loads(p.read_text()) if p.exists() else {}


def write_summary(summary: dict) -> None:
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")


def gate_k_wls(model, grid, sp: Table) -> tuple[Table, dict]:
    """PLAN_1b.md rule 1: ``a2744_uncover.imageplane`` (WLS family source) on the spectroscopic
    images ``sp``, decided by the unchanged ``m0417_images.gate_k``; rule 2's reported counts."""
    bt = lensmodel.backtrace_images(model, sp, {})
    ip = au.imageplane(model, grid, bt)
    d = np.asarray(ip["dtheta_arcsec"], float)
    bad = ~(np.isfinite(d) & (d <= mi.MATCH_ARCSEC))
    counts = {}
    for fam in dict.fromkeys(np.asarray(ip["system"], str)):
        m = np.asarray(ip["system"], str) == fam
        counts[fam] = (int(m.sum()), int(np.max(ip["n_predicted"][m])))
    return ip, {
        "gate_K": mi.gate_k(ip),
        "estimator": "wls_source (a2744_uncover.imageplane)",
        "not_reproduced": {
            str(i): float(x) for i, x in zip(ip["image_id"][bad], d[bad], strict=True)
        },
        "shared_match": int(np.sum(ip["shared_match"])),
        "extra_predicted": {f: c for f, c in counts.items() if c[1] > c[0]},
        "fewer_predicted": {f: c for f, c in counts.items() if c[1] < c[0]},
    }


def stage_1b() -> dict:
    """PLAN_1b.md: gate K-WLS on the committed image list, same Solver grid as stage 1."""
    images = mi.load_images(OUT / IMAGES)
    sp = images[images["z_spec"]]
    solver = mi.Solver(MODEL, images, GRID_STEP_K)
    ip, k = gate_k_wls(solver.model, solver.grid, sp)
    ip.meta.update(solver.model._meta())
    ip.meta["provenance"] = "model_prediction"
    ip.meta["source"] = f"{MODEL}, WLS family source (PLAN_1b.md)"
    lc._write(ip, OUT / f"imageplane_wls_{FIELD}_{MODEL}.ecsv")
    k["model_sha256"] = solver.model.sha256
    if k["gate_K"]["passes"]:
        k["next"] = "gate K-WLS passed (deviation): freeze stage 2 with the WLS source"
    else:
        k["stopped"] = "gate K-WLS failed: model not used (PLAN_1b.md rule 1)"
    return k


def y_axis_pa_deg(x, y, ra, dec) -> float:
    """Sky position angle (east of north) of the pixel grid's +y axis, from a linear fit of the
    catalogue's pixel to sky positions (0 for a north-up grid)."""
    x, y, ra, dec = (np.asarray(v, float) for v in (x, y, ra, dec))
    d0 = np.median(dec)
    a = np.c_[x - np.median(x), y - np.median(y), np.ones_like(x)]
    e = np.linalg.lstsq(a, (ra - np.median(ra)) * np.cos(np.radians(d0)), rcond=None)[0]
    n = np.linalg.lstsq(a, dec - d0, rcond=None)[0]
    return float(np.degrees(np.arctan2(e[1], n[1])))


def theta_sample(t: Table) -> np.ndarray:
    """PLAN_2a.md rule 1 before the segment-size cut: row indices, largest ``a_image`` first
    (ties by ``id``)."""
    a = np.asarray(t["a_image"], float)
    b = np.asarray(t["b_image"], float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ell = 1.0 - b / a
    keep = (
        (np.asarray(t["use_phot"]) == 1)
        & (np.asarray(t["flag_kron"]) == 0)
        & (np.asarray(t["flag_nearbcg"]) == 0)
        & (a >= MIN_A_PIX)
        & (ell >= MIN_ELL)
    )
    idx = np.flatnonzero(keep)
    return idx[np.lexsort((np.asarray(t["id"])[idx], -a[idx]))]


def segment_angle(ys, xs) -> float:
    """PLAN_2a.md rule 2: orientation of a segment's pixels from unweighted second moments
    (radians, counter-clockwise from +x)."""
    x = np.asarray(xs, float) - np.mean(xs)
    y = np.asarray(ys, float) - np.mean(ys)
    return float(0.5 * np.arctan2(2.0 * np.mean(x * y), np.mean(x * x) - np.mean(y * y)))


def fold_deg(d) -> np.ndarray:
    """Axis-angle difference (radians) folded into [0, 90] degrees."""
    r = np.mod(np.degrees(np.asarray(d, float)), 180.0)
    return np.minimum(r, 180.0 - r)


def theta_check(theta, phi_seg) -> dict:
    """PLAN_2a.md rules 2-3 on catalogue ``theta`` and segment angles ``phi_seg`` (radians)."""
    theta = np.asarray(theta, float)
    phi = np.asarray(phi_seg, float)
    out: dict = {"n": int(len(theta))}
    for name, pred in (("H_a_ccw_from_x", theta), ("H_b_cw_from_y", np.pi / 2 - theta)):
        r = fold_deg(phi - pred)
        med = float(np.median(r)) if len(r) else float("nan")
        frac = float(np.mean(r <= WITHIN_DEG)) if len(r) else 0.0
        out[name] = {
            "median_deg": med,
            "frac_within_20deg": frac,
            "passes": bool(med <= MEDIAN_MAX_DEG and frac >= WITHIN_FRAC),
        }
    ok = [k for k in ("H_a_ccw_from_x", "H_b_cw_from_y") if out[k]["passes"]]
    out["adopted"] = min(ok, key=lambda k: out[k]["median_deg"]) if ok else None
    return out


def stage_2a(cat: Path, segmap: Path) -> dict:
    """PLAN_2a.md: the ``theta_J2000`` convention of the DR3 SUPER catalogue on its segmap."""
    if file_sha256(cat) != SUPER_SHA or file_sha256(segmap) != SEGMAP_SHA:
        raise SystemExit("error: --cat / --segmap are not the pinned UNCOVER DR3 files")
    from scipy import ndimage

    t = Table.read(cat)
    with fits.open(segmap) as hl:
        seg = np.asarray(next(h for h in hl if h.data is not None).data)
    order = theta_sample(t)
    ids = np.asarray(t["id"], int)
    slices = ndimage.find_objects(seg)
    used, phis, thetas = [], [], []
    for i in order:
        sid = int(ids[i])
        sl = slices[sid - 1] if 0 < sid <= len(slices) else None
        if sl is None:
            continue
        ys, xs = np.nonzero(seg[sl] == sid)
        if len(ys) < MIN_SEG_PIX:
            continue
        used.append(sid)
        phis.append(segment_angle(ys, xs))
        thetas.append(float(t["theta_J2000"][i]))
        if len(used) == N_MAX:
            break
    res = theta_check(thetas, phis)
    res["inputs"] = {
        "catalogue": SUPER_CAT,
        "catalogue_sha256": SUPER_SHA,
        "segmap": SEGMAP,
        "segmap_sha256": SEGMAP_SHA,
        "n_rule1_before_segment_cut": int(len(order)),
        "theta_range_rad": [float(np.min(t["theta_J2000"])), float(np.max(t["theta_J2000"]))],
        "pixel_scale_arcsec": float(cr.pixel_scale(t["x"], t["y"], t["ra"], t["dec"])),
        "y_axis_pa_deg": y_axis_pa_deg(t["x"], t["y"], t["ra"], t["dec"]),
    }
    res["ids"] = used
    return res


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--stage", choices=("1", "1b", "2a"), default="1")
    ap.add_argument("--cat", type=Path, help="stage 2a: the DR3 SUPER catalogue")
    ap.add_argument("--segmap", type=Path, help="stage 2a: the DR3 segmentation map")
    args = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage == "2a":
        if args.cat is None or args.segmap is None:
            ap.error("--stage 2a needs --cat and --segmap")
        summary = read_summary()
        summary["stage2a_theta_convention"] = _clean(stage_2a(args.cat, args.segmap))
        write_summary(summary)
        r = summary["stage2a_theta_convention"]
        print(json.dumps({k: r[k] for k in ("n", "H_a_ccw_from_x", "H_b_cw_from_y", "adopted")}))
        return 0
    if args.stage == "1b":
        summary = read_summary()
        summary["stage1b_gate_K_wls"] = _clean(stage_1b())
        summary["status"] = (
            "stage 1 stopped (mean source); stage 1b passed (WLS source, deviation): stage 2 next"
            if summary["stage1b_gate_K_wls"]["gate_K"]["passes"]
            else "stage 1 and stage 1b failed: Abell 2744 extra-image route stopped with this model"
        )
        write_summary(summary)
        print(json.dumps(summary["stage1b_gate_K_wls"]["gate_K"]))
        return 0
    raw = Table.read(fetch_catalog(*au.UNCOVER_V2["images"]))
    image_list(raw).write(OUT / IMAGES, format="ascii.ecsv", overwrite=True)
    mi.FIELDS[FIELD] = {"images": IMAGES, "models": (MODEL,)}
    k = mi.run_field(FIELD, OUT, GRID_STEP_K)
    # a stage-1 rerun keeps the later stages' keys (stage 1b adds to the same file)
    summary: dict = {
        x: v for x, v in read_summary().items() if x.startswith("stage1b") or x == "status"
    }
    summary["stage1_gate_K"] = k
    if MODEL not in k["passing_models"]:
        summary["stopped"] = "gate K failed: model not used (PLAN.md rule 1)"
    else:
        summary["next"] = "gate K passed: pin an UNCOVER reference, run the extra-image stages"
    write_summary(_clean(summary))
    print(json.dumps({"gate_K": k[MODEL]["gate_K"], **{x: summary.get(x) for x in ("stopped",)}}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
