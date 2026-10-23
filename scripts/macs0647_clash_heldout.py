"""MACS0647: the HST-era CLASH Zitrin v2 maps on the systems first identified in JWST data.

Frozen plan: ``results/lens_consistency/macs0647_clash_heldout/PLAN.md``. Images: Meena+2023
(``macs0647_inputs/meena23_images.ecsv``, observed; ``uncertain`` rows left out). Per model: each
family at its fitted z (rule 1), gate K0 and frame check F on the pre-JWST families (rules 2-3),
then the held-out reading of families 13-23 (rule 4). Writes ``summary.json`` and per-model
image-plane tables.

    python scripts/macs0647_clash_heldout.py [--grid-step 0.25]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lens_consistency as lc
import m0417_images as mi
import numpy as np
from astropy.table import Table

from jwst_anomaly import lensmodel

OUT = Path("results/lens_consistency/macs0647_clash_heldout")
IMAGES = Path("results/lens_consistency/macs0647_inputs/meena23_images.ecsv")
MODELS = ("macs0647-clash-nfw", "macs0647-clash-ltm")
# Frozen rules (PLAN.md; ASSUMPTIONs)
K0_MAX_RMS_ARCSEC = 2.0  # rule 2
K0_MIN_MATCHED = 0.8  # rule 2
MATCH_ARCSEC = 3.0  # rules 1, 2, 4 (as m0417_images)
FRAME_MAX_ARCSEC = 0.5  # rule 3, per axis
Z_STEPS = 60  # rule 1
Z_MIN_OFFSET = 0.2
Z_MAX = 12.0  # rule 1 deviation from m0417_images (10)


def load_images(path: Path = IMAGES) -> Table:
    """Certain images of families with >= 2 certain images; ``held_out`` = ``new_in_jwst``."""
    t = Table.read(path, format="ascii.ecsv")
    t = t[~np.asarray(t["uncertain"], bool)]
    fam = np.asarray(t["system"], str)
    keep = np.array([np.sum(fam == f) >= 2 for f in fam])
    t = t[keep]
    t["image_id"] = np.asarray(t["image_id"], str)
    t["system"] = np.asarray(t["system"], str)
    t["held_out"] = np.asarray(t["new_in_jwst"], bool)
    t["a"] = np.ones(len(t))
    t.meta["source"] = f"{t.meta.get('source', path)}; uncertain rows and one-image families out"
    return t


def gate_k0(d: np.ndarray) -> dict:
    """Rule 2: rms <= K0_MAX_RMS_ARCSEC and fraction within MATCH_ARCSEC >= K0_MIN_MATCHED
    (NaN offsets count as unmatched and are left out of the rms)."""
    d = np.asarray(d, float)
    rms = mi.rms_of(d)
    frac = float(np.mean(np.isfinite(d) & (d <= MATCH_ARCSEC))) if len(d) else 0.0
    return {
        "n_images": len(d),
        "rms_arcsec": rms,
        "matched_fraction": frac,
        "passes": bool(np.isfinite(rms) and rms <= K0_MAX_RMS_ARCSEC and frac >= K0_MIN_MATCHED),
    }


def frame_check(dx: np.ndarray, dy: np.ndarray) -> dict:
    """Rule 3: |median dx| <= FRAME_MAX_ARCSEC and |median dy| <= FRAME_MAX_ARCSEC (finite rows)."""
    dx, dy = np.asarray(dx, float), np.asarray(dy, float)
    ok = np.isfinite(dx) & np.isfinite(dy)
    if not ok.any():
        return {"n": 0, "median_dx": float("nan"), "median_dy": float("nan"), "passes": False}
    mx, my = float(np.median(dx[ok])), float(np.median(dy[ok]))
    return {
        "n": int(ok.sum()),
        "median_dx": mx,
        "median_dy": my,
        "passes": bool(abs(mx) <= FRAME_MAX_ARCSEC and abs(my) <= FRAME_MAX_ARCSEC),
    }


def readings(d_by_model: dict[str, np.ndarray], usable: list[str]) -> np.ndarray:
    """Rule 4: ``lensing_residual_input`` where every usable model misses (> MATCH_ARCSEC or NaN);
    no usable model -> no reading."""
    n = len(next(iter(d_by_model.values())))
    if not usable:
        return np.full(n, "no_flag")
    miss = np.ones(n, bool)
    for name in usable:
        d = np.asarray(d_by_model[name], float)
        miss &= ~(np.isfinite(d) & (d <= MATCH_ARCSEC))
    return np.where(miss, "lensing_residual_input", "no_flag")


def signed_offsets(solver: mi.Solver, ip: Table, images: Table) -> tuple[np.ndarray, np.ndarray]:
    """Observed minus nearest predicted image (arcsec, x = West, y = North) at each row's ``z_used``
    from ``ip`` (the fitted z); NaN where the family has no prediction."""
    zmap = dict(zip(np.asarray(ip["image_id"], str), np.asarray(ip["z_used"], float), strict=True))
    t = images.copy()
    t["z"] = [zmap[i] for i in np.asarray(t["image_id"], str)]
    bt = lensmodel.backtrace_images(solver.model, t, {})
    dx = np.full(len(t), np.nan)
    dy = np.full(len(t), np.nan)
    fam = np.asarray(t["system"], str)
    for f in dict.fromkeys(fam):
        sel = np.where((fam == f) & np.isfinite(bt["beta_x"]))[0]
        if len(sel) < 2:
            continue
        pred = lensmodel.find_images(
            solver.model,
            solver.grid,
            float(np.mean(bt["beta_x"][sel])),
            float(np.mean(bt["beta_y"][sel])),
            float(bt["z_used"][sel[0]]),
        )
        if not len(pred):
            continue
        ox, oy = solver.model.to_frame(
            np.asarray(t["ra"], float)[sel], np.asarray(t["dec"], float)[sel]
        )
        px, py = np.asarray(pred["x"], float), np.asarray(pred["y"], float)
        k = np.argmin(np.hypot(ox[:, None] - px[None, :], oy[:, None] - py[None, :]), axis=1)
        dx[sel], dy[sel] = ox - px[k], oy - py[k]
    # model frame x = West; RA cos(dec) offset (East positive) = -dx
    return -dx, dy


def at_zphot(solver: mi.Solver, images: Table) -> dict:
    """Post hoc (not in PLAN.md): every family at the median of its photometric z instead of the
    fitted z (families without one left out), the same statistics for both sets."""
    zp = np.asarray(images["z_phot"], float)
    fam = np.asarray(images["system"], str)
    zf = {
        f: np.nanmedian(zp[fam == f]) for f in dict.fromkeys(fam) if np.isfinite(zp[fam == f]).any()
    }
    t = images[np.isin(fam, list(zf))]
    t["z"] = [zf[f] for f in np.asarray(t["system"], str)]
    ip = solver.image_plane(t)
    d = np.asarray(ip["dtheta_arcsec"], float)
    held = np.asarray(t["held_out"], bool)
    out = {"pre_jwst": gate_k0(d[~held]), "held_out": gate_k0(d[held])}
    for v in out.values():
        v.pop("passes")
    out["held_out_unmatched"] = [
        str(i)
        for i, x in zip(t["image_id"][held], d[held], strict=True)
        if not (np.isfinite(x) and x <= MATCH_ARCSEC)
    ]
    out["families_left_out"] = sorted(set(fam) - set(zf))
    return out


def run_model(name: str, images: Table, grid_step: float) -> tuple[dict, Table]:
    mi.Z_STEPS, mi.Z_MIN_OFFSET, mi.Z_MAX, mi.MATCH_ARCSEC = (
        Z_STEPS,
        Z_MIN_OFFSET,
        Z_MAX,
        MATCH_ARCSEC,
    )
    solver = mi.Solver(name, images, grid_step)
    ip = mi.fitted_z(solver, images)
    order = {i: k for k, i in enumerate(np.asarray(ip["image_id"], str))}
    ip = ip[[order[i] for i in np.asarray(images["image_id"], str)]]
    held = np.asarray(images["held_out"], bool)
    d = np.asarray(ip["dtheta_arcsec"], float)
    dra, ddec = signed_offsets(solver, ip, images)
    ip["dra_cosdec_arcsec"] = dra
    ip["ddec_arcsec"] = ddec
    ip["held_out"] = held
    zp = np.asarray(images["z_phot"], float)
    fam = np.asarray(images["system"], str)
    z_fit = {f: float(ip["z_fit"][fam == f][0]) for f in dict.fromkeys(fam)}
    z_phot = {
        f: (float(np.nanmedian(zp[fam == f])) if np.isfinite(zp[fam == f]).any() else None)
        for f in z_fit
    }
    counts = {f: (int(np.sum(fam == f)), int(np.max(ip["n_predicted"][fam == f]))) for f in z_fit}
    entry = {
        "model_sha256": solver.model.sha256,
        "gate_K0": gate_k0(d[~held]),
        "frame_F": frame_check(dra[~held], ddec[~held]),
        "held_out": gate_k0(d[held]),  # same statistics; "passes" is not a rule for this set
        "z_fit": z_fit,
        "z_phot_median": z_phot,
        "extra_predicted": {f: c for f, c in counts.items() if c[1] > c[0]},
        "fewer_predicted": {f: c for f, c in counts.items() if c[1] < c[0]},
    }
    entry["held_out"].pop("passes")
    entry["post_hoc_zphot"] = at_zphot(solver, images)
    entry["usable"] = entry["gate_K0"]["passes"] and entry["frame_F"]["passes"]
    return entry, ip


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grid-step", type=float, default=0.25)
    args = ap.parse_args(argv)
    images = load_images()
    summary: dict = {
        "plan": str(OUT / "PLAN.md"),
        "images": len(images),
        "held_out_images": int(np.sum(images["held_out"])),
        "rules": {
            "K0_MAX_RMS_ARCSEC": K0_MAX_RMS_ARCSEC,
            "K0_MIN_MATCHED": K0_MIN_MATCHED,
            "MATCH_ARCSEC": MATCH_ARCSEC,
            "FRAME_MAX_ARCSEC": FRAME_MAX_ARCSEC,
            "z_grid": [Z_STEPS, Z_MIN_OFFSET, Z_MAX],
            "grid_step": args.grid_step,
        },
    }
    d_by_model = {}
    for name in MODELS:
        entry, ip = run_model(name, images, args.grid_step)
        summary[name] = entry
        d_by_model[name] = np.asarray(ip["dtheta_arcsec"], float)
        ip.meta["model"] = name
        lc._write(ip, OUT / f"imageplane_{name}.ecsv")
    usable = [n for n in MODELS if summary[n]["usable"]]
    held = np.asarray(images["held_out"], bool)
    r = readings({n: d[held] for n, d in d_by_model.items()}, usable)
    summary["usable_models"] = usable
    summary["flags"] = [
        str(i) for i, v in zip(images["image_id"][held], r, strict=True) if v != "no_flag"
    ]
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: summary[k] for k in ("usable_models", "flags")}, indent=1))
    for n in MODELS:
        e = summary[n]
        print(n, e["gate_K0"], e["frame_F"], e["held_out"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
