"""MARS HFF maps vs JWST-era CANUCS images absent from the MARS catalogue (held-out stage).

Frozen plan: ``results/lens_consistency/mars_hff_heldout/PLAN.md``. Per field (Abell 370,
MACS0416): frame offset O from CANUCS to MARS positions, seen / held-out split, per-family anchors
(type E / N / M), exact image-plane residuals of the held-out images, gate H. No flag, flux or
extra-image reading.

    python scripts/mars_hff_heldout.py --zip DIR

``DIR`` receives the MARS zip (md5 checked); it is deleted after the run unless ``--keep``.
"""

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

import lens_consistency as lc
import m0417_images as mi
import mars_hff_inputs as mh
import numpy as np
import requests
from astropy.table import Table
from m0417_canucsref_posthoc import _clean

from jwst_anomaly import lensmodel

OUT = Path("results/lens_consistency/mars_hff_heldout")
LC = Path("results/lens_consistency")
# MARS catalogue field -> (CANUCS image list, training gate K rms from D-034 addendum ao)
FIELDS = {
    "Abell 370": (LC / "abell370_canucsref/canucs_a370_images.ecsv", 0.850),
    "MACS J0416.1-2403": (LC / "macs0416_canucsref/canucs_m0416_images.ecsv", 0.513),
}
OFFSET_PAIR_ARCSEC = 1.0  # ASSUMPTION (PLAN.md rule 1)
OFFSET_MIN_PAIRS = 5
SEEN_ARCSEC = 0.5  # ASSUMPTION (PLAN.md rule 2)
MIN_TESTED = 5  # ASSUMPTION (PLAN.md rule 5)


def _offsets_arcsec(ra, dec, ra_ref, dec_ref) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Tangent-plane offsets (arcsec) from every (ra, dec) to every reference: dx, dy, distance."""
    ra, dec, ra_ref, dec_ref = (np.asarray(v, float) for v in (ra, dec, ra_ref, dec_ref))
    cosd = np.cos(np.deg2rad(dec))[:, None]
    dx = (ra_ref[None, :] - ra[:, None]) * cosd * 3600.0
    dy = (dec_ref[None, :] - dec[:, None]) * 3600.0
    return dx, dy, np.hypot(dx, dy)


def frame_offset(canucs: Table, mars: Table, pair_arcsec: float = OFFSET_PAIR_ARCSEC) -> dict:
    """PLAN.md rule 1: median offset (arcsec, east, north) from CANUCS to the nearest MARS row."""
    if not len(mars) or not len(canucs):
        return {"n_pairs": 0, "dx_arcsec": 0.0, "dy_arcsec": 0.0, "applied": False}
    dx, dy, d = _offsets_arcsec(canucs["ra"], canucs["dec"], mars["ra"], mars["dec"])
    k = np.argmin(d, axis=1)
    i = np.arange(len(canucs))
    ok = d[i, k] <= pair_arcsec
    n = int(ok.sum())
    if n < OFFSET_MIN_PAIRS:
        return {"n_pairs": n, "dx_arcsec": 0.0, "dy_arcsec": 0.0, "applied": False}
    return {
        "n_pairs": n,
        "dx_arcsec": float(np.median(dx[i, k][ok])),
        "dy_arcsec": float(np.median(dy[i, k][ok])),
        "applied": True,
    }


def shift(canucs: Table, off: dict) -> Table:
    """CANUCS positions moved by the frame offset into the MARS frame."""
    t = canucs.copy()
    dec = np.asarray(t["dec"], float)
    t["ra"] = np.asarray(t["ra"], float) + off["dx_arcsec"] / 3600.0 / np.cos(np.deg2rad(dec))
    t["dec"] = dec + off["dy_arcsec"] / 3600.0
    return t


def split_seen(canucs: Table, mars: Table, seen_arcsec: float = SEEN_ARCSEC) -> Table:
    """PLAN.md rule 2: images with listed z > 0; seen when a MARS row is within ``seen_arcsec``."""
    t = canucs[np.asarray(canucs["z"], float) > 0]
    _, _, d = _offsets_arcsec(t["ra"], t["dec"], mars["ra"], mars["dec"])
    t["sep_mars_arcsec"] = d.min(axis=1) if len(mars) else np.full(len(t), np.inf)
    t["seen"] = t["sep_mars_arcsec"] <= seen_arcsec
    t["system"] = [lensmodel.image_family(i) for i in t["image_id"]]
    return t


def anchors(t: Table) -> dict[str, tuple[str, np.ndarray]]:
    """PLAN.md rule 3: per family with a held-out image, (type, mask of the anchor rows)."""
    sysv = np.asarray(t["system"], str)
    seen = np.asarray(t["seen"], bool)
    out = {}
    for s in dict.fromkeys(sysv[~seen]):
        fam = sysv == s
        ns = int((fam & seen).sum())
        if ns >= 2:
            out[s] = ("E", fam & seen)
        elif fam.sum() >= 2:
            out[s] = ("N" if ns == 0 else "M", fam)
    return out


def heldout_residuals(model, grid, t: Table) -> Table:
    """PLAN.md rule 3: image-plane residual of every held-out image in a family with an anchor."""
    bt = lensmodel.backtrace_images(model, t, {})
    ox, oy = model.to_frame(np.asarray(t["ra"], float), np.asarray(t["dec"], float))
    rows = []
    for s, (kind, m) in anchors(t).items():
        idx_h = np.where((np.asarray(t["system"], str) == s) & ~np.asarray(t["seen"], bool))[0]
        bx, by = np.asarray(bt["beta_x"])[m], np.asarray(bt["beta_y"])[m]
        if not (np.isfinite(bx).all() and np.isfinite(by).all()):
            for i in idx_h:
                rows.append((t["image_id"][i], s, kind, float(t["z"][i]), np.nan, 0, False))
            continue  # an anchor image outside the map: dropped and counted
        pred = lensmodel.find_images(
            model, grid, float(bx.mean()), float(by.mean()), float(t["z"][idx_h[0]])
        )
        for i in idx_h:
            d = np.hypot(ox[i] - pred["x"], oy[i] - pred["y"]) if len(pred) else np.array([])
            inside = bool(np.isfinite(bt["beta_x"][i]))
            dmin = float(d.min()) if len(d) and inside else (np.inf if inside else np.nan)
            rows.append((t["image_id"][i], s, kind, float(t["z"][i]), dmin, len(pred), inside))
    names = ("image_id", "system", "type", "z", "dtheta_arcsec", "n_predicted", "inside")
    if not rows:
        return Table(names=names, dtype=(str, str, str, float, float, int, bool))
    return Table(rows=rows, names=names)


def gate_h(r: Table) -> dict:
    """PLAN.md rule 4: rms <= 1.0 arcsec and >= 80 % within 3.0 arcsec, on the tested held-out
    images inside the map; no prediction (inf) counts as not reproduced and is out of the rms."""
    r = r[np.asarray(r["inside"], bool)] if len(r) else r
    d = np.asarray(r["dtheta_arcsec"], float) if len(r) else np.array([])
    rms = mi.rms_of(d)
    frac = float(np.mean(np.isfinite(d) & (d <= mi.MATCH_ARCSEC))) if len(d) else 0.0
    return {
        "n_tested": len(d),
        "rms_arcsec": rms,
        "matched_fraction": frac,
        "passes": bool(
            np.isfinite(rms) and rms <= mi.K_MAX_RMS_ARCSEC and frac >= mi.K_MIN_MATCHED
        ),
    }


def outcome(g: dict, min_tested: int = MIN_TESTED) -> str:
    """PLAN.md rule 5."""
    if g["n_tested"] < min_tested:
        return "not testable"
    return "reproduces" if g["passes"] else "does not reproduce"


def field_block(zf, field: str, cat: Table, out: Path) -> dict:
    folder, tag, z_lens = mh.FIELDS[field]
    path, k_train = FIELDS[field]
    base = f"{mh.ROOT}{folder}/deflection_angle_{tag}_"
    (ax, hx), (ay, _) = (mh._fits(zf, base + f"{c}_w_header.fits") for c in "xy")
    model = mh.build_model(
        ax, ay, hx, False, z_lens, "see mars_hff_inputs", f"mars-{tag} (A_arcsec)"
    )
    mars = cat[np.asarray(cat["field"]) == field]
    canucs = Table.read(path)
    off = frame_offset(canucs, mars)
    t = split_seen(shift(canucs, off), mars)
    grid = lensmodel.DeflectionGrid.compute(model, lc.grid_half_width(model, t), mh.GRID_STEP)
    r = heldout_residuals(model, grid, t)
    r["dtheta_arcsec"].info.format = ".3f"
    r.write(out / f"heldout_{tag}.ecsv", format="ascii.ecsv", overwrite=True)
    g = gate_h(r)
    ge = gate_h(r[np.asarray(r["type"]) == "E"]) if len(r) else gate_h(r)
    held = ~np.asarray(t["seen"], bool)
    return {
        "canucs_images": len(canucs),
        "frame_offset": off,
        "seen": int((~held).sum()),
        "held_out": int(held.sum()),
        "held_out_families": len(set(np.asarray(t["system"], str)[held])),
        "types": {k: int(np.sum(np.asarray(r["type"]) == k)) for k in "ENM"} if len(r) else {},
        "outside_map": int(np.sum(~np.asarray(r["inside"], bool))) if len(r) else 0,
        "gate_H": g,
        "reported": {"type_E_only": ge, "training_gate_K_rms_arcsec": k_train},
        "outcome": outcome(g),
        "over_3arcsec": [str(i) for i in r["image_id"][np.asarray(r["dtheta_arcsec"], float) > 3.0]]
        if len(r)
        else [],
    }


def run(zdir: Path, out: Path = OUT) -> dict:
    meta = mh.fetch_zip(requests.Session(), zdir)
    out.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zdir / mh.ZIP) as zf:
        cat = mh.parse_catalogue(zf.read(mh.CATALOGUE).decode())
        fields = {f: field_block(zf, f, cat, out) for f in FIELDS}
    summary = {"plan": str(out / "PLAN.md"), "zip": meta, "fields": fields}
    (out / "summary.json").write_text(json.dumps(_clean(summary), indent=1, allow_nan=False) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--zip", type=Path, required=True, help="directory for the zip")
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--keep", action="store_true", help="keep the downloaded zip")
    a = ap.parse_args(argv)
    try:
        s = run(a.zip, a.out)
    finally:
        if not a.keep:
            (a.zip / mh.ZIP).unlink(missing_ok=True)
    brief = {
        f: {k: d[k] for k in ("frame_offset", "seen", "held_out", "gate_H", "outcome")}
        for f, d in s["fields"].items()
    }
    print(json.dumps(_clean(brief), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
