"""MACS0417: extra-image test with an arc-aware reference flux (isophotal segment totals of the
JWST pipeline catalogue) in place of the D-034 compact 0.2" aperture.

Frozen plan: ``results/lens_consistency/macs0417_arcref/PLAN.md``. Reads the committed predicted and
control tables of ``../macs0417_frame`` (the lens model is not re-run) and the shifted Jauzac+19 list.

    python scripts/m0417_arcref.py [--cache data/cache/m0417x]
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import lens_consistency as lc
import m0417_extra as mx
import m0417_frame as mf
import numpy as np
from astropy.table import Table

from jwst_anomaly import lensmodel

FRAME = Path("results/lens_consistency/macs0417_frame")
OUT = Path("results/lens_consistency/macs0417_arcref")
# Frozen rules (PLAN.md; ASSUMPTIONs)
ELLIPSE_NSIGMA = 3.0  # counterpart: position inside the source's 3-sigma ellipse
MIN_REF_SNR = 5.0
MAX_REF_MU = 50.0  # D-034 cap, per member image
MAX_REF_SPREAD = 3.0
FLUX_RANGE = (1 / 3, 3.0)
TEN_SIGMA_DMAG = 2.5 * np.log10(2.0)  # depth is ~5 sigma; 10 sigma is 0.753 mag brighter
BASE_SEARCH = 1.0
C2_MIN_TESTABLE, C2_MIN_FRAC, C2_MIN_ARC = 8, 0.8, 3
ARC_MIN_ELL = 0.5
TESTABLE = ("recovered", "confused", "absent", "ambiguous")


def pixel_scale(cat: Table) -> float:
    """Arcsec per pixel from the catalogue's own pixel and sky centroids (median over pairs)."""
    sky = cat["sky_centroid"]
    x, y = np.asarray(cat["xcentroid"], float), np.asarray(cat["ycentroid"], float)
    i, j = np.arange(0, len(cat) - 1, 2), np.arange(1, len(cat), 2)
    dpix = np.hypot(x[i] - x[j], y[i] - y[j])
    ok = dpix > 50
    return float(np.median(sky[i][ok].separation(sky[j][ok]).arcsec / dpix[ok]))


def ellipse_radius(dx, dy, a: float, b: float, pa_deg: float) -> float:
    """Normalised elliptical radius of offset (dx East, dy North; arcsec) for an ellipse with
    semi-axes ``a``, ``b`` (arcsec) and major-axis position angle ``pa_deg`` (East of North)."""
    t = np.radians(pa_deg)
    u = dx * np.sin(t) + dy * np.cos(t)  # along the major axis
    v = dx * np.cos(t) - dy * np.sin(t)
    return float(np.hypot(u / a, v / b))


def counterpart(ra: float, dec: float, cat: dict) -> int:
    """Index of the catalogue source whose 3-sigma ellipse contains (ra, dec), -1 if none (rule 1)."""
    dx = (cat["ra"] - ra) * np.cos(np.radians(dec)) * 3600.0  # source minus position, East
    dy = (cat["dec"] - dec) * 3600.0
    near = np.flatnonzero(np.hypot(dx, dy) < 10.0)
    best, best_r = -1, np.inf
    for k in near:
        r = ellipse_radius(
            -dx[k], -dy[k], ELLIPSE_NSIGMA * cat["a"][k], ELLIPSE_NSIGMA * cat["b"][k], cat["pa"][k]
        )
        if r <= 1.0 and r < best_r:
            best, best_r = int(k), r
    return best


def units(images: list[dict]) -> list[dict]:
    """Merge one family's catalogued images that share a counterpart (rule 2)."""
    out: dict[int, dict] = {}
    solo = []
    for im in images:
        if im["cp"] < 0:
            solo.append({"ids": [im["id"]], "cp": -1, "mus": [abs(im["mu"])]})
            continue
        u = out.setdefault(im["cp"], {"ids": [], "cp": im["cp"], "mus": []})
        u["ids"].append(im["id"])
        u["mus"].append(abs(im["mu"]))
    return list(out.values()) + solo


def reference(cands: list[dict], cat: dict) -> tuple[float, str]:
    """(L_ref, state) from candidate units (rule 3): state '' or 'no_reference' /
    'inconsistent_reference'."""
    usable = []
    for u in cands:
        k = u["cp"]
        if k < 0 or not (cat["snr"][k] > MIN_REF_SNR) or max(u["mus"]) > MAX_REF_MU:
            continue
        usable.append((cat["flux"][k] / sum(u["mus"]), sum(u["mus"])))
    if not usable:
        return float("nan"), "no_reference"
    lum = [v[0] for v in usable]
    if max(lum) > MAX_REF_SPREAD * min(lum):
        return float("nan"), "inconsistent_reference"
    return min(usable, key=lambda v: v[1])[0], ""


def classify(f_pred: float, depth: float, ratios) -> str:
    """Frozen per-target class (rule 6). ``ratios`` = F / F_pred of the sources in the radius."""
    m_pred = -2.5 * np.log10(f_pred / 3631.0)
    if not np.isfinite(m_pred) or m_pred > depth:
        return "undetectable"
    r = np.asarray(ratios, float)
    if ((r >= FLUX_RANGE[0]) & (r <= FLUX_RANGE[1])).any():
        return "recovered"
    if (r > FLUX_RANGE[1]).any():
        return "confused"
    if m_pred <= depth - TEN_SIGMA_DMAG:
        return "absent"
    return "ambiguous"


def gate_c2(classes, arc_ref) -> dict:
    """Frozen gate C2 (rule 7) on control classes and whether each control's reference is an arc."""
    c = np.asarray(classes, str)
    arc = np.asarray(arc_ref, bool)

    def part(m):
        t = c[m & np.isin(c, TESTABLE)]
        ok = int(np.isin(t, ("recovered", "confused")).sum())
        return {"n_testable": int(len(t)), "n_ok": ok, "fraction": ok / len(t) if len(t) else None}

    allp, arcp = part(np.ones(len(c), bool)), part(arc)
    passes = (
        allp["n_testable"] >= C2_MIN_TESTABLE
        and allp["fraction"] >= C2_MIN_FRAC
        and arcp["n_testable"] >= C2_MIN_ARC
        and arcp["fraction"] >= C2_MIN_FRAC
    )
    return {"all": allp, "arc_referenced": arcp, "passes": bool(passes)}


def load_cat(path: Path) -> dict:
    t = Table.read(path)
    scale = pixel_scale(t)
    sky = t["sky_centroid"]
    flux = np.asarray(t["isophotal_flux"], float)
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = flux / np.asarray(t["isophotal_flux_err"], float)
    return {
        "label": np.asarray(t["label"]),
        "ra": sky.ra.deg,
        "dec": sky.dec.deg,
        "a": np.asarray(t["semimajor_sigma"], float) * scale,
        "b": np.asarray(t["semiminor_sigma"], float) * scale,
        "pa": np.asarray(t["sky_orientation"], float),
        "flux": flux,
        "snr": snr,
        "ext": np.asarray(t["is_extended"], bool),
        "ell": np.asarray(t["ellipticity"], float),
        "scale": scale,
    }


def in_radius(ra: float, dec: float, radius: float, cat: dict) -> np.ndarray:
    dx = (cat["ra"] - ra) * np.cos(np.radians(dec)) * 3600.0
    dy = (cat["dec"] - dec) * 3600.0
    return np.flatnonzero(np.hypot(dx, dy) <= radius)


def is_arc(k: int, cat: dict) -> bool:
    return bool(k >= 0 and cat["ext"][k] and cat["ell"][k] >= ARC_MIN_ELL)


def run(cache: Path) -> dict:
    _, _, obs, tile = mx.FIELDS["macs0417"]
    cat_path = cache / Path(mx.S3.format(obs, tile, "cat.ecsv")).name
    cat = load_cat(cat_path)
    depth = lc.depth_mag(lc.load_shapes(cat_path))
    pred = Table.read(FRAME / "predicted_macs0417.ecsv")
    ctrl = Table.read(FRAME / "control_macs0417.ecsv")
    mu = dict(
        zip(
            np.asarray(ctrl["matched_image"], str),
            np.asarray(ctrl["magnification"], float),
            strict=True,
        )
    )
    unpredicted = json.loads((FRAME / "summary.json").read_text())["macs0417"]["unpredicted"]
    dra, ddec = (
        json.loads((FRAME / "summary.json").read_text())["measurement"]["offset"][k]
        for k in ("dra", "ddec")
    )
    lst = mf.shifted_loader(dra, ddec)("macs0417")
    fams: dict[str, list[dict]] = {}
    for r in lst:
        iid = str(r["image_id"])
        if iid not in mu:
            continue  # no matched prediction (unpredicted or outside), no model |mu|
        fams.setdefault(str(r["system"]), []).append(
            {
                "id": iid,
                "ra": float(r["ra"]),
                "dec": float(r["dec"]),
                "mu": mu[iid],
                "cp": counterpart(float(r["ra"]), float(r["dec"]), cat),
            }
        )
    fam_units = {s: units(v) for s, v in fams.items()}

    def radius(sys_id: str) -> float:
        has_unpred = any(lensmodel.image_family(u) == sys_id for u in unpredicted)
        return lc.system_search_radius(pred, sys_id, BASE_SEARCH, has_unpred)

    def test(ra, dec, mu_t, sys_id, cands, exclude) -> dict:
        lref, state = reference(cands, cat)
        rad = radius(sys_id)
        row = {"L_ref": lref, "search_arcsec": rad, "f_pred": float("nan"), "ratios": []}
        if state:
            return {**row, "class": state}
        f_pred = lref * mu_t
        idx = [k for k in in_radius(ra, dec, rad, cat) if k not in exclude]
        ratios = [float(cat["flux"][k] / f_pred) for k in idx]
        return {
            **row,
            "f_pred": f_pred,
            "m_pred": float(-2.5 * np.log10(f_pred / 3631.0)),
            "ratios": ratios,
            "class": classify(f_pred, depth, ratios),
        }

    def ref_cp(cands):
        usable = [
            u
            for u in cands
            if u["cp"] >= 0 and cat["snr"][u["cp"]] > MIN_REF_SNR and max(u["mus"]) <= MAX_REF_MU
        ]
        return min(usable, key=lambda u: sum(u["mus"]))["cp"] if usable else -1

    controls = []
    for s, us in fam_units.items():
        for u in us:
            first = next(im for im in fams[s] if im["id"] == u["ids"][0])
            cands = [v for v in us if v is not u]
            rcp = ref_cp(cands)
            res = test(first["ra"], first["dec"], sum(u["mus"]), s, cands, {rcp})
            controls.append(
                {
                    "system": s,
                    "images": "+".join(u["ids"]),
                    "sum_abs_mu": sum(u["mus"]),
                    "counterpart": int(cat["label"][u["cp"]]) if u["cp"] >= 0 else -1,
                    "ref_counterpart": int(cat["label"][rcp]) if rcp >= 0 else -1,
                    "arc_ref": is_arc(rcp, cat),
                    **{k: v for k, v in res.items() if k != "ratios"},
                    "n_sources": len(res["ratios"]),
                    "best_ratio": (
                        min(res["ratios"], key=lambda q: abs(np.log(q)) if q > 0 else np.inf)
                        if res["ratios"]
                        else float("nan")
                    ),
                }
            )
    c2 = gate_c2([c["class"] for c in controls], [c["arc_ref"] for c in controls])

    targets = []
    cls = np.asarray(pred["image_class"], str)
    for r in pred[~np.isin(cls, mx.NOT_TARGETS)]:
        s = str(r["system"])
        us = fam_units.get(s, [])
        rcp = ref_cp(us)
        res = test(float(r["ra"]), float(r["dec"]), abs(float(r["magnification"])), s, us, {rcp})
        targets.append(
            {
                "system": s,
                "image_class": str(r["image_class"]),
                "ra": float(r["ra"]),
                "dec": float(r["dec"]),
                "abs_mu": abs(float(r["magnification"])),
                "ref_counterpart": int(cat["label"][rcp]) if rcp >= 0 else -1,
                "arc_ref": is_arc(rcp, cat),
                **{k: v for k, v in res.items() if k != "ratios"},
                "n_sources": len(res["ratios"]),
                "max_ratio": max(res["ratios"]) if res["ratios"] else float("nan"),
                "reading": (
                    "lensing_residual_input"
                    if res["class"] == "absent" and c2["passes"]
                    else "no_flag"
                ),
            }
        )
    OUT.mkdir(parents=True, exist_ok=True)
    Table(rows=controls).write(OUT / "control_macs0417.ecsv", overwrite=True)
    Table(rows=targets).write(OUT / "predicted_macs0417.ecsv", overwrite=True)
    return {
        "inputs": {
            "jwst_cat": cat_path.name,
            "jwst_cat_sha256": hashlib.sha256(cat_path.read_bytes()).hexdigest(),
            "frame_offset_arcsec": [dra, ddec],
            "pixel_scale_arcsec": cat["scale"],
            "depth_mag": depth,
        },
        "gate_C2": c2,
        "control_classes": _count([c["class"] for c in controls]),
        "target_classes": _count([t["class"] for t in targets]),
        "flags": [
            f"{t['system']}@{t['ra']:.6f},{t['dec']:.6f}"
            for t in targets
            if t["reading"] != "no_flag"
        ],
        "flags_if_gate_ignored": sum(t["class"] == "absent" for t in targets),
        "rules": {
            "ELLIPSE_NSIGMA": ELLIPSE_NSIGMA,
            "MIN_REF_SNR": MIN_REF_SNR,
            "MAX_REF_MU": MAX_REF_MU,
            "MAX_REF_SPREAD": MAX_REF_SPREAD,
            "FLUX_RANGE": list(FLUX_RANGE),
            "TEN_SIGMA_DMAG": TEN_SIGMA_DMAG,
            "BASE_SEARCH": BASE_SEARCH,
            "C2": [C2_MIN_TESTABLE, C2_MIN_FRAC, C2_MIN_ARC],
            "ARC_MIN_ELL": ARC_MIN_ELL,
        },
    }


def _count(v) -> dict:
    return {c: v.count(c) for c in sorted(set(v))}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cache", type=Path, default=Path("data/cache/m0417x"))
    args = ap.parse_args(argv)
    summary = run(args.cache)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
