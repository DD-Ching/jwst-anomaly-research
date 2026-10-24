"""SPT-CL J0546-5345: input check, unit check U and free-z gate K of the LTM deflection maps.

Frozen plan: ``results/lens_consistency/spt0546_inputs/PLAN.md``. The Allingham+2025 LTM model
(Zenodo 10.5281/zenodo.16534337, CC-BY-4.0) releases alpha_x / alpha_y at z_s = 3.5, a Lenstool
format image list without redshifts and a spec-z catalogue. Checks: C (coverage), Z (spec-z),
U (arcsec vs pixel reading, per-system z free), K (``m0417_images.gate_k`` at fitted z), R (a DJA v7
imaging field). No extra-image, flux or flag reading.

    python scripts/spt0546_inputs.py --maps DIR [--grid-step 0.25]

``DIR`` receives the two maps and the two catalogues (md5 checked against the record); the maps are
deleted after the run unless ``--keep``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import warnings
from pathlib import Path

import bullet_gatek as bg
import lens_consistency as lc
import m0417_images as mi
import macs1931_inputs as m19
import numpy as np
import requests
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.wcs import WCS, FITSFixedWarning

from jwst_anomaly import lensmodel

OUT = Path("results/lens_consistency/spt0546_inputs")
RECORD = "https://zenodo.org/api/records/16534337"
FILE_URL = "https://zenodo.org/records/16534337/files/{}?download=1"
MAPS = (
    "alpha_x_sptcl0546_z3p5_20250611_LTM_Gc.fits",
    "alpha_y_sptcl0546_z3p5_20250611_LTM_Gc.fits",
)
MUL = "actcl0546_mulfile_v1.mul"
SPECZ = "actcl0546_spec-z.cat"
Z_LENS = 1.07  # record description
Z_NORM = 3.5  # record: every z-dependent map is computed for z_s = 3.5
COSMOLOGY = (70.0, 0.3)  # flat; ASSUMPTION (PLAN.md)
PIXEL_ARCSEC = 0.04  # header |CDELT| 1.111e-5 deg
READINGS = {"A_arcsec": 1.0, "P_pixels": PIXEL_ARCSEC}  # map value -> arcsec
U_MAX_MEDIAN_ARCSEC = 0.5  # ASSUMPTION (PLAN.md rule 3)
SPECZ_MATCH_ARCSEC = 0.5  # PLAN.md rule 2
SPECZ_CONF = (3, 4)
SPECZ_MIN_DZ = 0.1  # z > z_lens + 0.1
GRID_STEP = 0.25
DJA_NAMES = ("spt0546", "spt-0546", "j0546", "sptcl0546", "actcl0546")


def z_grid() -> np.ndarray:
    """The ``m0417_images.fitted_z`` grid: 60 log-spaced z from z_lens + 0.2 to 10."""
    return np.geomspace(Z_LENS + mi.Z_MIN_OFFSET, mi.Z_MAX, mi.Z_STEPS)


def cosmology() -> FlatLambdaCDM:
    return FlatLambdaCDM(H0=COSMOLOGY[0], Om0=COSMOLOGY[1])


def dls_ds(z_s: float) -> float:
    c = cosmology()
    return float(1.0 - c.comoving_transverse_distance(Z_LENS) / c.comoving_transverse_distance(z_s))


def parse_mul(text: str) -> Table:
    """Lenstool-format image list: id, ra, dec (the z column is 0 throughout, i.e. none).
    ``candidate`` marks ``c``-prefixed ids; ``system`` is the id before the first dot, without
    the ``c``."""
    rows = []
    for line in text.splitlines():
        f = line.split()
        if not f or f[0].startswith("#"):
            continue
        if len(f) < 3:
            raise ValueError(f"mul row too short: {line!r}")
        iid = f[0]
        rows.append(
            (iid, iid.removeprefix("c").split(".")[0], float(f[1]), float(f[2]), iid[0] == "c")
        )
    if not rows:
        raise ValueError("mul file: no rows")
    t = Table(rows=rows, names=("image_id", "system", "ra", "dec", "candidate"))
    t["image_id"] = np.asarray(t["image_id"], str)
    t["system"] = np.asarray(t["system"], str)
    return t


def parse_specz(text: str) -> Table:
    rows = [
        [float(v) for v in ln.split()[:9]]
        for ln in text.splitlines()
        if ln.strip() and not ln.startswith("#")
    ]
    t = Table(rows=rows, names=("id", "ra", "dec", "a", "b", "theta", "z", "mag", "conf"))
    t["conf"] = np.asarray(t["conf"], int)
    return t


def _sep_arcsec(ra1, dec1, ra2, dec2):
    cosd = np.cos(np.radians(0.5 * (dec1 + dec2)))
    return np.hypot((ra1 - ra2) * cosd, dec1 - dec2) * 3600.0


def specz_systems(images: Table, cat: Table) -> dict:
    """PLAN.md rule 2: systems with a main image <= 0.5" from a conf 3 / 4 entry behind the lens."""
    good = cat[np.isin(cat["conf"], SPECZ_CONF) & (np.asarray(cat["z"]) > Z_LENS + SPECZ_MIN_DZ)]
    hits = {}
    for r in images[~np.asarray(images["candidate"], bool)]:
        if not len(good):
            break
        d = _sep_arcsec(float(r["ra"]), float(r["dec"]), np.asarray(good["ra"]), good["dec"])
        k = int(np.argmin(d))
        if d[k] <= SPECZ_MATCH_ARCSEC:
            hits.setdefault(str(r["system"]), []).append(
                {"image": str(r["image_id"]), "z": float(good["z"][k]), "sep": float(d[k])}
            )
    return {"entries_behind_lens_conf34": len(good), "systems": hits, "n_systems": len(hits)}


def coverage(model: lensmodel.MapLensModel, images: Table) -> dict:
    """PLAN.md rule 1: every main image inside the map (deflection finite there)."""
    x, y = model.to_frame(np.asarray(images["ra"], float), np.asarray(images["dec"], float))
    ax, _ = model.deflection_xy(x, y)
    inside = np.isfinite(np.asarray(ax, float))
    return {
        "n": len(images),
        "inside": int(inside.sum()),
        "outside": [str(i) for i in images["image_id"][~inside]],
        "passes": bool(inside.all()),
    }


def unit_minima(model, images: Table, zs: np.ndarray) -> dict[str, dict]:
    """Per system: the smallest rms (over the z grid) of the backtrace ``dtheta_arcsec``."""
    sysv = np.asarray(images["system"], str)
    per_z = []
    for z in zs:
        t = images.copy()
        t["z"] = z
        per_z.append(np.asarray(lensmodel.backtrace_images(model, t, {})["dtheta_arcsec"], float))
    out = {}
    for s in dict.fromkeys(sysv):
        m = sysv == s
        r = np.array([mi.rms_of(d[m]) for d in per_z])
        r = np.where(np.isfinite(r), r, np.inf)
        k = int(np.argmin(r))
        out[s] = {"min_rms_arcsec": float(r[k]), "z_at_min": float(zs[k])}
    return out


def unit_check(minima: dict[str, dict[str, float]]) -> dict:
    """PLAN.md rule 3: a reading is accepted when the median of its per-system minima is
    <= U_MAX_MEDIAN_ARCSEC; U passes when exactly one reading is accepted."""
    med = {}
    for r, d in minima.items():
        v = np.array(list(d.values()), float)
        med[r] = float(np.median(v)) if v.size and np.all(np.isfinite(v)) else float("nan")
    accepted = [r for r, v in med.items() if np.isfinite(v) and v <= U_MAX_MEDIAN_ARCSEC]
    return {
        "median_arcsec": med,
        "accepted": accepted,
        "passes": len(accepted) == 1,
        "reading": accepted[0] if len(accepted) == 1 else None,
    }


def fitted_z(model, grid, images: Table, zs: np.ndarray) -> Table:
    """Exact image-plane residuals of each system at its fitted z (``m0417_images.best_z``)."""
    per_z = []
    for z in zs:
        t = images.copy()
        t["z"] = z
        per_z.append(
            lensmodel.imageplane_residuals(model, grid, lensmodel.backtrace_images(model, t, {}))
        )
    out = []
    for fam in dict.fromkeys(np.asarray(images["system"], str)):
        sel = [np.asarray(ip["system"], str) == fam for ip in per_z]
        k = mi.best_z(
            [np.asarray(ip["dtheta_arcsec"], float)[s] for ip, s in zip(per_z, sel, strict=True)]
        )
        rows = per_z[k][sel[k]]
        rows["z_fit"] = np.full(len(rows), zs[k])
        rows["z_fit_on_edge"] = np.full(len(rows), k in (0, len(zs) - 1))
        out.append(rows)
    return vstack(out)


def candidate_offsets(model, grid, main_ip: Table, images: Table) -> dict:
    """Reported: each ``c`` candidate's distance to the nearest image predicted for its system's
    mean source at the system's fitted z (systems with main images only)."""
    out = {}
    cands = images[np.asarray(images["candidate"], bool)]
    sysv = np.asarray(main_ip["system"], str)
    for r in cands:
        m = sysv == str(r["system"])
        if not m.any():
            continue
        z = float(main_ip["z_fit"][m][0])
        pred = lensmodel.find_images(
            model,
            grid,
            float(np.mean(main_ip["beta_x"][m])),
            float(np.mean(main_ip["beta_y"][m])),
            z,
        )
        ox, oy = model.to_frame(np.atleast_1d(float(r["ra"])), np.atleast_1d(float(r["dec"])))
        d = np.hypot(np.asarray(pred["x"]) - ox[0], np.asarray(pred["y"]) - oy[0])
        out[str(r["image_id"])] = {
            "z_fit": z,
            "nearest_predicted_arcsec": float(d.min()) if len(d) else float("nan"),
        }
    return out


def dja_check(index_html: str) -> dict:
    fields = sorted(set(re.findall(r"([a-z0-9][a-z0-9_-]*?)-grizli-v7", index_html.lower())))
    if not fields:
        raise ValueError("no field parsed from the DJA index")
    hits = [f for f in fields if any(n in f for n in DJA_NAMES)]
    return {
        "index": m19.DJA_INDEX,
        "fields_listed": len(fields),
        "names_searched": list(DJA_NAMES),
        "matches": hits,
        "passes": bool(hits),
    }


def fetch_inputs(session: requests.Session, maps: Path) -> dict[str, str]:
    """Download (if absent) and md5-check every input against the record; sha256 per file."""
    record = m19._get(session, RECORD).json()
    md5 = {f["key"]: f["checksum"].removeprefix("md5:") for f in record["files"]}
    maps.mkdir(parents=True, exist_ok=True)
    shas = {}
    for name in (*MAPS, MUL, SPECZ):
        p = maps / name
        if not p.exists():
            p.write_bytes(m19._get(session, FILE_URL.format(name)).content)
        blob = p.read_bytes()
        if hashlib.md5(blob).hexdigest() != md5[name]:
            raise RuntimeError(f"{name}: md5 differs from the record")
        shas[name] = hashlib.sha256(blob).hexdigest()
    return shas


def read_map(path: Path) -> tuple[np.ndarray, fits.Header]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # non-standard CDELT cards (fixed by astropy)
        with fits.open(path) as hdul:
            return np.asarray(hdul[0].data, float), hdul[0].header.copy()


def build_model(ax, ay, header, scale: float, sha: str, name: str) -> lensmodel.MapLensModel:
    """Map values x ``scale`` (arcsec) / (D_LS/D_S)(Z_NORM): deflection at D_LS/D_S = 1."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(header)
    k = scale / dls_ds(Z_NORM)
    model = lensmodel.MapLensModel(
        ax * k, ay * k, wcs, Z_LENS, cosmology(), source=name, sha256=sha
    )
    if not np.isclose(model.pixel_arcsec, PIXEL_ARCSEC, rtol=1e-6):
        raise ValueError(f'{name}: map pixel {model.pixel_arcsec}" != {PIXEL_ARCSEC}"')
    return model


def run(maps: Path, out: Path = OUT, grid_step: float = GRID_STEP) -> dict:
    session = requests.Session()
    shas = fetch_inputs(session, maps)
    images = parse_mul((maps / MUL).read_text())
    main = images[~np.asarray(images["candidate"], bool)]
    main["a"] = np.ones(len(main))
    cat = parse_specz((maps / SPECZ).read_text())
    (ax, hx), (ay, _) = (read_map(maps / f) for f in MAPS)
    model_sha = hashlib.sha256((shas[MAPS[0]] + shas[MAPS[1]]).encode()).hexdigest()
    models = {
        r: build_model(ax, ay, hx, s, model_sha, f"spt0546-ltm ({r})") for r, s in READINGS.items()
    }
    zs = z_grid()
    summary = {
        "plan": str(out / "PLAN.md"),
        "inputs_sha256": shas,
        "cosmology": {"H0": COSMOLOGY[0], "Om0": COSMOLOGY[1], "z_lens": Z_LENS, "z_norm": Z_NORM},
        "image_list": {
            "rows": len(images),
            "main_images": len(main),
            "main_systems": len(set(main["system"])),
            "candidate_rows": int(np.sum(images["candidate"])),
        },
        "map_header": {
            "shape": list(ax.shape),
            "pixel_arcsec": round(models["A_arcsec"].pixel_arcsec, 6),
            "ctype": [hx["CTYPE1"], hx["CTYPE2"]],
            "bunit": hx.get("BUNIT"),
        },
        "check_C": coverage(models["A_arcsec"], main),
        "check_Z": specz_systems(images, cat),
        "z_grid": [float(zs[0]), float(zs[-1]), len(zs)],
    }
    if summary["check_C"]["passes"]:
        minima = {r: unit_minima(m, main, zs) for r, m in models.items()}
        u = unit_check(
            {r: {s: v["min_rms_arcsec"] for s, v in d.items()} for r, d in minima.items()}
        )
        summary["check_U"] = {**u, "per_system": minima}
        if u["passes"]:
            model = models[u["reading"]]
            half = lc.grid_half_width(model, images)
            grid = lensmodel.DeflectionGrid.compute(model, half, grid_step)
            ip = fitted_z(model, grid, main, zs)
            lc._write(ip, out / "imageplane_fitz.ecsv")
            systems = bg.per_system(ip)
            for s in systems:
                m = np.asarray(ip["system"], str) == s
                systems[s]["z_fit"] = float(ip["z_fit"][m][0])
                systems[s]["z_fit_on_edge"] = bool(ip["z_fit_on_edge"][m][0])
            summary.update(
                {
                    "grid_step_arcsec": grid_step,
                    "grid_half_width_arcsec": half,
                    "gate_K_free_z": mi.gate_k(ip),
                    "reported": {
                        "per_system": systems,
                        "candidates": candidate_offsets(model, grid, ip, images),
                    },
                }
            )
    summary["check_R"] = dja_check(m19._get(session, m19.DJA_INDEX).text)
    out.mkdir(parents=True, exist_ok=True)
    img = images.copy()
    img.meta = {
        "provenance": "observed (positions, as listed by the model authors); no redshifts",
        "source": FILE_URL.format(MUL),
    }
    img.write(out / "spt0546_images.ecsv", format="ascii.ecsv", overwrite=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--maps", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--grid-step", type=float, default=GRID_STEP)
    ap.add_argument("--keep", action="store_true", help="keep the downloaded maps")
    a = ap.parse_args(argv)
    try:
        s = run(a.maps, a.out, a.grid_step)
    finally:
        if not a.keep:
            for f in MAPS:
                (a.maps / f).unlink(missing_ok=True)
    print(json.dumps({k: v for k, v in s.items() if k != "reported"}, indent=1)[:4000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
