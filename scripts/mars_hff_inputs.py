"""MARS HFF free-form maps (Cha & Jee 2023): input check C / U and gate K on six fields.

Frozen plan: ``results/lens_consistency/mars_hff_inputs/PLAN.md``. Zenodo 10.5281/zenodo.7804575
(CC-BY-4.0) releases one zip with 100 x 100 deflection maps per Frontier Field (D_ds / D_s = 1) and
one image catalogue (RA, Dec, z, class, Delta_rms). Per field: C (coverage), U (arcsec vs pixel
reading at catalogue z), K (``m0417_images.gate_k``, mean source), K-WLS reported. No extra-image,
flux or flag reading.

    python scripts/mars_hff_inputs.py --zip DIR

``DIR`` receives the zip (md5 checked against the record); it is deleted after the run unless
``--keep``.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import warnings
import zipfile
from pathlib import Path

import a2744_uncover as au
import bullet_gatek as bg
import lens_consistency as lc
import m0417_images as mi
import macs1931_inputs as m19
import numpy as np
import requests
import spt0546_inputs as sp
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning

from jwst_anomaly import lensmodel

OUT = Path("results/lens_consistency/mars_hff_inputs")
RECORD = "https://zenodo.org/api/records/7804575"
ZIP = "MARS_HFF_clusters_result.zip"
ZIP_URL = f"https://zenodo.org/api/records/7804575/files/{ZIP}/content"
ROOT = "MARS_HFF_clusters_result/"
CATALOGUE = ROOT + "MARS_HFF_Multiple_image_catalog.txt"
COSMOLOGY = (70.0, 0.27)  # flat; the authors' Function_torch.py (omegam 0.27, vac 0.73)
# catalogue section name -> (zip folder, map tag, z_lens from load_parameter_torch.py)
FIELDS = {
    "Abell 370": ("Abell 370", "A370", 0.375),
    "Abell 2744": ("Abell 2744", "A2744", 0.308),
    "Abell S1063": ("Abell S1063", "AS1063", 0.348),
    "MACS J0416.1-2403": ("MACS J0416.1-2403", "MACSJ0416", 0.396),
    "MACS J0717.5+3745": ("MACS J0717.5+3745", "MACSJ0717", 0.545),
    "MACS J1149.5-2223": ("MACS J1149.5+2223", "MACSJ1149", 0.543),
}
C_MIN_INSIDE = 0.9  # ASSUMPTION (PLAN.md rule 1)
U_MAX_MEDIAN_ARCSEC = 0.5  # ASSUMPTION (PLAN.md rule 2)
MIN_SPEC_PER_SYSTEM = 2
GRID_STEP = 0.25


def parse_catalogue(text: str) -> Table:
    """All rows with their field (from the ``# Multiple image(s) of <field>`` headers); ``fitted_z``
    marks a ``^{**}`` redshift (fitted by MARS)."""
    rows, field = [], None
    for line in text.splitlines():
        s = line.strip()
        m = re.match(r"#\s*Multiple images? of (.+?)\s*$", s)
        if m:
            field = m.group(1)
            continue
        if not s or s.startswith("#"):
            continue
        f = s.split()
        if len(f) < 6 or field is None:
            raise ValueError(f"catalogue row without field or too short: {line!r}")
        ztxt = f[3]
        rows.append(
            (field, f[0], f[0].split(".")[0], float(f[1]), float(f[2]),
             float(ztxt.split("^")[0]), "**" in ztxt, f[4], float(f[5]))
        )  # fmt: skip
    if not rows:
        raise ValueError("catalogue: no rows")
    t = Table(
        rows=rows,
        names=("field", "image_id", "system", "ra", "dec", "z", "fitted_z", "cls", "delta_rms"),
    )
    for c in ("field", "image_id", "system", "cls"):
        t[c] = np.asarray(t[c], str)
    return t


def spec_images(cat: Table, field: str) -> Table:
    """PLAN.md: gold rows of ``field`` without a fitted z, in systems with >= 2 such images."""
    t = cat[
        (np.asarray(cat["field"]) == field)
        & (np.asarray(cat["cls"]) == "gold")
        & ~np.asarray(cat["fitted_z"], bool)
    ]
    sysv = np.asarray(t["system"], str)
    keep = np.array([np.sum(sysv == s) >= MIN_SPEC_PER_SYSTEM for s in sysv], bool)
    t = t[keep] if len(t) else t
    t["a"] = np.ones(len(t))
    return t


def coverage_field(model, images: Table, min_inside: float = C_MIN_INSIDE) -> dict:
    """PLAN.md rule 1: fraction of images inside the map; passes when >= ``min_inside``."""
    x, y = model.to_frame(np.asarray(images["ra"], float), np.asarray(images["dec"], float))
    ax, _ = model.deflection_xy(x, y)
    inside = np.isfinite(np.asarray(ax, float))
    frac = float(inside.mean()) if len(inside) else 0.0
    return {
        "n": len(images),
        "inside": int(inside.sum()),
        "fraction": frac,
        "outside": [str(i) for i in images["image_id"][~inside]],
        "passes": bool(len(inside) and frac >= min_inside),
        "mask": inside,
    }


def record_outcome(fields: dict[str, dict]) -> dict:
    """PLAN.md rule 4: usable fields are those whose C, U and K all pass."""
    ok = [
        f
        for f, d in fields.items()
        if d.get("check_C", {}).get("passes")
        and d.get("check_U", {}).get("passes")
        and d.get("gate_K", {}).get("passes")
    ]
    return {"usable_fields": ok, "n_usable": len(ok), "n_fields": len(fields)}


def build_model(ax, ay, header, scale_px: bool, z_lens: float, sha: str, name: str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(header)
    cosmo = FlatLambdaCDM(H0=COSMOLOGY[0], Om0=COSMOLOGY[1])
    probe = lensmodel.MapLensModel(ax, ay, wcs, z_lens, cosmo, source=name, sha256=sha)
    k = probe.pixel_arcsec if scale_px else 1.0
    return lensmodel.MapLensModel(ax * k, ay * k, wcs, z_lens, cosmo, source=name, sha256=sha)


def _fits(zf: zipfile.ZipFile, member: str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with fits.open(io.BytesIO(zf.read(member))) as hdul:
            return np.asarray(hdul[0].data, float), hdul[0].header.copy()


def fetch_zip(session: requests.Session, dest: Path) -> dict:
    record = m19._get(session, RECORD).json()
    md5 = {f["key"]: f["checksum"].removeprefix("md5:") for f in record["files"]}[ZIP]
    dest.mkdir(parents=True, exist_ok=True)
    p = dest / ZIP
    if not p.exists():
        p.write_bytes(m19._get(session, ZIP_URL).content)
    blob = p.read_bytes()
    if hashlib.md5(blob).hexdigest() != md5:
        raise RuntimeError(f"{ZIP}: md5 differs from the record")
    return {"bytes": len(blob), "md5": md5, "sha256": hashlib.sha256(blob).hexdigest()}


def per_system_u(model, images: Table) -> dict[str, float | None]:
    d = np.asarray(lensmodel.backtrace_images(model, images, {})["dtheta_arcsec"], float)
    sysv = np.asarray(images["system"], str)
    out = {}
    for s in dict.fromkeys(sysv):
        r = mi.rms_of(d[sysv == s])
        out[s] = r if np.isfinite(r) else None
    return out


def field_block(zf, field: str, cat: Table, out: Path, grid_step: float) -> dict:
    folder, tag, z_lens = FIELDS[field]
    base = f"{ROOT}{folder}/deflection_angle_{tag}_"
    (ax, hx), (ay, _) = (_fits(zf, base + f"{c}_w_header.fits") for c in "xy")
    sha = hashlib.sha256(zf.read(base + "x_w_header.fits") + zf.read(base + "y_w_header.fits"))
    sha = sha.hexdigest()
    models = {
        r: build_model(ax, ay, hx, r == "P_pixels", z_lens, sha, f"mars-{tag} ({r})")
        for r in ("A_arcsec", "P_pixels")
    }
    images = spec_images(cat, field)
    res = {
        "z_lens": z_lens,
        "map_shape": list(models["A_arcsec"].shape),
        "pixel_arcsec": round(models["A_arcsec"].pixel_arcsec, 4),
        "bunit": hx.get("BUNIT"),
        "spec_images": len(images),
        "spec_systems": len(set(images["system"])),
    }
    c = coverage_field(models["A_arcsec"], images)
    mask = c.pop("mask")
    res["check_C"] = c
    if not c["passes"]:
        return res
    images = images[mask]
    vals = {r: per_system_u(m, images) for r, m in models.items()}
    u = sp.unit_check(vals, U_MAX_MEDIAN_ARCSEC)
    res["check_U"] = {**u, "per_system": vals}
    if not u["passes"]:
        return res
    model = models[u["reading"]]
    grid = lensmodel.DeflectionGrid.compute(model, lc.grid_half_width(model, images), grid_step)
    bt = lensmodel.backtrace_images(model, images, {})
    ip = lensmodel.imageplane_residuals(model, grid, bt)
    ipw = au.imageplane(model, grid, bt)
    lc._write(ip, out / f"imageplane_{tag}.ecsv")
    res["gate_K"] = mi.gate_k(ip)
    res["reported"] = {
        "K_wls_not_a_gate": mi.gate_k(ipw),
        "per_system": bg.per_system(ip),
        "authors_median_delta_rms_arcsec": float(np.median(images["delta_rms"])),
    }
    return res


def run(zdir: Path, out: Path = OUT, grid_step: float = GRID_STEP) -> dict:
    meta = fetch_zip(requests.Session(), zdir)
    out.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zdir / ZIP) as zf:
        cat = parse_catalogue(zf.read(CATALOGUE).decode())
        missing = sorted(set(cat["field"]) - set(FIELDS))
        if missing:
            raise ValueError(f"catalogue fields without a map entry: {missing}")
        fields = {f: field_block(zf, f, cat, out, grid_step) for f in FIELDS}
    summary = {
        "plan": str(out / "PLAN.md"),
        "zip": meta,
        "cosmology": {"H0": COSMOLOGY[0], "Om0": COSMOLOGY[1]},
        "catalogue_rows": len(cat),
        "grid_step_arcsec": grid_step,
        "fields": fields,
        "record_outcome": record_outcome(fields),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--zip", type=Path, required=True, help="directory for the zip")
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--grid-step", type=float, default=GRID_STEP)
    ap.add_argument("--keep", action="store_true", help="keep the downloaded zip")
    a = ap.parse_args(argv)
    try:
        s = run(a.zip, a.out, a.grid_step)
    finally:
        if not a.keep:
            (a.zip / ZIP).unlink(missing_ok=True)
    brief = {
        f: {k: d.get(k) for k in ("check_C", "gate_K")} | {"U": d.get("check_U", {}).get("reading")}
        for f, d in s["fields"].items()
    }
    print(json.dumps({"fields": brief, "outcome": s["record_outcome"]}, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
