"""MACS J1931.8-2635: unit check U and gate K of the VENUS AstroLensPy v2.1 deflection maps.

Frozen plan: ``results/lens_consistency/macs1931_gatek/PLAN.md``. The α normalisation (z_s 5.079)
is checked against the authors' ``k_dlds`` values before any map is read (``k_dlds_check``). U
decides between the README unit (arcsec) and the header unit (pixels, x 0.04") by the distance of
our per-system mean sources to the authors' ``sourcepos.txt`` (Bary columns); gate K is
``m0417_images.gate_k`` on the 26 spec-z images with the accepted reading. No extra-image, flux or
flag reading.

    python scripts/macs1931_gatek.py --maps DIR [--grid-step 0.25] [--out DIR]

``DIR`` holds ``Alpha_x.fits.zip``, ``Alpha_y.fits.zip``, ``sourcepos.txt`` and ``chires.txt``
from Zenodo 10.5281/zenodo.18646525 (md5 checked against the record before use).
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import a2744_uncover as au
import lens_consistency as lc
import m0417_images as mi
import macs1931_inputs as m19
import numpy as np
import requests
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

from jwst_anomaly import lensmodel

OUT = Path("results/lens_consistency/macs1931_gatek")
IMAGES = m19.OUT / m19.IMAGES
Z_LENS = 0.3525  # header Z_CL and README
COSMOLOGY = (70.0, 0.3)  # flat; reproduces the authors' k_dlds (k_dlds_check)
Z_NORM = 5.079  # README "DEFAULT SOURCE REDSHIFT"; the record: maps computed for zs = 5.079
# Readme_Magnifications_summary.txt: (z_s, k_dlds)
K_DLDS = ((2.0, 0.89043), (5.0, 0.99885), (9.0, 1.03257))
K_DLDS_TOL = 5e-5  # ASSUMPTION: the README prints five decimals
PIXEL_ARCSEC = 0.04  # README and header CDELT
READINGS = {"A_arcsec": 1.0, "P_pixels": PIXEL_ARCSEC}  # map value -> arcsec
U_MAX_MEDIAN_ARCSEC = 0.5  # ASSUMPTION (PLAN.md rule 1)
GRID_STEP = 0.25  # m0417_images default (PLAN.md rule 2)
MAPS = ("Alpha_x.fits.zip", "Alpha_y.fits.zip")
CHIRES_BARY = 0.4444  # ../macs1931_inputs/summary.json chires_rms_arcsec.Bary_RMSi


def dls_ds(z_s: float) -> float:
    c = FlatLambdaCDM(H0=COSMOLOGY[0], Om0=COSMOLOGY[1])
    return float(
        (c.angular_diameter_distance(Z_LENS, z_s) / c.angular_diameter_distance(z_s)).value
    )


def k_dlds_check() -> dict:
    """(D_LS/D_S)(z) / (D_LS/D_S)(Z_NORM) against the authors' printed ``k_dlds``."""
    rows = [
        {"z_s": z, "authors": k, "ours": round(dls_ds(z) / dls_ds(Z_NORM), 6)} for z, k in K_DLDS
    ]
    return {
        "z_norm": Z_NORM,
        "rows": rows,
        "passes": all(abs(r["ours"] - r["authors"]) <= K_DLDS_TOL for r in rows),
    }


def md5_ok(path: Path, record_files: dict[str, str]) -> str:
    """sha256 of ``path`` after its md5 matched the Zenodo record (raises otherwise)."""
    md5, sha = hashlib.md5(), hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 22):
            md5.update(chunk)
            sha.update(chunk)
    want = record_files[path.name]
    if f"md5:{md5.hexdigest()}" != want:
        raise ValueError(f"{path.name}: md5 {md5.hexdigest()} != record {want}")
    return sha.hexdigest()


def record_checksums() -> dict[str, str]:
    r = requests.get(m19.RECORD, timeout=60)
    r.raise_for_status()
    return {f["key"]: f["checksum"] for f in r.json()["files"]}


def read_zip_map(path: Path) -> tuple[np.ndarray, object]:
    """The single FITS member of a map zip, its one 2-D HDU: (float64 data, header)."""
    from astropy.io import fits

    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if n.lower().endswith(".fits")]
        if len(names) != 1:
            raise ValueError(f"{path.name}: expected one FITS member, found {names}")
        blob = z.read(names[0])
    with fits.open(io.BytesIO(blob)) as hdul:
        hdus = [h for h in hdul if h.data is not None and h.data.ndim == 2]
        if len(hdus) != 1:  # the VENUS maps sit in one image extension after an empty primary
            raise ValueError(f"{path.name}: expected one 2-D HDU, found {len(hdus)}")
        return np.asarray(hdus[0].data, np.float64), hdus[0].header.copy()


def build_model(ax, ay, header, scale: float, sha: str, name: str) -> lensmodel.MapLensModel:
    """Map values x ``scale`` (arcsec) / (D_LS/D_S)(Z_NORM): deflection at D_LS/D_S = 1."""
    import warnings

    from astropy.wcs import WCS, FITSFixedWarning

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(header)
    k = scale / dls_ds(Z_NORM)
    return lensmodel.MapLensModel(
        ax * k,
        ay * k,
        wcs,
        Z_LENS,
        FlatLambdaCDM(H0=COSMOLOGY[0], Om0=COSMOLOGY[1]),
        source=name,
        sha256=sha,
    )


def parse_sourcepos(text: str) -> dict[str, tuple[float, float]]:
    """``sourcepos.txt``: system -> (Bary_Pos_x, Bary_Pos_y) in 1-based map pixels."""
    out = {}
    for line in text.splitlines():
        f = line.split()
        if len(f) == 7 and f[0].isdigit():
            out[f[0]] = (float(f[5]), float(f[6]))
    if not out:
        raise ValueError("sourcepos.txt: no rows")
    return out


def authors_sources(model: lensmodel.MapLensModel, src: dict) -> dict[str, tuple[float, float]]:
    """The authors' source pixels in our model frame (arcsec), through the map WCS."""
    out = {}
    for s, (i1, j1) in src.items():
        ra, dec = model.wcs.pixel_to_world_values(i1 - 1.0, j1 - 1.0)
        x, y = model.to_frame(np.atleast_1d(float(ra)), np.atleast_1d(float(dec)))
        out[s] = (float(x[0]), float(y[0]))
    return out


def unit_distances(model, images: Table, src: dict) -> dict[str, float]:
    """Per system: our mean back-traced source to the authors' Bary source (arcsec)."""
    bt = lensmodel.backtrace_images(model, images, {})
    ref = authors_sources(model, src)
    sysv = np.asarray(bt["system"], str)
    out = {}
    for s in dict.fromkeys(sysv):
        m = sysv == s
        bx, by = float(np.mean(bt["beta_x"][m])), float(np.mean(bt["beta_y"][m]))
        out[s] = float(np.hypot(bx - ref[s][0], by - ref[s][1])) if s in ref else float("nan")
    return out


def unit_check(distances: dict[str, dict[str, float]]) -> dict:
    """PLAN.md rule 1: a reading is accepted when its median per-system distance is
    <= U_MAX_MEDIAN_ARCSEC; U passes when exactly one reading is accepted."""
    med = {
        r: float(np.median(list(d.values())))
        if d and np.all(np.isfinite(list(d.values())))
        else float("nan")
        for r, d in distances.items()
    }
    accepted = [r for r, v in med.items() if np.isfinite(v) and v <= U_MAX_MEDIAN_ARCSEC]
    return {
        "median_arcsec": med,
        "accepted": accepted,
        "passes": len(accepted) == 1,
        "reading": accepted[0] if len(accepted) == 1 else None,
    }


def load_images(path: Path = IMAGES) -> Table:
    t = Table.read(path, format="ascii.ecsv")
    t = t[np.asarray(t["used"], bool)]
    t["image_id"] = np.asarray(t["image_id"], str)
    t["system"] = np.asarray(t["system"], str)
    t["a"] = np.ones(len(t))
    t.meta.setdefault("source", str(path))
    return t


def per_system(ip: Table) -> dict:
    out = {}
    sysv = np.asarray(ip["system"], str)
    d = np.asarray(ip["dtheta_arcsec"], float)
    for s in dict.fromkeys(sysv):
        m = sysv == s
        out[s] = {
            "rms_arcsec": mi.rms_of(d[m]),
            "n_catalogued": int(m.sum()),
            "n_predicted": int(np.max(ip["n_predicted"][m])),
        }
    return out


def run(maps: Path, out: Path = OUT, grid_step: float = GRID_STEP) -> dict:
    norm = k_dlds_check()
    if not norm["passes"]:
        raise SystemExit(f"error: k_dlds check fails: {norm}")
    checks = record_checksums()
    shas = {f: md5_ok(maps / f, checks) for f in (*MAPS, "sourcepos.txt", "chires.txt")}
    src = parse_sourcepos((maps / "sourcepos.txt").read_text())
    (ax, hx), (ay, _) = (read_zip_map(maps / f) for f in MAPS)
    images = load_images()
    spec = images[np.asarray(images["z_spec"], bool)]
    model_sha = hashlib.sha256((shas[MAPS[0]] + shas[MAPS[1]]).encode()).hexdigest()
    models = {
        r: build_model(ax, ay, hx, s, model_sha, f"m1931-venus-v2.1 ({r})")
        for r, s in READINGS.items()
    }
    del ax, ay
    dist = {r: unit_distances(m, spec, src) for r, m in models.items()}
    u = unit_check(dist)
    summary = {
        "plan": str(out / "PLAN.md"),
        "inputs_sha256": shas,
        "normalisation_check": norm,
        "unit_check_U": {**u, "per_system_arcsec": dist},
    }
    if u["passes"]:
        model = models[u["reading"]]
        half = lc.grid_half_width(model, images)
        grid = lensmodel.DeflectionGrid.compute(model, half, grid_step)
        ip = lensmodel.imageplane_residuals(
            model, grid, lensmodel.backtrace_images(model, spec, {})
        )
        lc._write(ip, out / "imageplane_spec.ecsv")
        bt_spec = lensmodel.backtrace_images(model, spec, {})
        ip_wls = au.imageplane(model, grid, bt_spec)
        nz = images[np.asarray(images["z"], float) > 0]
        ip_all = lensmodel.imageplane_residuals(
            model, grid, lensmodel.backtrace_images(model, nz, {})
        )
        systems = per_system(ip)
        summary.update(
            {
                "grid_step_arcsec": grid_step,
                "grid_half_width_arcsec": half,
                "gate_K": mi.gate_k(ip),
                "reported": {
                    "gate_K_numbers_wls": mi.gate_k(ip_wls),
                    "gate_K_numbers_nonzero_z": mi.gate_k(ip_all),
                    "per_system_spec": systems,
                    "per_system_nonzero_z": per_system(ip_all),
                    "extra_predicted_spec": [
                        s for s, v in systems.items() if v["n_predicted"] > v["n_catalogued"]
                    ],
                    "fewer_predicted_spec": [
                        s for s, v in systems.items() if v["n_predicted"] < v["n_catalogued"]
                    ],
                    "authors_chires_bary_rms_arcsec": CHIRES_BARY,
                },
            }
        )
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--maps", type=Path, required=True)
    ap.add_argument("--grid-step", type=float, default=GRID_STEP)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    s = run(args.maps, args.out, args.grid_step)
    print(
        json.dumps(
            {k: s.get(k) for k in ("normalisation_check", "unit_check_U", "gate_K")}, indent=1
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
