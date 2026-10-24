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
import bullet_gatek as bg
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
    model = lensmodel.MapLensModel(
        ax * k,
        ay * k,
        wcs,
        Z_LENS,
        FlatLambdaCDM(H0=COSMOLOGY[0], Om0=COSMOLOGY[1]),
        source=name,
        sha256=sha,
    )
    if not np.isclose(model.pixel_arcsec, PIXEL_ARCSEC, rtol=1e-6):  # the pixel reading's scale
        raise ValueError(f'{name}: map pixel {model.pixel_arcsec}" != {PIXEL_ARCSEC}"')
    return model


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


def verified_inputs(maps: Path) -> tuple[dict[str, str], str]:
    """sha256 of every input after its md5 matched the Zenodo record, and the model's sha256."""
    checks = record_checksums()
    shas = {f: md5_ok(maps / f, checks) for f in (*MAPS, "sourcepos.txt", "chires.txt")}
    return shas, hashlib.sha256((shas[MAPS[0]] + shas[MAPS[1]]).encode()).hexdigest()


def load_models(maps: Path, model_sha: str, readings=READINGS) -> dict:
    (ax, hx), (ay, _) = (read_zip_map(maps / f) for f in MAPS)
    return {
        r: build_model(ax, ay, hx, READINGS[r], model_sha, f"m1931-venus-v2.1 ({r})")
        for r in readings
    }


def run(maps: Path, out: Path = OUT, grid_step: float = GRID_STEP) -> dict:
    norm = k_dlds_check()
    if not norm["passes"]:
        raise SystemExit(f"error: k_dlds check fails: {norm}")
    shas, model_sha = verified_inputs(maps)
    src = parse_sourcepos((maps / "sourcepos.txt").read_text())
    images = load_images()
    spec = images[np.asarray(images["z_spec"], bool)]
    models = load_models(maps, model_sha)
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
        bt_spec = lensmodel.backtrace_images(model, spec, {})
        ip = lensmodel.imageplane_residuals(model, grid, bt_spec)
        lc._write(ip, out / "imageplane_spec.ecsv")
        ip_wls = au.imageplane(model, grid, bt_spec)
        nz = images[np.asarray(images["z"], float) > 0]
        ip_all = lensmodel.imageplane_residuals(
            model, grid, lensmodel.backtrace_images(model, nz, {})
        )
        systems = bg.per_system(ip)
        summary.update(
            {
                "grid_step_arcsec": grid_step,
                "grid_half_width_arcsec": half,
                "gate_K": mi.gate_k(ip),
                "reported": {
                    "gate_K_numbers_wls": mi.gate_k(ip_wls),
                    "gate_K_numbers_nonzero_z": mi.gate_k(ip_all),
                    "per_system_spec": systems,
                    "per_system_nonzero_z": bg.per_system(ip_all),
                    "extra_predicted_spec": [
                        s for s, v in systems.items() if v["n_predicted"] > v["n_catalogued"]
                    ],
                    "fewer_predicted_spec": [
                        s for s, v in systems.items() if v["n_predicted"] < v["n_catalogued"]
                    ],
                    "authors_chires_rms_arcsec_all_used": m19.chires_rms(
                        (maps / "chires.txt").read_text(), list(images["image_id"])
                    ),
                },
            }
        )
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def parse_chires_bary(text: str) -> dict[str, float]:
    """``chires.txt``: image -> the authors' Bary image-plane residual (arcsec)."""
    rows = [ln.split() for ln in text.splitlines()[1:]]
    return {r[0]: float(r[3]) for r in rows if len(r) == 4}


def posthoc(maps: Path, out: Path = OUT, grid_step: float = GRID_STEP) -> dict:
    """Post hoc (not in PLAN.md; after gate K): for the 44 non-zero-z images, the U distance of
    every system (its mixed per-image z as catalogued) and our image-plane residual next to the
    authors' ``chires.txt`` Bary residual. Run :func:`run` first (it checks the inputs)."""
    s = json.loads((out / "summary.json").read_text())
    reading = s["unit_check_U"]["reading"]
    if reading is None:
        raise SystemExit("error: unit check U did not pass; no reading to use")
    shas, model_sha = verified_inputs(maps)
    if shas != s["inputs_sha256"]:
        raise SystemExit("error: inputs differ from the plan run's (summary.json)")
    src = parse_sourcepos((maps / "sourcepos.txt").read_text())
    model = load_models(maps, model_sha, (reading,))[reading]
    images = load_images()
    nz = images[np.asarray(images["z"], float) > 0]
    dist = unit_distances(model, nz, src)
    grid = lensmodel.DeflectionGrid.compute(model, lc.grid_half_width(model, images), grid_step)
    ip = lensmodel.imageplane_residuals(model, grid, lensmodel.backtrace_images(model, nz, {}))
    bary = parse_chires_bary((maps / "chires.txt").read_text())
    ids = np.asarray(ip["image_id"], str)
    d = np.asarray(ip["dtheta_arcsec"], float)
    sysv = np.asarray(nz["system"], str)
    zs = np.asarray(nz["z"], float)
    spec_bary = np.array([bary.get(i, np.nan) for i in ids[np.asarray(nz["z_spec"], bool)]])
    allsys = np.asarray(images["system"], str)
    partial = sorted({k for k in sysv if (allsys == k).sum() > (sysv == k).sum()}, key=int)
    mixed = sorted({k for k in sysv if len(set(zs[sysv == k])) > 1}, key=int)
    res = {
        "label": "post hoc (not in PLAN.md)",
        "reading": reading,
        "u_distance_arcsec_all_systems": dist,
        "mixed_z_systems": mixed,
        "per_image": {
            i: {"ours_arcsec": float(v), "authors_bary_arcsec": bary.get(i)}
            for i, v in zip(ids, d, strict=True)
        },
        "rms_ours_without_mixed_z": mi.rms_of(d[~np.isin(sysv, mixed)]),
        "rms_authors_bary_spec": mi.rms_of(spec_bary),
        "n_spec": len(spec_bary),
        "n_spec_with_authors_value": int(np.isfinite(spec_bary).sum()),
        "partial_families": partial,
        "n_without_mixed_z": int(np.sum(~np.isin(sysv, mixed))),
    }
    (out / "posthoc.json").write_text(json.dumps(res, indent=1) + "\n")
    return res


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--maps", type=Path, required=True)
    ap.add_argument("--grid-step", type=float, default=GRID_STEP)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--posthoc", action="store_true", help="post hoc checks after the plan run")
    args = ap.parse_args(argv)
    if args.posthoc:
        r = posthoc(args.maps, args.out, args.grid_step)
        print(json.dumps(r["u_distance_arcsec_all_systems"], indent=1))
        return 0
    s = run(args.maps, args.out, args.grid_step)
    print(
        json.dumps(
            {k: s.get(k) for k in ("normalisation_check", "unit_check_U", "gate_K")}, indent=1
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
