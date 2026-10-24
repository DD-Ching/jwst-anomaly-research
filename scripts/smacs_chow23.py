"""SMACS 0723: input check, unit / x-sign check U and fixed-z gate K of the Chow+2023 alpha maps.

Frozen plan: ``results/lens_consistency/smacs_chow23/PLAN.md``. The Chow et al. 2023 GLAFIC model
(Zenodo 10.5281/zenodo.10158248, arXiv:2310.09790) releases alpha maps at D_ds / D_s = 1 on a 0.06"
CAR grid; only the two alpha files are fetched (in parallel). Checks: C (coverage), U (arcsec vs
pixel reading, alpha_x as given vs negated), K (``m0417_images.gate_k``, mean source), K-WLS
(reported).
Reported: per-system counts for Chow and for ICLv2 computed the same way, and the systems whose
predicted counts differ. Image list: the fixed-z systems of the Mahler+22 ICLv2 list. No flux, flag
or extra-image reading.

    python scripts/smacs_chow23.py --maps DIR [--grid-step 0.25]

The maps are deleted after the run unless ``--keep``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import lens_consistency as lc
import numpy as np
import requests
import s1063_wslap as sw
import spt0546_inputs as sp
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning

from jwst_anomaly import lensmodel

OUT = Path("results/lens_consistency/smacs_chow23")
IMAGES = Path("results/lens_consistency/smacs_djaref/mahler22_images.ecsv")
RECORD = "https://zenodo.org/api/records/10158248"
MAPS = ("Chow23+_LensModel_alpha_x.fits", "Chow23+_LensModel_alpha_y.fits")
Z_LENS = 0.39  # as smacs_djaref (ICLv2 best.par); ASSUMPTION for this model (PLAN.md)
COSMOLOGY = (70.0, 0.3)  # flat; ASSUMPTION (PLAN.md), only D_ds / D_s enters
PIXEL_ARCSEC = 0.06  # header |CDELT| 1.66667e-5 deg
# reading -> (map value -> arcsec, sign applied to alpha_x)
READINGS = {
    "A": (1.0, 1.0),
    "A_neg_x": (1.0, -1.0),
    "P": (PIXEL_ARCSEC, 1.0),
    "P_neg_x": (PIXEL_ARCSEC, -1.0),
}
U_MAX_MEDIAN_ARCSEC = 1.0  # ASSUMPTION (PLAN.md rule 2)
GRID_STEP = 0.25
ICLV2 = "smacs0723-iclv2"
ICLV2_K_RMS_D034R = 0.321  # D-034 addendum r, results/lens_consistency/smacs_djaref/summary.json
CHUNK = 1 << 20


def fixed_z_images(path: Path = IMAGES) -> Table:
    """The ``z_spec`` systems of the Mahler+22 list, each at its z."""
    t = Table.read(path)
    t = t[np.asarray(t["z_spec"], bool)]
    t["system"] = np.asarray(t["system"], str)
    return t


def fetch(session: requests.Session, name: str, dest: Path) -> str:
    """Stream one record file to ``dest``; return its sha256."""
    h = hashlib.sha256()
    with session.get(f"{RECORD}/files/{name}/content", stream=True, timeout=120) as r:
        r.raise_for_status()
        with dest.open("wb") as f:
            for block in r.iter_content(CHUNK):
                h.update(block)
                f.write(block)
    return h.hexdigest()


def build_model(ax, ay, header, scale: float, sign_x: float, sha: str, name: str):
    """Map values x ``scale`` (arcsec), alpha_x x ``sign_x``, at D_ds / D_s = 1 (readme)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(header)
    model = lensmodel.MapLensModel(
        ax * scale * sign_x,
        ay * scale,
        wcs,
        Z_LENS,
        FlatLambdaCDM(H0=COSMOLOGY[0], Om0=COSMOLOGY[1]),
        source=name,
        sha256=sha,
    )
    if not np.isclose(model.pixel_arcsec, PIXEL_ARCSEC, rtol=1e-4):
        raise ValueError(f'{name}: map pixel {model.pixel_arcsec}" != {PIXEL_ARCSEC}"')
    return model


def count_differences(chow: dict, iclv2: dict) -> list[str]:
    """Reported: systems whose predicted image counts differ between the two models."""
    return [s for s in chow if s in iclv2 and chow[s]["n_predicted"] != iclv2[s]["n_predicted"]]


def run(maps: Path, out: Path = OUT, grid_step: float = GRID_STEP) -> dict:
    maps.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(len(MAPS)) as pool:  # two independent downloads
        futs = {m: pool.submit(fetch, requests.Session(), m, maps / m) for m in MAPS}
        shas = {m: f.result() for m, f in futs.items()}
    images = fixed_z_images()
    (ax, hx), (ay, _) = (sw._read(maps / m) for m in MAPS)
    model_sha = hashlib.sha256((shas[MAPS[0]] + shas[MAPS[1]]).encode()).hexdigest()
    models = {
        r: build_model(ax, ay, hx, s, sx, model_sha, f"smacs0723-chow23 ({r})")
        for r, (s, sx) in READINGS.items()
    }
    summary = {
        "plan": str(out / "PLAN.md"),
        "files": {m: {"sha256": shas[m], "bytes": (maps / m).stat().st_size} for m in MAPS},
        "cosmology": {"H0": COSMOLOGY[0], "Om0": COSMOLOGY[1], "z_lens": Z_LENS},
        "image_list": {
            "source": str(IMAGES),
            "images": len(images),
            "systems": sorted(set(images["system"])),
        },
        "map_header": {
            "shape": list(models["A"].shape),
            "pixel_arcsec": round(models["A"].pixel_arcsec, 6),
            "ctype": [hx.get("CTYPE1"), hx.get("CTYPE2")],
            "bunit": hx.get("BUNIT"),
        },
        "check_C": sp.coverage(models["A"], images),
    }
    if summary["check_C"]["passes"]:
        vals = {r: sw.unit_values(m, images) for r, m in models.items()}
        u = sp.unit_check(vals, U_MAX_MEDIAN_ARCSEC)
        summary["check_U"] = {**u, "per_system": vals}
        if u["passes"]:
            k = sw.k_block(models[u["reading"]], images, grid_step)
            lc._write(k["ip"], out / "imageplane_fixedz.ecsv")
            lc._write(k["ip_wls"], out / "imageplane_fixedz_wls.ecsv")
            ki = sw.k_block(lc.load_model(ICLV2)[0], images, grid_step)
            summary["gate_K"] = k["gate_K"]
            summary["reported"] = {
                "grid_step_arcsec": grid_step,
                "grid_half_width_arcsec": k["grid_half_width_arcsec"],
                "K_wls_not_a_gate": k["K_wls"],
                "per_system": k["per_system"],
                "per_system_wls": k["per_system_wls"],
                "iclv2_same_run": {
                    "gate_K": ki["gate_K"],
                    "K_wls": ki["K_wls"],
                    "per_system": ki["per_system"],
                },
                "predicted_count_differs": count_differences(k["per_system"], ki["per_system"]),
                "iclv2_K_rms_arcsec_d034r": ICLV2_K_RMS_D034R,
            }
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--maps", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--grid-step", type=float, default=GRID_STEP)
    ap.add_argument("--keep", action="store_true", help="keep the fetched maps")
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
