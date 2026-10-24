"""RXCJ2211-0350: input check for the fitted-model route (no gate-K residual or prediction read).

The VENUS GLAFIC v0.91 model (Zenodo 10.5281/zenodo.20312327, CC0; paper arXiv:2512.05180 per the
record) releases kappa, gamma1 and gamma2 maps at z_s = 3 and ``full_output_v91.csv`` (images, their
deflections and back-traced sources), but no deflection maps. Checks frozen in
``results/lens_consistency/rxcj2211_inputs/PLAN.md``:

- S: the CSV source equals the image minus its deflection (sign convention read from the data);
- D: deflections from the kappa map (FFT, zero-padded) reproduce the CSV deflections;
- D-affine: only if D fails, D with an affine term fitted on odd systems, tested on even ones;
- R: a DJA v7 imaging field for the cluster.

Writes ``summary.json`` and ``rj2211_images.ecsv``. The zip is held in memory only.

    python scripts/rxcj2211_inputs.py
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import macs1931_inputs as m1931
import numpy as np
import requests
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.table import Table
from scipy.ndimage import map_coordinates
from scipy.signal import fftconvolve

OUT = Path("results/lens_consistency/rxcj2211_inputs")
RECORD = "https://zenodo.org/api/records/20312327"
ZIP = "RXCJ2211_lens_model_GLAFIC.zip"
ZIP_URL = f"https://zenodo.org/records/20312327/files/{ZIP}?download=1"
DIR = "RXCJ2211_lens_model_GLAFIC/"
CSV = DIR + "full_output_v91.csv"
README = DIR + "glafic_rj2211_v0.91_readme.txt"
KAPPA = DIR + "hlsp_venus_model_rj2211_glafic_v0.91_kappa.fits"
DJA_NAMES = ("rxj2211", "rxcj2211", "j2211", "rj2211")
Z_LENS = 0.397  # README zl = 0.3970
COSMO = FlatLambdaCDM(H0=70.0, Om0=0.3)  # README: omega 0.3, lambda 0.7, hubble 0.7
S_MAX = 0.01  # arcsec, check S (PLAN)
D_MEDIAN_MAX = 0.10  # arcsec, check D (ASSUMPTION, PLAN)
D_P90_MAX = 0.25  # arcsec, check D (ASSUMPTION, PLAN)


def offsets_arcsec(ra, dec, ra0, dec0):
    """Tangent-plane (east, north) offsets in arcsec of (ra, dec) from (ra0, dec0); small fields."""
    cosd = np.cos(np.radians(dec0))
    return (np.asarray(ra) - ra0) * cosd * 3600.0, (np.asarray(dec) - dec0) * 3600.0


def check_s(t: Table) -> dict:
    """Image minus deflection against the CSV source, for both readings of the x sign."""
    ok = np.isfinite(np.asarray(t["deflect_x_best"], float))
    r = t[ok]
    ex, ny = offsets_arcsec(r["ra"], r["dec"], r["ra_source_best"], r["dec_source_best"])
    out = {"rows": int(ok.sum())}
    for name, sx in (("x_along_minus_ra", -1.0), ("x_along_plus_ra", 1.0)):
        # image - source = deflection; east offset = sx * deflect_x
        d = np.hypot(ex - sx * r["deflect_x_best"], ny - r["deflect_y_best"])
        out[name] = round(float(np.median(d)), 5)
    best = min(("x_along_minus_ra", "x_along_plus_ra"), key=lambda k: out[k])
    out.update({"convention": best, "passes": out[best] <= S_MAX})
    return out


def alpha_from_kappa(kappa: np.ndarray, pix: float) -> tuple[np.ndarray, np.ndarray]:
    """Deflection (arcsec, along the array's +x / +y) of a convergence map on a square grid of
    ``pix`` arcsec, alpha = (1/pi) * int kappa(t') (t - t') / |t - t'|^2 d^2t', with the kernel on a
    grid twice the map size (zero padding; mass outside the map is absent)."""
    ny, nx = kappa.shape
    yy, xx = np.mgrid[-ny + 1 : ny, -nx + 1 : nx].astype(float) * pix
    r2 = xx**2 + yy**2
    r2[ny - 1, nx - 1] = np.inf  # the self term carries no deflection
    kx, ky = xx / (np.pi * r2), yy / (np.pi * r2)
    area = pix * pix
    ax = fftconvolve(kappa, kx, mode="same") * area
    ay = fftconvolve(kappa, ky, mode="same") * area
    return ax, ay


def dls_ds(z) -> np.ndarray:
    z = np.atleast_1d(np.asarray(z, float))
    return (
        np.array([COSMO.angular_diameter_distance(Z_LENS, zi).value for zi in z])
        / COSMO.angular_diameter_distance(z).value
    )


def sample(ax, ay, wcs_hdr, ra, dec, z, zs_map):
    """Map deflections at (ra, dec), rescaled from the map's z_s to each image z, as (east, north)
    arcsec. The map is RA---CAR / DEC--CAR with CDELT1 < 0: array +x is −RA (west)."""
    cd1, cd2 = float(wcs_hdr["CDELT1"]), float(wcs_hdr["CDELT2"])
    e, n = offsets_arcsec(ra, dec, wcs_hdr["CRVAL1"], wcs_hdr["CRVAL2"])
    # FITS pixels are 1-based; numpy 0-based
    px = wcs_hdr["CRPIX1"] - 1 + (e / 3600.0) / cd1
    py = wcs_hdr["CRPIX2"] - 1 + (n / 3600.0) / cd2
    vx = map_coordinates(ax, [py, px], order=1, mode="nearest")
    vy = map_coordinates(ay, [py, px], order=1, mode="nearest")
    scale = dls_ds(z) / dls_ds(zs_map)[0]
    return np.sign(cd1) * vx * scale, vy * scale  # array +x maps to sign(CDELT1) east


def compare(pred_e, pred_n, csv_e, csv_n) -> dict:
    d = np.hypot(pred_e - csv_e, pred_n - csv_n)
    med, p90 = float(np.median(d)), float(np.percentile(d, 90))
    return {
        "n": int(d.size),
        "median_arcsec": round(med, 4),
        "p90_arcsec": round(p90, 4),
        "max_arcsec": round(float(d.max()), 4),
        "passes": bool(med <= D_MEDIAN_MAX and p90 <= D_P90_MAX),
    }


def affine_fit(theta_e, theta_n, res_e, res_n):
    """Least-squares a + B theta for the residual deflection (6 parameters)."""
    A = np.column_stack([np.ones_like(theta_e), theta_e, theta_n])
    ce, *_ = np.linalg.lstsq(A, res_e, rcond=None)
    cn, *_ = np.linalg.lstsq(A, res_n, rcond=None)
    return ce, cn


def affine_apply(ce, cn, theta_e, theta_n):
    A = np.column_stack([np.ones_like(theta_e), theta_e, theta_n])
    return A @ ce, A @ cn


def dja_check(index_html: str) -> dict:
    """Fields of the DJA v7 imaging index whose name holds one of ``DJA_NAMES``."""
    fields = sorted(set(m1931.re.findall(r"([a-z0-9][a-z0-9_-]*?)-grizli-v7", index_html.lower())))
    if not fields:
        raise ValueError("no field parsed from the DJA index")
    hits = [f for f in fields if any(n in f for n in DJA_NAMES)]
    return {
        "index": m1931.DJA_INDEX,
        "fields_listed": len(fields),
        "names_searched": list(DJA_NAMES),
        "matches": hits,
        "passes": bool(hits),
    }


def read_csv(text: str) -> Table:
    t = Table.read(text, format="ascii.csv", fill_values=[("", "nan")])
    for c in t.colnames:
        if c != "id" and hasattr(t[c], "filled"):
            t[c] = t[c].filled(np.nan)
    t["id"] = [f"{float(v):.1f}" if not isinstance(v, str) else v for v in t["id"]]
    t["system"] = [s.split(".")[0] for s in t["id"]]
    return t


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.parse_args(argv)
    session = requests.Session()
    record = m1931._get(session, RECORD).json()
    listed = {f["key"]: f for f in record["files"]}
    blob = m1931._get(session, ZIP_URL).content
    if "md5:" + hashlib.md5(blob).hexdigest() != listed[ZIP]["checksum"]:
        raise RuntimeError(f"{ZIP}: md5 differs from the record")
    zf = zipfile.ZipFile(io.BytesIO(blob))
    readme = zf.read(README).decode()
    t = read_csv(zf.read(CSV).decode())
    with fits.open(io.BytesIO(zf.read(KAPPA))) as hdul:
        hdr = hdul[0].header
        kappa = np.asarray(hdul[0].data, float)
    pix = abs(float(hdr["CDELT2"])) * 3600.0
    zs_map = float(hdr["ZS"])

    s = check_s(t)
    ok = np.isfinite(np.asarray(t["deflect_x_best"], float)) & (np.asarray(t["z_lensing"]) > 0)
    r = t[ok]
    sx = -1.0 if s["convention"] == "x_along_minus_ra" else 1.0
    csv_e, csv_n = sx * np.asarray(r["deflect_x_best"]), np.asarray(r["deflect_y_best"])
    ax, ay = alpha_from_kappa(kappa, pix)
    pe, pn = sample(ax, ay, hdr, r["ra"], r["dec"], r["z_lensing"], zs_map)
    d = compare(pe, pn, csv_e, csv_n)
    d_aff = None
    if not d["passes"]:
        te, tn = offsets_arcsec(r["ra"], r["dec"], hdr["CRVAL1"], hdr["CRVAL2"])
        odd = np.array([int(x) % 2 == 1 for x in r["system"]])
        ce, cn = affine_fit(te[odd], tn[odd], (csv_e - pe)[odd], (csv_n - pn)[odd])
        fe, fn = affine_apply(ce, cn, te, tn)
        d_aff = {
            "fit_on": "odd systems",
            "fit_rows": int(odd.sum()),
            "coef_east": [round(float(c), 6) for c in ce],
            "coef_north": [round(float(c), 6) for c in cn],
            "test_even": compare((pe + fe)[~odd], (pn + fn)[~odd], csv_e[~odd], csv_n[~odd]),
            "train_odd_reported": compare((pe + fe)[odd], (pn + fn)[odd], csv_e[odd], csv_n[odd]),
        }
        d_aff["passes"] = d_aff["test_even"]["passes"]
    # post hoc (added after the frozen checks failed; no verdict depends on it): one scale factor
    s_fit = float(np.sum(pe * csv_e + pn * csv_n) / np.sum(pe**2 + pn**2))
    posthoc = {
        "label": "post hoc diagnostic, added after D and D-affine failed",
        "scale_fit_all_rows": round(s_fit, 4),
        "after_scale": compare(s_fit * pe, s_fit * pn, csv_e, csv_n),
        "dls_ds_zs3": round(float(dls_ds(zs_map)[0]), 4),
        "median_csv_deflection_arcsec": round(float(np.median(np.hypot(csv_e, csv_n))), 3),
        "median_map_deflection_arcsec": round(float(np.median(np.hypot(pe, pn))), 3),
    }
    posthoc["after_scale"].pop("passes")
    dja = dja_check(m1931._get(session, m1931.DJA_INDEX).text)

    zs = np.asarray(t["z_lensing"], float)
    spec = sorted(
        {ln.split()[2] for ln in readme.splitlines() if ln.startswith("# ZS_REF")}, key=int
    )
    summary = {
        "sources": {
            "record": f"{RECORD} (doi {record.get('doi')}, published "
            f"{record['metadata'].get('publication_date')}, license "
            f"{record['metadata'].get('license', {}).get('id')})",
            "zip_md5": listed[ZIP]["checksum"],
            "zip_bytes": len(blob),
            "members_sha256": {
                n: hashlib.sha256(zf.read(n)).hexdigest() for n in (CSV, README, KAPPA)
            },
        },
        "readme": {
            ln[2:].split("=")[0].strip(): ln.split("=", 1)[1].strip()
            for ln in readme.splitlines()[:3]
            if "=" in ln and "," not in ln
        }
        | {"fit": readme.splitlines()[2][2:].strip()},
        "image_list": {
            "rows": len(t),
            "systems": len(set(t["system"])),
            "rows_z_free": int((zs < 0).sum()),
            "systems_z_free": sorted(
                {s for s, z in zip(t["system"], zs, strict=True) if z < 0}, key=int
            ),
            "systems_with_zs_ref": spec,
            "rows_with_deflection": int(np.isfinite(np.asarray(t["deflect_x_best"], float)).sum()),
        },
        "maps": {
            "released": sorted(n.removeprefix(DIR) for n in zf.namelist() if n.endswith(".fits")),
            "deflection_maps_released": False,
            "shape": list(kappa.shape),
            "pixel_arcsec": round(pix, 6),
            "ctype": [hdr["CTYPE1"], hdr["CTYPE2"]],
            "zs": zs_map,
            "kappa_edge_median": round(
                float(np.median(np.r_[kappa[0], kappa[-1], kappa[:, 0], kappa[:, -1]])), 4
            ),
        },
        "check_S": s,
        "check_D": d,
        "check_D_affine": d_aff,
        "reference_catalogue_dja_v7": dja,
        "posthoc": posthoc,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    keep = ["id", "system", "ra", "dec", "z_lensing", "mu_best", "mu_16", "mu_84"]
    img = t[keep]
    img.meta = {
        "provenance": "observed (positions); z_lensing: spec-z for systems with a ZS_REF in the "
        "README, else the model's (photo-z / fitted); -1 = free; mu: model_prediction",
        "source": f"{ZIP_URL} member {CSV}",
    }
    img.write(OUT / "rj2211_images.ecsv", format="ascii.ecsv", overwrite=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
