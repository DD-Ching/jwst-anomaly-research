"""MACS J1931 stage 4: MUSE Lyα line-flux reference for continuum-free systems 5-7, 100, 101.

Frozen plan: ``results/lens_consistency/macs1931_muse/PLAN.md``. The cube is read only as ESO SODA
``BAND`` cutouts (full field), streamed to a scratch directory and deleted after use.

    python scripts/m1931_muse.py run --scratch DIR      # PLAN.md
    python scripts/m1931_muse.py run4b --scratch DIR    # PLAN_4b.md (Gaia registration)
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

import m0417_canucsref_posthoc as cp
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.wcs import WCS

OUT = Path("results/lens_consistency/macs1931_muse")
IMAGES = Path("results/lens_consistency/macs1931_inputs/m1931_images.ecsv")
MU_OBS = Path("results/lens_consistency/macs1931_gatek/imageplane_spec.ecsv")
TARGETS = Path("results/lens_consistency/macs1931_stage3d/predicted_macs1931.ecsv")
CUBE_ID = "ADP.2017-03-27T15:26:07.707"
SODA = "https://dataportal.eso.org/dataPortal/soda/sync"
SYSTEMS = (5, 6, 7, 100, 101)
LYA = 1215.67
# Frozen constants (PLAN.md rules 1-8; ASSUMPTIONs)
HALF_BAND, CONT_GAP, PEAK_SEARCH, WIN_HALF = 80.0, 30.0, 20.0, 4
PIX, AP_ARCSEC = 0.2, 0.8
N_RAND, R_RAND, SRC_EXCL, SEED = 200, (2.0, 6.0), 1.0, 20261024
BLEND_ARCSEC, CTRL_SEARCH = 1.6, 1.0
MIN_SIGMA, ABSENT_SIGMA, FLUX_RANGE = 5.0, 10.0, (1 / 3, 3.0)
REF_SNR, GATE_MIN, GATE_FRAC = 5.0, 4, 0.75
TESTABLE = ("recovered", "confused", "absent", "ambiguous")


def lam_nom(z: float) -> float:
    return LYA * (1.0 + z)


def continuum_mask(lam: np.ndarray, lnom: float) -> np.ndarray:
    """Rule 1: continuum slices are those more than CONT_GAP from the nominal line."""
    return np.abs(lam - lnom) > CONT_GAP


def subtract_continuum(cube: np.ndarray, cmask: np.ndarray) -> np.ndarray:
    """Data minus each spaxel's median over the continuum slices."""
    return cube - np.nanmedian(cube[cmask], axis=0)[None]


def line_window(spec: np.ndarray, lam: np.ndarray, lnom: float) -> tuple[int, int, int]:
    """Rule 2: (lo, hi, peak) slice indices; hi exclusive."""
    near = np.where(np.abs(lam - lnom) <= PEAK_SEARCH)[0]
    k = int(near[np.nanargmax(spec[near])])
    return max(k - WIN_HALF, 0), min(k + WIN_HALF + 1, len(lam)), k


def narrow_band(sub: np.ndarray, lo: int, hi: int, dlam: float) -> np.ndarray:
    """Rule 3: continuum-subtracted window sum times the slice width."""
    return dlam * sub[lo:hi].sum(axis=0)


def aperture_spectrum(sub: np.ndarray, x: float, y: float) -> np.ndarray:
    r = AP_ARCSEC / PIX
    return np.array([cp.aperture_sum(np.nan_to_num(s), x, y, r) for s in sub])


def local_noise(nb: np.ndarray, x: float, y: float, ex_x, ex_y, rng) -> float:
    """Rule 4: sigma-clipped std of N_RAND random apertures R_RAND away, avoiding listed images."""
    r = AP_ARCSEC / PIX
    vals, tries = [], 0
    while len(vals) < N_RAND and tries < 50 * N_RAND:
        tries += 1
        a, s = rng.uniform(0, 2 * np.pi), rng.uniform(*R_RAND) / PIX
        rx, ry = x + s * np.cos(a), y + s * np.sin(a)
        if (
            len(ex_x)
            and np.hypot(np.asarray(ex_x) - rx, np.asarray(ex_y) - ry).min() * PIX < SRC_EXCL
        ):
            continue
        v = cp.aperture_sum(nb, rx, ry, r)
        if np.isfinite(v):
            vals.append(v)
    if len(vals) < N_RAND:
        raise RuntimeError(f"only {len(vals)} of {N_RAND} random apertures accepted")
    return float(sigma_clipped_stats(vals)[2])


def max_in_radius(nb: np.ndarray, x: float, y: float, radius_arcsec: float) -> float:
    """Rule 7: max aperture sum over spaxel centres within the search radius (NaN if none)."""
    r, rs = AP_ARCSEC / PIX, radius_arcsec / PIX
    best = float("nan")
    for yy in range(int(np.floor(y - rs)), int(np.ceil(y + rs)) + 1):
        for xx in range(int(np.floor(x - rs)), int(np.ceil(x + rs)) + 1):
            if (xx - x) ** 2 + (yy - y) ** 2 > rs * rs:
                continue
            v = cp.aperture_sum(nb, xx, yy, r)
            if np.isfinite(v) and not v <= best:
                best = v
    return best


def is_blended(sep_other_arcsec) -> bool:
    """Rule 5: within BLEND_ARCSEC of an observed image of another system, overlapping window."""
    s = np.asarray(sep_other_arcsec, float)
    return bool(s.size and (s < BLEND_ARCSEC).any())


def classify(in_cube: bool, blended: bool, f_pred: float, f_meas: float, sigma: float) -> str:
    """Rule 7, classes in order."""
    if not in_cube or not np.isfinite(f_meas):
        return "off_cube"
    if blended:
        return "blended"
    if not f_pred >= MIN_SIGMA * sigma:
        return "undetectable_local"
    r = f_meas / f_pred
    if FLUX_RANGE[0] <= r <= FLUX_RANGE[1]:
        return "recovered"
    if r > FLUX_RANGE[1]:
        return "confused"
    if f_pred >= ABSENT_SIGMA * sigma:
        return "absent"
    return "ambiguous"


def gate_l(classes) -> dict:
    """Rule 8."""
    c = [str(v) for v in classes]
    n = sum(v in TESTABLE for v in c)
    ok = sum(v in ("recovered", "confused") for v in c)
    frac = ok / n if n else float("nan")
    return {
        "n_testable": n,
        "n_ok": ok,
        "frac": frac,
        "passes": bool(n >= GATE_MIN and frac >= GATE_FRAC),
    }


def windows_overlap(a: tuple[float, float], b: tuple[float, float]) -> bool:
    return a[0] < b[1] and b[0] < a[1]


# ----------------------------------------------------------------------------------- I/O and run


def fetch_band(z: float, scratch: Path) -> tuple[Path, dict]:
    lnom = lam_nom(z)
    lo, hi = (lnom - HALF_BAND) * 1e-10, (lnom + HALF_BAND) * 1e-10
    url = f"{SODA}?ID=ivo://eso.org/ID?{CUBE_ID}&BAND={lo:.4e}%20{hi:.4e}"
    p = scratch / f"m1931_muse_z{z:.3f}.fits"
    for attempt in range(4):  # back-off on transient 429 / 5xx / network errors
        try:
            h = hashlib.sha256()
            with urllib.request.urlopen(url, timeout=600) as r, open(p, "wb") as f:
                while b := r.read(1 << 20):
                    h.update(b)
                    f.write(b)
            break
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            p.unlink(missing_ok=True)
            code = getattr(e, "code", None)
            if attempt == 3 or (code is not None and code != 429 and code < 500):
                raise
            time.sleep(2 ** (attempt + 1))
    # The file sha256 differs between identical requests (SODA rewrites header cards), so the
    # DATA HDU's array bytes are hashed too; that hash identifies the cube data.
    with fits.open(p) as hd:
        data_sha = hashlib.sha256(np.ascontiguousarray(hd["DATA"].data).tobytes()).hexdigest()
    meta = {"url": url, "bytes": p.stat().st_size, "sha256": h.hexdigest(), "data_sha256": data_sha}
    return p, meta


def load_band(p: Path):
    with fits.open(p) as hd:
        data = hd["DATA"].data.astype(float)
        hdr = hd["DATA"].header
    w = WCS(hdr).celestial
    lam = hdr["CRVAL3"] + (np.arange(hdr["NAXIS3"]) + 1 - hdr["CRPIX3"]) * hdr["CD3_3"]
    return data, w, lam, float(hdr["CD3_3"])


def shifted_wcs(w: WCS, dra: float, ddec: float) -> WCS:
    """A WCS that maps catalogue (Gaia-frame) positions to the cube pixels where MUSE shows them:
    a source at (ra, dec) appears at (ra + dra / cos dec, dec + ddec) in the cube's own WCS."""
    w2 = w.deepcopy()
    w2.wcs.crval = [
        w.wcs.crval[0] - dra / 3600.0 / np.cos(np.radians(w.wcs.crval[1])),
        w.wcs.crval[1] - ddec / 3600.0,
    ]
    return w2


def in_cube(w: WCS, shape, ra: float, dec: float) -> tuple[bool, float, float]:
    x, y = (float(v) for v in w.world_to_pixel_values(ra, dec))
    return bool(0 <= x < shape[-1] and 0 <= y < shape[-2]), x, y


GAIA = OUT / "gaia_dr3.csv"
GAIA_EPOCH, MUSE_EPOCH, REG_MIN, REG_MAD, STAMP = 2016.0, 2015.5, 4, 0.2, 7  # PLAN_4b rule 0


def registration_ok(dra, ddec) -> dict:
    """PLAN_4b rule 0: per-axis median MUSE - Gaia offset (arcsec) and whether it may be applied."""
    a, b = np.asarray(dra, float), np.asarray(ddec, float)
    mad = [float(np.median(np.abs(v - np.median(v)))) if v.size else np.nan for v in (a, b)]
    return {
        "n": int(a.size),
        "dra": float(np.median(a)) if a.size else np.nan,
        "ddec": float(np.median(b)) if b.size else np.nan,
        "mad": mad,
        "applied": bool(a.size >= REG_MIN and max(mad) <= REG_MAD),
    }


def measure_registration(data: np.ndarray, w: WCS) -> dict:
    from photutils.centroids import centroid_2dg

    img = np.nanmedian(data, axis=0)
    g = Table.read(GAIA, format="ascii.csv")
    dra, ddec, stars = [], [], []
    for r in g:
        if np.ma.is_masked(r["pmra"]) or not np.isfinite(float(r["pmra"])):
            continue
        dec0 = float(r["dec"])
        dt = MUSE_EPOCH - GAIA_EPOCH
        ra = float(r["ra"]) + float(r["pmra"]) * dt / 3.6e6 / np.cos(np.radians(dec0))
        dec = dec0 + float(r["pmdec"]) * dt / 3.6e6
        x, y = (float(v) for v in w.world_to_pixel_values(ra, dec))
        xi, yi = int(round(x)), int(round(y))
        st = img[yi - STAMP : yi + STAMP + 1, xi - STAMP : xi + STAMP + 1]
        if (
            yi - STAMP < 0
            or xi - STAMP < 0
            or st.shape != (2 * STAMP + 1,) * 2
            or not np.isfinite(st).all()
        ):
            continue
        cx, cy = centroid_2dg(st - np.median(st))
        rr, dd = w.pixel_to_world_values(xi - STAMP + cx, yi - STAMP + cy)
        dra.append((float(rr) - ra) * np.cos(np.radians(dec)) * 3600.0)
        ddec.append((float(dd) - dec) * 3600.0)
        stars.append(
            {"source_id": str(r["source_id"]), "dra_arcsec": dra[-1], "ddec_arcsec": ddec[-1]}
        )
    return {**registration_ok(dra, ddec), "stars": stars}


def run(scratch: Path, register: bool = False, out: Path = OUT) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    imgs = Table.read(IMAGES)
    mu = {str(r["image_id"]): abs(float(r["magnification"])) for r in Table.read(MU_OBS)}
    tg = Table.read(TARGETS)
    zsys = {int(s): float(imgs["z"][imgs["system"] == str(s)][0]) for s in SYSTEMS}
    lams = {s: (lam_nom(z) - PEAK_SEARCH, lam_nom(z) + PEAK_SEARCH) for s, z in zsys.items()}
    all_ra, all_dec = np.asarray(imgs["ra"], float), np.asarray(imgs["dec"], float)
    rows, cutouts, windows, gate_rows = [], [], {}, []
    rng = np.random.default_rng(SEED)
    for z in sorted(set(zsys.values())):
        p, meta = fetch_band(z, scratch)
        cutouts.append({"z": z, **meta})
        try:
            data, w, lam, dlam = load_band(p)
        finally:
            p.unlink()  # cloud disk rule: the band is held in memory only
        if register:
            reg = measure_registration(data, w)
            cutouts[-1]["registration"] = reg
            if not reg["applied"]:
                raise RuntimeError(f"PLAN_4b rule 0 fails: {reg}")
            w = shifted_wcs(w, reg["dra"], reg["ddec"])
        sub = subtract_continuum(data, continuum_mask(lam, lam_nom(z)))
        ex = w.world_to_pixel_values(all_ra, all_dec)
        for s in [k for k, v in zsys.items() if v == z]:
            others = [o for o in SYSTEMS if o != s and windows_overlap(lams[s], lams[o])]
            oth = imgs[np.isin(imgs["system"], [str(o) for o in others])]

            def sep_other(ra, dec, oth=oth):
                d = np.hypot(
                    (np.asarray(oth["ra"]) - ra) * np.cos(np.radians(dec)),
                    np.asarray(oth["dec"]) - dec,
                )
                return d * 3600.0

            own = imgs[imgs["system"] == str(s)]
            pos = []
            for r in own:
                ok, x, y = in_cube(w, data.shape, float(r["ra"]), float(r["dec"]))
                pos.append(
                    (
                        str(r["image_id"]),
                        ok,
                        x,
                        y,
                        is_blended(sep_other(float(r["ra"]), float(r["dec"]))),
                    )
                )
            specs = [aperture_spectrum(sub, x, y) for _, ok, x, y, _ in pos if ok]
            near = np.abs(lam - lam_nom(z)) <= PEAK_SEARCH
            spec = np.sum(specs, axis=0) if specs else np.full(lam.shape, np.nan)
            if not np.isfinite(spec[near]).any():
                # no in-cube image with a usable spectrum: no window, hence no reference
                windows[s] = {"lam_nom": lam_nom(z), "reference": None, "ref_snr": float("nan")}
                for iid, ok, x, y, bl in pos:
                    gate_rows.append(
                        _row(
                            s,
                            "control",
                            iid,
                            ok,
                            x,
                            y,
                            bl,
                            mu.get(iid),
                            CTRL_SEARCH,
                            None,
                            None,
                            rng,
                            False,
                            {},
                            None,
                            mu,
                        )
                    )
                for r in tg[tg["system"] == str(s)]:
                    rows.append(
                        _row(
                            s,
                            "target",
                            str(r["image_class"]),
                            False,
                            np.nan,
                            np.nan,
                            False,
                            float(r["abs_mu"]),
                            float(r["search_arcsec"]),
                            None,
                            None,
                            rng,
                            False,
                            {},
                            None,
                            mu,
                            ra=float(r["ra"]),
                            dec=float(r["dec"]),
                        )
                    )
                continue
            lo, hi, pk = line_window(spec, lam, lam_nom(z))
            windows[s] = {
                "lam_nom": lam_nom(z),
                "peak": float(lam[pk]),
                "lo": float(lam[lo]),
                "hi": float(lam[hi - 1]),
            }
            nb = narrow_band(sub, lo, hi, dlam)
            meas = {}
            for iid, ok, x, y, bl in pos:
                f = cp.aperture_sum(nb, x, y, AP_ARCSEC / PIX) if ok else float("nan")
                sg = (
                    local_noise(nb, x, y, ex[0], ex[1], rng)
                    if ok and np.isfinite(f)
                    else float("nan")
                )
                meas[iid] = (ok, x, y, bl, f, sg)
            cand = {
                k: v
                for k, v in meas.items()
                if v[0] and not v[3] and np.isfinite(v[4]) and v[5] > 0
            }
            ref = max(cand, key=lambda k: cand[k][4] / cand[k][5]) if cand else None
            ref_snr = cand[ref][4] / cand[ref][5] if ref else float("nan")
            windows[s].update({"reference": ref, "ref_snr": ref_snr})
            usable = ref is not None and ref_snr >= REF_SNR
            for iid, (ok, x, y, bl, _, _) in meas.items():
                if iid == ref:
                    continue
                gate_rows.append(
                    _row(
                        s,
                        "control",
                        iid,
                        ok,
                        x,
                        y,
                        bl,
                        mu.get(iid),
                        CTRL_SEARCH,
                        nb,
                        ex,
                        rng,
                        usable,
                        cand,
                        ref,
                        mu,
                    )
                )
            for r in tg[tg["system"] == str(s)]:
                ok, x, y = in_cube(w, data.shape, float(r["ra"]), float(r["dec"]))
                bl = is_blended(sep_other(float(r["ra"]), float(r["dec"])))
                rows.append(
                    _row(
                        s,
                        "target",
                        str(r["image_class"]),
                        ok,
                        x,
                        y,
                        bl,
                        float(r["abs_mu"]),
                        float(r["search_arcsec"]),
                        nb,
                        ex,
                        rng,
                        usable,
                        cand,
                        ref,
                        mu,
                        ra=float(r["ra"]),
                        dec=float(r["dec"]),
                    )
                )
    g = gate_l([r["class"] for r in gate_rows])
    for r in rows:
        r["flag"] = bool(g["passes"] and r["class"] == "absent")
    Table(rows=gate_rows).write(out / "controls.ecsv", overwrite=True)
    Table(rows=rows).write(out / "targets.ecsv", overwrite=True)
    summary = {
        "cube": CUBE_ID,
        "cutouts_streamed_deleted": cutouts,
        "windows": windows,
        "gate_l": g,
        "target_classes": {
            c: sum(r["class"] == c for r in rows) for c in sorted({r["class"] for r in rows})
        },
        "flags": sum(r["flag"] for r in rows),
    }
    text = json.dumps(cp._clean(summary), indent=1, allow_nan=False, default=float)
    (out / "summary.json").write_text(text + "\n")
    return summary


def _row(
    s,
    kind,
    iid,
    ok,
    x,
    y,
    bl,
    mu_t,
    search,
    nb,
    ex,
    rng,
    usable,
    cand,
    ref,
    mu,
    ra=np.nan,
    dec=np.nan,
):
    out = {
        "system": s,
        "kind": kind,
        "id": iid,
        "ra": ra,
        "dec": dec,
        "x": x,
        "y": y,
        "blended": bl,
        "abs_mu": np.nan if mu_t is None else mu_t,
        "reference": ref or "",
        "f_pred": np.nan,
        "f_meas": np.nan,
        "sigma_loc": np.nan,
        "pred_over_sigma": np.nan,
        "ratio": np.nan,
    }
    if not usable:
        out["class"] = "no_line_reference"
        return out
    if mu_t is None or ref not in mu:
        raise KeyError(f"no model magnification for {iid} or reference {ref}")
    f_pred = cand[ref][4] * out["abs_mu"] / mu[ref]
    f_meas = max_in_radius(nb, x, y, search) if ok else float("nan")
    sg = local_noise(nb, x, y, ex[0], ex[1], rng) if ok and np.isfinite(f_meas) else float("nan")
    out.update(
        {
            "f_pred": f_pred,
            "f_meas": f_meas,
            "sigma_loc": sg,
            "pred_over_sigma": f_pred / sg,
            "ratio": f_meas / f_pred,
        }
    )
    out["class"] = classify(ok, bl, f_pred, f_meas, sg)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cmd", choices=["run", "run4b"])
    ap.add_argument("--scratch", type=Path, required=True)
    a = ap.parse_args(argv)
    a.scratch.mkdir(parents=True, exist_ok=True)
    if a.cmd == "run4b":
        res = run(a.scratch, register=True, out=OUT / "stage4b")
    else:
        res = run(a.scratch)
    print(json.dumps(res, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
