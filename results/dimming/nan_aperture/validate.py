"""PLAN.md G1: fill accuracy on known on-sky apertures. Range-reads the SCI / ERR stamp of every
finite-flux, S/N >= 5 exposure of the committed NEXUS e5 run (<= 150 by file), masks 1 and 2
random non-centre aperture pixels (seed 0; the two are not each other's reflection), fills with
``fill_nan_aperture`` (max_nan 2) and compares with the unmasked flux.
Usage (repo root): validate.py  -> results/dimming/nan_aperture/g1.json"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, "scripts")
import dimming_exposure_check as dec  # noqa: E402
import dimming_screen as ds  # noqa: E402
import persistence_check as pc  # noqa: E402

EPOCH, N_MAX, MEDIAN_MAX, CALL_MIN = 5, 150, 0.05, 0.95  # PLAN.md G1 (ASSUMPTIONs)
t = Table.read(f"results/dimming/persistence_run/nexus/exposures_e{EPOCH}.ecsv")
t = t[np.isfinite(t["flux"]) & (t["snr"] >= pc.MIN_SNR)]
t.sort(["file", "uid"])
t = t[:N_MAX]
obs = ds.load_field(ds.DEFAULT_CONFIG, "nexus")["epochs"][EPOCH]["obs"]
uri = {f.rsplit("/", 1)[-1]: f for d in dec.cal_dirs(obs) for f in pc.list_cal(d)}
rng = np.random.default_rng(0)


def stamp(row):
    x, y = float(row["x"]), float(row["y"])
    xi, yi = int(round(x)), int(round(y))
    with pc._open(uri[str(row["file"])]) as h:
        ny, nx = h["SCI"].header["NAXIS2"], h["SCI"].header["NAXIS1"]
        y0, y1 = max(0, yi - pc.HALF), min(ny, yi + pc.HALF + 1)
        x0, x1 = max(0, xi - pc.HALF), min(nx, xi + pc.HALF + 1)
        sci = np.asarray(h["SCI"].section[y0:y1, x0:x1], dtype=float)
        err = np.asarray(h["ERR"].section[y0:y1, x0:x1], dtype=float)
    return sci, err, x - x0, y - y0, 2.0 if "LONG" in str(row["detector"]).upper() else 3.0


with ThreadPoolExecutor(16) as ex:
    stamps = list(ex.map(stamp, t))


def fill_bilinear(s, e, x, y, r):
    """POST HOC alternative: value at the reflection through the true (sub-pixel) centre."""
    from scipy.ndimage import map_coordinates

    s2, e2 = s.copy(), e.copy()
    bad = np.argwhere(
        ~np.isfinite(s)
        & (
            np.hypot(*(np.mgrid[0 : s.shape[0], 0 : s.shape[1]] - np.array([y, x])[:, None, None]))
            <= r
        )
    )
    good = np.where(np.isfinite(s), s, np.nan)
    for by, bx in bad:
        c = [[2 * y - by], [2 * x - bx]]
        s2[by, bx] = map_coordinates(good, c, order=1, cval=np.nan)[0]
        e2[by, bx] = map_coordinates(np.where(np.isfinite(e), e, np.nan), c, order=1, cval=np.nan)[
            0
        ]
    return s2, e2


def area_corrected(s, e, x, y, r):
    """POST HOC alternative: sum of finite aperture pixels x n_ap / n_finite."""
    yy, xx = np.mgrid[0 : s.shape[0], 0 : s.shape[1]]
    rr = np.hypot(xx - x, yy - y)
    ap, ann = rr <= r, (rr >= 2 * r + 1) & (rr <= 4 * r) & np.isfinite(s)
    ok = ap & np.isfinite(s)
    bkg = float(np.median(s[ann]))
    return float(np.sum(s[ok] - bkg) * ap.sum() / ok.sum())


rows = []
for row, (sci, err, x, y, r) in zip(t, stamps, strict=True):
    f0, e0 = pc.aperture_flux(sci, err, x, y, r, 2 * r + 1, 4 * r)
    yy, xx = np.mgrid[0 : sci.shape[0], 0 : sci.shape[1]]
    cx, cy = int(round(x)), int(round(y))
    ap = (np.hypot(xx - x, yy - y) <= r) & ~((xx == cx) & (yy == cy)) & np.isfinite(sci)
    cand = list(zip(*np.nonzero(ap), strict=True))
    first = cand[rng.integers(len(cand))]
    others = [c for c in cand if c != first and c != (2 * cy - first[0], 2 * cx - first[1])]
    second = others[rng.integers(len(others))]
    out = {"uid": str(row["uid"]), "file": str(row["file"]), "snr0": float(f0 / e0)}
    for k, holes in ((1, [first]), (2, [first, second])):
        s = sci.copy()
        for hy, hx in holes:
            s[hy, hx] = np.nan
        s2, e2, _ = pc.fill_nan_aperture(s, err, x, y, r, 2)
        f, e = pc.aperture_flux(s2, e2, x, y, r, 2 * r + 1, 4 * r)
        out[f"ratio{k}"] = float(f / f0)
        sb, eb = fill_bilinear(s, err, x, y, r)
        out[f"ratio{k}_bilinear"] = float(
            pc.aperture_flux(sb, eb, x, y, r, 2 * r + 1, 4 * r)[0] / f0
        )
        out[f"ratio{k}_area"] = float(area_corrected(s, err, x, y, r) / f0)
        out[f"same_call{k}"] = bool(
            np.isfinite(f) and (f / e >= pc.MIN_SNR) == (f0 / e0 >= pc.MIN_SNR)
        )
    rows.append(out)

res = {
    "n": len(rows),
    "provenance": "derived",
    "source": f"NEXUS e{EPOCH} _cal stamps (S3 range reads)",
}
for k in (1, 2):
    dev = np.abs(np.array([q[f"ratio{k}"] for q in rows]) - 1)
    same = float(np.mean([q[f"same_call{k}"] for q in rows]))
    res[f"mask{k}"] = {
        "median_abs_dev": round(float(np.nanmedian(dev)), 4),
        "p90_abs_dev": round(float(np.nanpercentile(dev, 90)), 4),
        "same_call_frac": round(same, 3),
        "n_finite": int(np.isfinite(dev).sum()),  # fills refused (NaN) are counted here, not hidden
        "pass": bool(np.nanmedian(dev) <= MEDIAN_MAX and same >= CALL_MIN),  # refused = call lost
    }
    res[f"mask{k}"]["median_signed_dev"] = round(
        float(np.nanmedian([q[f"ratio{k}"] - 1 for q in rows])), 4
    )
res["G1_pass"] = res["mask1"]["pass"] and res["mask2"]["pass"]
res["post_hoc_alternatives"] = {  # not frozen in PLAN.md; inputs to a later plan only
    f"mask{k}_{alt}": {
        "median_abs_dev": round(float(np.nanmedian(np.abs(v - 1))), 4),
        "median_signed_dev": round(float(np.nanmedian(v - 1)), 4),
        "p90_abs_dev": round(float(np.nanpercentile(np.abs(v - 1), 90)), 4),
    }
    for k in (1, 2)
    for alt in ("bilinear", "area")
    for v in [np.array([q[f"ratio{k}_{alt}"] for q in rows])]
}
Path("results/dimming/nan_aperture/g1.json").write_text(json.dumps(res, indent=1) + "\n")
print(json.dumps(res, indent=1))
