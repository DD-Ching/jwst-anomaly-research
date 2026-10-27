"""Post hoc checks of the MACS J0308 two-epoch grid hits (not in the frozen plan).

1. Coverage: for every hit and band, the number of contributing exposures (bits set in ``CON``)
   and the relative weight (``WHT`` over the band's median positive weight) in each epoch.
2. Nearest bright star: distance and position angle to the nearest grid row whose epoch-2 F410M
   flux is above the 99.5th percentile of measured rows (ASSUMPTION, as in the MACS J1931 post hoc).
   ``classify`` reads each group from 1–2 (star wings, the fixed-offset star image found in this
   run, single-exposure coverage); its thresholds are post hoc ASSUMPTIONs.
3. Known case: distance of each listed SN Ares target (programme 9478) to its nearest grid row,
   and that row's per-band significance and Δm.
4. Detector level (``--sw-groups`` too): each such hit mapped into the six epoch-1 LW ``_cal``
   exposures per band (header WCS; SCI and DQ read by HTTP byte range from MAST): detector, pixel,
   the 0.15″ aperture sum over the local 13 × 13 median and the DQ bits there. The same detector
   pixel in every exposure would mark persistence; one exposure only would mark a cosmic ray.
   Controls (group 100 × g + k): for each such group, the three nearest steady rows 2–20″ away
   (epoch-1 F300M flux 1–3 µJy, |F300M σ| < 1; ASSUMPTION), measured the same way, show whether a
   source is seen in every exposure.
5. Fading compact sources: for the groups given by ``--sw-groups``, the background-subtracted
   0.15″-radius aperture sum in the epoch-1 SW mosaics (F115W / F150W / F200W, same visit), read
   by HTTP byte range from MAST with ``cutouts.make_cutouts``, against the local rms of the stamp.

All numbers are ``derived``.

    python scripts/m0308_grid_posthoc.py --data data/transients/macs0308_grid \\
        --out results/transients/macs0308_grid/posthoc.json --sw-groups 1 3 4
"""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.wcs import WCS

E1 = "jw05293-o005_t001_nircam_clear-"
E2 = "jw06882-o015_t015_nircam_clear-"
BANDS = ("f300m", "f410m")
SW1 = ("f115w", "f150w", "f200w")
ARES = {  # programme 9478 target list (www.stsci.edu/jwst-program-info/download/jwst/pdf/9478/)
    "Ares": (47.2225746, 26.76629),
    "Ares_Image2": (47.2342638, 26.76035),
    "Ares_Image3": (47.2386733, 26.76029),
}
STAR_PCT = 99.5  # ASSUMPTION: bright-star flux percentile
APER_ARCSEC = 0.15  # the grid's aperture radius
STAR_NEAR = 3.0  # ASSUMPTION: arcsec from a bright star where wings / spikes dominate
# post hoc: sky-frame offset of groups 1, 3, 4 from their nearest bright star, measured at the
# epoch-1 roll; an instrument-frame artefact moves with PA_APER (stored as ``pa_aper``)
GHOST_SEP, GHOST_PA = 8.33, 120.1
MAST = "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/"


def n_bits(v: int) -> int:
    return bin(int(v) & 0xFFFFFFFF).count("1")  # CON is int32: bit 31 makes it negative


def n_exposures(con: np.ndarray, y: int, x: int) -> int:
    """Bits set over every CON context plane (32 exposures per plane)."""
    return sum(n_bits(con[k, y, x]) for k in range(con.shape[0]))


def coverage(data: Path, hits: Table) -> dict:
    out: dict = {}
    for b in BANDS:
        for tag, stem in (("e1", E1), ("e2", E2)):
            with fits.open(data / f"{stem}{b}_i2d.fits", memmap=True) as h:
                w = WCS(h["SCI"].header)
                con, wht = h["CON"].data, h["WHT"].data
                rows = wht[:: max(wht.shape[0] // 32, 1)]
                ref = float(np.median(rows[rows > 0]))
                x, y = w.world_to_pixel_values(np.asarray(hits["ra"]), np.asarray(hits["dec"]))
                seen = set()
                for g, xi, yi in zip(hits["group"], x, y, strict=True):
                    if int(g) in seen:  # one row per group: its first member
                        continue
                    seen.add(int(g))
                    xi, yi = int(round(float(xi))), int(round(float(yi)))
                    inside = 0 <= yi < wht.shape[0] and 0 <= xi < wht.shape[1]
                    rec = out.setdefault(str(int(g)), {})
                    rec[f"{b}_{tag}_n_exp"] = n_exposures(con, yi, xi) if inside else 0
                    rec[f"{b}_{tag}_wht_rel"] = (
                        round(float(wht[yi, xi]) / ref, 3) if inside else 0.0
                    )
    return out


def bright_stars(grid: Table, hits: Table) -> dict:
    f = np.asarray(grid["F410M_flux2"], float)
    cut = np.nanpercentile(f[np.isfinite(f)], STAR_PCT)
    stars = grid[np.isfinite(f) & (f >= cut)]
    sc = SkyCoord(stars["ra"], stars["dec"], unit="deg")
    out = {"flux_cut_uJy": round(float(cut), 2), "n_stars": len(stars), "groups": {}}
    for r in hits:
        if str(int(r["group"])) in out["groups"] or len(stars) == 0:
            continue
        d = SkyCoord(r["ra"], r["dec"], unit="deg").separation(sc).arcsec
        d[d <= 0.05] = np.inf  # the hit itself when it is bright
        i = int(np.argmin(d))
        pa = sc[i].position_angle(SkyCoord(r["ra"], r["dec"], unit="deg")).deg
        out["groups"][str(int(r["group"]))] = {
            "sep_arcsec": round(float(d[i]), 2),
            "pa_deg": round(float(pa), 1),
            "star_ra": round(float(stars["ra"][i]), 6),
            "star_dec": round(float(stars["dec"][i]), 6),
            "star_f410m_e1_uJy": round(float(stars["F410M_flux1"][i]), 1),
        }
    return out


def classify(cov: dict, star: dict) -> dict:
    """Post hoc reading of each group: every matching label (thresholds are ASSUMPTIONs).

    The fixed-offset rule is a sky-frame offset valid only at the epoch-1 roll (``pa_aper``).
    """
    out = {}
    for g, c in cov.items():
        s = star["groups"].get(g)
        n_exp = min(v for k, v in c.items() if k.endswith("_n_exp"))
        labels = []
        if n_exp == 0:
            labels.append("off_footprint")
        if s and s["sep_arcsec"] < STAR_NEAR:
            labels.append("star_wings")
        if s and abs(s["sep_arcsec"] - GHOST_SEP) < 0.2 and abs(s["pa_deg"] - GHOST_PA) < 1.0:
            labels.append("fixed_offset_star_image")
        if n_exp == 1:
            labels.append("single_exposure")
        out[g] = labels or ["unexplained"]
    return out


def known_case(grid: Table) -> dict:
    c = SkyCoord(grid["ra"], grid["dec"], unit="deg")
    out = {}
    for name, (ra, dec) in ARES.items():
        s = c.separation(SkyCoord(ra, dec, unit="deg")).arcsec
        i = int(np.argmin(s))
        r = grid[i]
        out[name] = {
            "pos_id": str(r["pos_id"]),
            "origin": int(r["origin"]),
            "sep_arcsec": round(float(s[i]), 3),
        }
        for b in ("F300M", "F410M"):
            for k in ("flux1", "flux2", "dmag", "sigma"):
                v = float(r[f"{b}_{k}"])
                out[name][f"{b}_{k}"] = None if not np.isfinite(v) else round(v, 3)
    return out


def sw_check(hits: Table, groups: list[int]) -> dict:
    from jwst_anomaly.cutouts import make_cutouts

    sel = hits[np.isin(hits["group"], groups)]
    targets = Table(
        {"source_uid": [f"g{int(g)}" for g in sel["group"]], "ra": sel["ra"], "dec": sel["dec"]}
    )
    out: dict = {}
    with tempfile.TemporaryDirectory() as tmp:
        for b in SW1:
            cut = make_cutouts(
                f"mast:JWST/product/{E1}{b}_i2d.fits", targets, 2.0, Path(tmp), band=b
            )
            for row in cut:
                rec = out.setdefault(row["source_uid"], {})
                if not row["path"]:
                    rec[b] = None
                    continue
                with fits.open(row["path"]) as h:
                    img = h["SCI"].data.astype(float)
                    scale = abs(
                        WCS(h["SCI"].header).proj_plane_pixel_scales()[0].to_value("arcsec")
                    )
                yy, xx = np.indices(img.shape)
                cy, cx = (img.shape[0] - 1) / 2, (img.shape[1] - 1) / 2
                rr = np.hypot(yy - cy, xx - cx) * scale
                ann = img[(rr > 0.4) & (rr < 0.8) & np.isfinite(img)]
                bkg, rms = (
                    float(np.median(ann)),
                    float(1.4826 * np.median(np.abs(ann - np.median(ann)))),
                )
                ap = (rr <= APER_ARCSEC) & np.isfinite(img)
                s = float(np.sum(img[ap] - bkg))
                rec[b] = {
                    "aper_sum": round(s, 4),
                    "snr_pixel_rms": round(s / (rms * np.sqrt(ap.sum())), 2),
                }
    return out


CAL1 = {
    b: [f"jw05293005001_031{v}_0000{d}_nrc{m}long_cal.fits" for d in (1, 2, 3) for m in ("a", "b")]
    for b, v in (("f300m", "03"), ("f410m", "01"))
}


def detector_level(hits: Table, groups: list[int]) -> list[dict]:
    import fsspec

    sel = hits[np.isin(hits["group"], groups)]
    out = []
    for b, names in CAL1.items():
        for name in names:
            with fsspec.open(MAST + name, "rb", block_size=2**16, cache_type="readahead") as fo:
                with fits.open(fo, lazy_load_hdus=True) as h:
                    w = WCS(h["SCI"].header)
                    ny, nx = h["SCI"].header["NAXIS2"], h["SCI"].header["NAXIS1"]
                    for r in sel:
                        x, y = w.world_to_pixel_values(float(r["ra"]), float(r["dec"]))
                        xi, yi = int(round(float(x))), int(round(float(y)))
                        if not (6 <= xi < nx - 6 and 6 <= yi < ny - 6):
                            continue
                        sci = h["SCI"].section[yi - 6 : yi + 7, xi - 6 : xi + 7].astype(float)
                        dq = h["DQ"].section[yi - 2 : yi + 3, xi - 2 : xi + 3]
                        yy, xx = np.indices(sci.shape)
                        ap = np.hypot(yy - 6, xx - 6) <= APER_ARCSEC / 0.063 + 0.5
                        bkg = float(np.nanmedian(sci))
                        out.append(
                            {
                                "group": int(r["group"]),
                                "band": b,
                                "cal": name,
                                "x": xi,
                                "y": yi,
                                "aper_sum_over_median": round(float(np.nansum(sci[ap] - bkg)), 3),
                                "dq_or": int(np.bitwise_or.reduce(dq.ravel())),
                            }
                        )
    return out


def controls(grid: Table, hits: Table, groups: list[int]) -> Table:
    c = SkyCoord(grid["ra"], grid["dec"], unit="deg")
    f1, sig = np.asarray(grid["F300M_flux1"], float), np.asarray(grid["F300M_sigma"], float)
    rows = []
    for g in groups:
        m = hits[hits["group"] == g]
        if len(m) == 0:
            raise SystemExit(f"--sw-groups: group {g} not in hits")
        r = m[0]
        s = c.separation(SkyCoord(r["ra"], r["dec"], unit="deg")).arcsec
        ok = (s > 2) & (s < 20) & (np.abs(sig) < 1) & (f1 > 1) & (f1 < 3)
        idx = np.where(ok)[0][np.argsort(s[ok])][:3]
        rows += [
            (100 * g + k, float(grid["ra"][i]), float(grid["dec"][i])) for k, i in enumerate(idx)
        ]
    return Table(rows=rows, names=("group", "ra", "dec"), dtype=(int, float, float))


def pa_aper(data: Path) -> dict:
    out = {}
    for tag, stem in (("e1", E1), ("e2", E2)):
        hdr = fits.getheader(data / f"{stem}f300m_i2d.fits", "SCI")
        out[tag] = round(float(hdr["PA_APER"]), 2)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--sw-groups", type=int, nargs="*", default=[])
    a = ap.parse_args(argv)
    hits = Table.read(a.data / "out" / "hits.ecsv")
    grid = Table.read(a.data / "out" / "grid.ecsv")
    cov, star = coverage(a.data, hits), bright_stars(grid, hits)
    det: dict = {"hits": [], "controls": []}
    if a.sw_groups:  # one pass over the remote _cal files for hits and controls
        ctl = controls(grid, hits, a.sw_groups)
        both = vstack([hits[np.isin(hits["group"], a.sw_groups)]["group", "ra", "dec"], ctl])
        rows = detector_level(both, list(both["group"]))
        det["hits"] = [r for r in rows if r["group"] in a.sw_groups]
        det["controls"] = [r for r in rows if r["group"] not in a.sw_groups]
        if not det["controls"]:
            print("warning: no control rows measured")
    res = {
        "provenance": "derived",
        "coverage": cov,
        "bright_star": star,
        "classification": classify(cov, star),
        "known_case": known_case(grid),
        "sw_epoch1": sw_check(hits, a.sw_groups) if a.sw_groups else {},
        "detector_level": det["hits"],
        "detector_level_controls": det["controls"],
        "pa_aper": pa_aper(a.data),
        "mosaic_bytes": {p.name: p.stat().st_size for p in sorted(a.data.glob("*_i2d.fits"))},
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: res[k] for k in ("classification", "sw_epoch1")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
