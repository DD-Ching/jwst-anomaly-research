"""DJA v7.4 0.5-arcsec colours (CEERS control ``dja05``, D-022) vs the PSF-matched CEERS DR1.0.

Plan: ``results/ceers_dr1_colours/PLAN.md``. Reads the two catalogs from local paths (streamed into
scratch and deleted after use; pinned by sha256 in the result), matches one-to-one within 0.2 arcsec
and compares five adjacent colours. Writes ``summary.json`` and ``colour_offsets.png``.

Provenance: both inputs are ``observed`` catalog values; every difference and statistic is
``derived``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

from jwst_anomaly.photometry import BAD_FLAGS

# DJA ``{band}_tot_corr`` is one aperture-to-total factor per source, the same in every band, so
# it cannot change a colour (checked on this run); it is not used.
BANDS = ["f115w", "f150w", "f200w", "f277w", "f356w", "f444w"]
COLOURS = list(zip(BANDS[:-1], BANDS[1:], strict=True))
MATCH_ARCSEC = 0.2  # ASSUMPTION: the D-022 matching radius
MIN_SNR = 10.0  # ASSUMPTION: per band, in both catalogs
MAX_ABS_MEDIAN = 0.05  # ASSUMPTION (plan): |median colour offset| in mag
MAX_NORM_NMAD = 1.5  # ASSUMPTION (plan): NMAD of the error-normalised offset
# 0-based index of 0.5" in the readme's aperture list (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, ...)
CEERS_APER_05 = 5
# ASSUMPTION: DJA flux_radius (0.04"/pix detection image) splitting compact / extended
COMPACT_R50_PIX = 2.5
# t021 footprint (approximate box around jw01345-o001_t021; ASSUMPTION, a subsample check only)
T021_BOX = (214.80, 214.98, 52.78, 52.94)


def nmad(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if x.size else float("nan")


def load(dja_path: Path, ceers_path: Path) -> dict[str, np.ndarray]:
    d = fits.getdata(dja_path, 1)
    c = fits.getdata(ceers_path, 1)
    keep_c = (c["STAR_FLAG"] == 0) & (c["BAD_REGION_FLAG"] == 0)
    c = c[keep_c]
    sd = SkyCoord(d["ra"] * u.deg, d["dec"] * u.deg)
    sc = SkyCoord(c["RA"] * u.deg, c["DEC"] * u.deg)
    i, sep, _ = sc.match_to_catalog_sky(sd)
    j, _, _ = sd[i].match_to_catalog_sky(sc)
    ok = (sep.arcsec < MATCH_ARCSEC) & (j == np.arange(len(sc)))
    d, c, sep = d[i[ok]], c[ok], sep.arcsec[ok]
    out = {"ra": np.asarray(c["RA"], float), "dec": np.asarray(c["DEC"], float), "sep": sep}
    out["r50"] = np.asarray(d["flux_radius"], float)
    for b in BANDS:
        B = b.upper()
        fd = np.asarray(d[f"{b}_flux_aper_1"], float)
        ed = np.asarray(d[f"{b}_fluxerr_aper_1"], float)
        bad = (np.asarray(d[f"{b}_flag_aper_1"], int) & BAD_FLAGS) != 0
        fd[bad] = np.nan
        fc = np.asarray(c[f"{B}_FLUX_APER"][:, CEERS_APER_05], float)
        ec = np.asarray(c[f"{B}_FLUXERR_APER_EMP"][:, CEERS_APER_05], float)
        out[f"dja_{b}"], out[f"dja_{b}_err"] = fd, ed
        out[f"ceers_{b}"], out[f"ceers_{b}_err"] = fc, ec
        out[f"ceers_auto_{b}"] = np.asarray(c[f"{B}_FLUX"], float)
    return out


def colour(f1, f2):
    with np.errstate(divide="ignore", invalid="ignore"):
        return -2.5 * np.log10(f1 / f2)


def colour_err(f1, e1, f2, e2):
    with np.errstate(divide="ignore", invalid="ignore"):
        return 1.0857 * np.hypot(e1 / f1, e2 / f2)


def compare(m: dict[str, np.ndarray], sel: np.ndarray) -> dict[str, dict]:
    res = {}
    for b1, b2 in COLOURS:
        good = sel.copy()
        for cat in ("dja", "ceers"):
            for b in (b1, b2):
                f, e = m[f"{cat}_{b}"], m[f"{cat}_{b}_err"]
                good &= np.isfinite(f) & np.isfinite(e) & (e > 0) & (f / e >= MIN_SNR)
        cd = colour(m[f"dja_{b1}"][good], m[f"dja_{b2}"][good])
        cc = colour(m[f"ceers_{b1}"][good], m[f"ceers_{b2}"][good])
        cauto = colour(m[f"ceers_auto_{b1}"][good], m[f"ceers_auto_{b2}"][good])
        sig = np.hypot(
            colour_err(m[f"dja_{b1}"], m[f"dja_{b1}_err"], m[f"dja_{b2}"], m[f"dja_{b2}_err"])[
                good
            ],
            colour_err(
                m[f"ceers_{b1}"], m[f"ceers_{b1}_err"], m[f"ceers_{b2}"], m[f"ceers_{b2}_err"]
            )[good],
        )
        delta = cd - cc
        r = {
            "n": int(good.sum()),
            "median_delta": float(np.median(delta)) if good.any() else None,
            "nmad_delta": nmad(delta),
            "nmad_norm": nmad(delta / sig),
            "median_sigma": float(np.median(sig)) if good.any() else None,
            "median_delta_ceers_aper05_minus_auto": float(np.nanmedian(cc - cauto))
            if good.any()
            else None,
        }
        r["passes"] = bool(
            r["n"] > 0
            and abs(r["median_delta"]) <= MAX_ABS_MEDIAN
            and r["nmad_norm"] <= MAX_NORM_NMAD
        )
        res[f"{b1}-{b2}"] = r
    return res


MAG_BINS = ((18.0, 23.0), (23.0, 25.0))  # CEERS 0.5" F200W AB bins for the per-band table
R50_BINS = ((0.0, 2.5), (2.5, 4.0), (4.0, 99.0))  # DJA flux_radius bins [pix]


def per_band(m: dict[str, np.ndarray]) -> list[dict]:
    """Median DJA - CEERS 0.5" magnitude per band in F200W-magnitude and size bins."""
    with np.errstate(divide="ignore", invalid="ignore"):
        mag = 31.4 - 2.5 * np.log10(m["ceers_f200w"])  # nJy -> AB
    rows = []
    for lo, hi in MAG_BINS:
        for r0, r1 in R50_BINS:
            s = (mag >= lo) & (mag < hi) & (m["r50"] >= r0) & (m["r50"] < r1)
            row = {"f200w_mag": [lo, hi], "r50_pix": [r0, r1], "n": int(s.sum())}
            for b in BANDS:
                fd, fc = m[f"dja_{b}"][s] * 1000.0, m[f"ceers_{b}"][s]  # uJy -> nJy
                g = np.isfinite(fd) & np.isfinite(fc) & (fd > 0) & (fc > 0)
                row[b] = float(np.median(-2.5 * np.log10(fd[g] / fc[g]))) if g.any() else None
            rows.append(row)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dja", type=Path, required=True)
    ap.add_argument("--ceers", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("results/ceers_dr1_colours"))
    a = ap.parse_args()
    m = load(a.dja, a.ceers)
    n = len(m["ra"])
    allsel = np.ones(n, bool)
    compact = m["r50"] < COMPACT_R50_PIX
    ra0, ra1, de0, de1 = T021_BOX
    t021 = (m["ra"] > ra0) & (m["ra"] < ra1) & (m["dec"] > de0) & (m["dec"] < de1)
    with np.errstate(divide="ignore", invalid="ignore"):
        f200w_mag = 31.4 - 2.5 * np.log10(m["ceers_f200w"])  # nJy -> AB
    summary = {
        "provenance": "derived",
        "inputs": {"dja": a.dja.name, "ceers": a.ceers.name},
        "n_matched": int(n),
        "median_sep_arcsec": float(np.median(m["sep"])),
        "thresholds_assumption": {
            "match_arcsec": MATCH_ARCSEC,
            "min_snr": MIN_SNR,
            "max_abs_median_mag": MAX_ABS_MEDIAN,
            "max_norm_nmad": MAX_NORM_NMAD,
            "compact_r50_pix": COMPACT_R50_PIX,
            "t021_box_deg": T021_BOX,
        },
        "all": compare(m, allsel),
        "compact": compare(m, compact),
        "extended": compare(m, ~compact),
        "t021_box": compare(m, t021),
        "bright_f200w_lt_25": compare(m, f200w_mag < 25.0),
        "faint_f200w_ge_25": compare(m, f200w_mag >= 25.0),
        "per_band_dja_minus_ceers_mag": per_band(m),
    }
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    plot(m, a.out / "colour_offsets.png")
    for k in ("all", "compact", "extended", "t021_box", "bright_f200w_lt_25", "faint_f200w_ge_25"):
        print(k)
        for col, r in summary[k].items():
            if r["n"] == 0:
                print(f"  {col:12s} n=     0 (empty selection)")
                continue
            print(
                f"  {col:12s} n={r['n']:6d} med={r['median_delta']:+.3f} nmad={r['nmad_delta']:.3f}"
                f" norm={r['nmad_norm']:.2f} sig={r['median_sigma']:.3f}"
                f" ceers05-auto={r['median_delta_ceers_aper05_minus_auto']:+.3f}"
                f" pass={r['passes']}"
            )
    for row in summary["per_band_dja_minus_ceers_mag"]:
        print(row)


def plot(m: dict[str, np.ndarray], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(COLOURS), figsize=(18, 3.6), sharey=True)
    for ax, (b1, b2) in zip(axes, COLOURS, strict=True):
        good = np.ones(len(m["ra"]), bool)
        for cat in ("dja", "ceers"):
            for b in (b1, b2):
                f, e = m[f"{cat}_{b}"], m[f"{cat}_{b}_err"]
                good &= np.isfinite(f) & (e > 0) & (f / e >= MIN_SNR)
        delta = colour(m[f"dja_{b1}"], m[f"dja_{b2}"]) - colour(m[f"ceers_{b1}"], m[f"ceers_{b2}"])
        ax.scatter(m["r50"][good], delta[good], s=1, alpha=0.2, color="k", rasterized=True)
        ax.axhline(0, color="tab:red", lw=0.8)
        ax.set_xscale("log")
        ax.set_xlim(1, 20)
        ax.set_ylim(-0.8, 0.8)
        ax.set_title(f"{b1.upper()}-{b2.upper()}")
        ax.set_xlabel("DJA flux_radius [pix]")
    axes[0].set_ylabel('colour DJA 0.5" - CEERS DR1.0 0.5" [mag]')
    fig.tight_layout()
    fig.savefig(path, dpi=110)


if __name__ == "__main__":
    main()
