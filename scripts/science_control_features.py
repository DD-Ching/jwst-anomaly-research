"""Science (SMACS 0723) vs control (CEERS t021) ranked-galaxy feature distributions (D-022).

Plan: ``results/science_control_features/PLAN.md`` (decision rule frozen there before any fetch).
Reads one pipeline run directory, takes each sample's ranked galaxy stratum (the ``scores``
table of ``<sample>``: quality-gated, stars removed), joins its ``features`` rows and compares
every shared feature overall and in ``ref_mag`` bins with the two-sample KS D and the k-sample
Anderson-Darling p.
Writes ``summary.json`` and ``feature_cdfs.png``.

Provenance: catalog features are ``derived`` (from ``observed`` MAST/DJA catalogs); every statistic
here is ``derived``; the thresholds are ASSUMPTIONs from the plan.
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy import stats

SCIENCE = "smacs0723_nircam"
CONTROL = "ceers_t021"
MAG_BINS = [(-np.inf, 24.0), (24.0, 26.0), (26.0, 28.0)]  # ASSUMPTION (plan): ref_mag bins
MIN_PER_BIN = 30  # ASSUMPTION: fewer sources in either field -> the bin is not judged
D_OK = 0.1  # ASSUMPTION (plan): comparable
D_STRONG = 0.2  # ASSUMPTION (plan): "differ strongly"
BRIGHT_BIN = 0  # index of the bin where H2 (cluster members) may explain a difference
# Post hoc (deviation from the plan, added after the first pass): split SMACS by DJA v7.4 eazy
# z_phot to separate H2 (cluster members) from H3 (processing). ASSUMPTION ranges and radius.
Z_MEMBER = (0.25, 0.55)  # around the cluster z = 0.39
Z_BACKGROUND_MIN = 0.8
ZOUT_MATCH_ARCSEC = 0.2


def verdict(max_d: float) -> str:
    """Frozen plan rule on the largest binned KS D of one feature."""
    if not np.isfinite(max_d):
        return "not judged"
    if max_d <= D_OK:
        return "comparable"
    if max_d > D_STRONG:
        return "differ strongly"
    return "intermediate"


def ranked_features(run_dir: Path, sample: str) -> Table:
    feats = Table.read(run_dir / sample / "features.ecsv")
    scores = Table.read(run_dir / sample / "scores.ecsv")
    keep = np.isin(feats["source_uid"], scores["source_uid"])
    return feats[keep]


def compare(a: np.ndarray, b: np.ndarray) -> dict[str, float | int]:
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    out: dict[str, float | int] = {"n_science": int(a.size), "n_control": int(b.size)}
    if min(a.size, b.size) < MIN_PER_BIN:
        return out | {"ks_d": float("nan"), "ad_p": float("nan")}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # AD p is capped/floored at 0.25/0.001 by scipy
        ad = stats.anderson_ksamp([a, b])
    for name, x in (("science", a), ("control", b)):
        q16, q50, q84 = np.percentile(x, [16, 50, 84])
        out |= {f"{name}_q16": float(q16), f"{name}_median": float(q50), f"{name}_q84": float(q84)}
    return out | {"ks_d": float(stats.ks_2samp(a, b).statistic), "ad_p": float(ad.pvalue)}


def zphot_for(run_dir: Path, sample: str, uids: np.ndarray, zout_path: Path) -> np.ndarray:
    """DJA z_phot (``z_phot``) of each ranked source, NaN beyond ZOUT_MATCH_ARCSEC."""
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    src = Table.read(run_dir / sample / "sources.ecsv")
    row = {u_: i for i, u_ in enumerate(src["source_uid"])}  # source_uid is unique per run
    pos = np.array([row.get(x, -1) for x in uids])
    ra = np.where(pos >= 0, np.asarray(src["ra"], float)[pos], np.nan)
    dec = np.where(pos >= 0, np.asarray(src["dec"], float)[pos], np.nan)
    if not np.all(pos >= 0):  # a ranked uid without a source row gets NaN
        ra, dec = np.nan_to_num(ra), np.nan_to_num(dec)
    z = Table.read(zout_path)
    a = SkyCoord(ra * u.deg, dec * u.deg)
    b = SkyCoord(np.asarray(z["ra"], float) * u.deg, np.asarray(z["dec"], float) * u.deg)
    idx, sep, _ = a.match_to_catalog_sky(b)
    zp = np.asarray(z["z_phot"], float)[idx]
    return np.where((sep.arcsec <= ZOUT_MATCH_ARCSEC) & (pos >= 0), zp, np.nan)


def binned(a, ma, b, mb) -> list[dict]:
    out = []
    for lo, hi in MAG_BINS + [(MAG_BINS[-1][1], np.inf)]:  # last: post hoc ref_mag > 28
        sa, sb = (ma > lo) & (ma <= hi), (mb > lo) & (mb <= hi)
        out.append({"ref_mag": [float(lo), float(hi)]} | compare(a[sa], b[sb]))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out", type=Path, default=Path("results/science_control_features"))
    ap.add_argument("--zout", type=Path, help="DJA SMACS eazy zout (post hoc split)")
    args = ap.parse_args()
    sci, ctl = ranked_features(args.run_dir, SCIENCE), ranked_features(args.run_dir, CONTROL)
    shared = [c for c in sci.colnames if c in ctl.colnames and c != "source_uid"]
    msci, mctl = np.asarray(sci["ref_mag"], float), np.asarray(ctl["ref_mag"], float)
    result: dict = {
        "run_dir": args.run_dir.name,
        "provenance": "derived",
        "n_ranked": {SCIENCE: len(sci), CONTROL: len(ctl)},
        "ref_mag_coverage": {
            name: {
                "nan": int(np.isnan(m).sum()),
                "gt_28": int((m > 28).sum()),
                "in_plan_bins": int(sum(((m > lo) & (m <= hi)).sum() for lo, hi in MAG_BINS)),
            }
            for name, m in ((SCIENCE, msci), (CONTROL, mctl))
        },
        "only_science": [c for c in sci.colnames if c not in ctl.colnames],
        "only_control": [c for c in ctl.colnames if c not in sci.colnames],
        "assumptions": {
            "mag_bins": [list(map(float, b)) for b in MAG_BINS],
            "min_per_bin": MIN_PER_BIN,
            "d_ok": D_OK,
            "d_strong": D_STRONG,
        },
        "features": {},
    }
    for col in shared:
        a, b = np.asarray(sci[col], float), np.asarray(ctl[col], float)
        bins = binned(a, msci, b, mctl)[: len(MAG_BINS)]  # the plan's bins only
        ds = [x["ks_d"] for x in bins if np.isfinite(x["ks_d"])]
        faint = [
            x["ks_d"] for i, x in enumerate(bins) if i != BRIGHT_BIN and np.isfinite(x["ks_d"])
        ]
        max_d = max(ds) if ds else float("nan")
        max_faint = max(faint) if faint else float("nan")
        result["features"][col] = {
            "all": compare(a, b),
            "bins": bins,
            "max_binned_ks_d": max_d,
            "verdict": verdict(max_d),
            # H2 exemption (plan): a difference only in the bright bin does not count
            "verdict_excl_bright": verdict(max_faint),
        }
    # Post hoc: the same comparison with an extra ref_mag > 28 bin (half of each sample is there).
    result["post_hoc_faint_bin"] = {
        col: binned(np.asarray(sci[col], float), msci, np.asarray(ctl[col], float), mctl)[-1]
        for col in shared
    }
    if args.zout:
        zp = zphot_for(args.run_dir, SCIENCE, np.asarray(sci["source_uid"]), args.zout)
        groups = {
            "member": (zp >= Z_MEMBER[0]) & (zp <= Z_MEMBER[1]),
            "background": zp >= Z_BACKGROUND_MIN,
        }
        split: dict = {
            "provenance": "derived (post hoc)",
            "zout_sha256_pinned_in": "SOURCES.md (DJA v7.4 SMACS eazy photo-z)",
            "n_matched": int(np.isfinite(zp).sum()),
            "n": {k: int(v.sum()) for k, v in groups.items()},
            "features": {},
        }
        for col in shared:
            a, b = np.asarray(sci[col], float), np.asarray(ctl[col], float)
            split["features"][col] = {k: binned(a[g], msci[g], b, mctl) for k, g in groups.items()}
        result["post_hoc_zphot_split"] = split
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
    plot(sci, ctl, shared, args.out / "feature_cdfs.png")
    for col, r in result["features"].items():
        print(
            f"{col:28s} all D={r['all'].get('ks_d', float('nan')):.3f} "
            f"binned max={r['max_binned_ks_d']:.3f} -> {r['verdict']} / {r['verdict_excl_bright']}"
        )


def plot(sci: Table, ctl: Table, cols: list[str], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(cols)
    ncol = 4
    fig, axes = plt.subplots(-(-n // ncol), ncol, figsize=(3.2 * ncol, 2.4 * -(-n // ncol)))
    msci, mctl = np.asarray(sci["ref_mag"], float), np.asarray(ctl["ref_mag"], float)
    for ax, col in zip(np.ravel(axes), cols, strict=False):
        for (lo, hi), ls in zip(MAG_BINS, ["-", "--", ":"], strict=True):
            for tab, m, colr in ((sci, msci, "C3"), (ctl, mctl, "C0")):
                x = np.asarray(tab[col], float)[(m > lo) & (m <= hi)]
                x = np.sort(x[np.isfinite(x)])
                if x.size:
                    ax.plot(x, np.arange(1, x.size + 1) / x.size, ls, color=colr, lw=1)
        ax.set_title(col, fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in np.ravel(axes)[n:]:
        ax.axis("off")
    fig.suptitle(
        "SMACS (red) vs CEERS t021 (blue); ref_mag ≤24 solid, 24–26 dashed, 26–28 dotted",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=90)


if __name__ == "__main__":
    main()
