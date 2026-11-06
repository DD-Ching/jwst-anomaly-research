"""Science (SMACS 0723) vs control (CEERS t021) ranked-feature distributions (TASKS Then 8).

Plan: ``results/science_control_features/PLAN.md``. Reads one pipeline run directory, takes each
sample's ranked galaxy stratum (the rows of ``<sample>/scores``) with their ``<sample>/features``
values, and compares every feature both samples compute: medians, 16/84 % quantiles, two-sample KS
D and the k-sample Anderson-Darling p, overall and in ``ref_mag`` bins. Writes ``summary.json`` and
``feature_distributions.png``.

Provenance: the feature values are ``derived`` from ``observed`` catalogs; every statistic here is
``derived``; the thresholds are ASSUMPTIONs from the plan.
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy import stats

# ASSUMPTION (plan): ref_mag bins separating depth from population
MAG_BINS = [(-np.inf, 24.0), (24.0, 26.0), (26.0, 28.0)]
D_COMPARABLE = 0.1  # ASSUMPTION (plan): KS D at or below this in every bin -> comparable
D_STRONG = 0.2  # ASSUMPTION (plan): KS D above this in any bin -> "differ strongly"
MIN_N = 20  # ASSUMPTION: fewer finite values per sample in a bin -> statistic not reported
ID_COLUMNS = ("source_uid", "sample_id", "ra", "dec")
# Post hoc (chosen after seeing the 24-26 mag histograms): a compact, F277W-F356W-blue population.
COMPACT_LOG_CI = 0.23  # ASSUMPTION (post hoc): log10 CI_50_30 below this = compact
BLUE_F277W_F356W = -0.3  # ASSUMPTION (post hoc): stellar-like blue F277W-F356W
POSTHOC_BINS = [*MAG_BINS, (28.0, np.inf)]


def ranked_features(run_dir: Path, sample: str) -> Table:
    """Feature rows of the sources ranked in ``sample``'s stratum (rows of its scores table)."""
    feats = Table.read(run_dir / sample / "features.ecsv", format="ascii.ecsv")
    scores = Table.read(run_dir / sample / "scores.ecsv", format="ascii.ecsv")
    keep = np.isin(np.asarray(feats["source_uid"]), np.asarray(scores["source_uid"]))
    out = feats[keep]
    out.meta["provenance"] = "derived"
    out.meta["source"] = f"{run_dir.name}/{sample}: features of the ranked stratum"
    return out


def shared_features(a: Table, b: Table) -> list[str]:
    def numeric(t: Table) -> list[str]:
        return [
            c for c in t.colnames if c not in ID_COLUMNS and np.issubdtype(t[c].dtype, np.number)
        ]

    other = set(numeric(b))
    return [c for c in numeric(a) if c in other]


def as_float(t: Table, col: str) -> np.ndarray:
    c = t[col]
    return np.asarray(c.filled(np.nan) if hasattr(c, "filled") else c, dtype=float)


def values(t: Table, col: str, mask: np.ndarray | None = None) -> np.ndarray:
    x = as_float(t, col)
    if mask is not None:
        x = x[mask]
    return x[np.isfinite(x)]


def compare(x: np.ndarray, y: np.ndarray) -> dict[str, float | int | None]:
    out: dict[str, float | int | None] = {"n_science": int(x.size), "n_control": int(y.size)}
    for name, v in (("science", x), ("control", y)):
        q = np.percentile(v, [16, 50, 84]) if v.size else [np.nan] * 3
        out[f"{name}_q16"], out[f"{name}_median"], out[f"{name}_q84"] = (
            round(float(qq), 4) for qq in q
        )
    if min(x.size, y.size) < MIN_N:
        out["ks_d"] = out["ad_p"] = None
        return out
    out["ks_d"] = round(float(stats.ks_2samp(x, y).statistic), 4)
    with warnings.catch_warnings():  # AD p is capped to [0.001, 0.25] by scipy; keep the cap
        warnings.simplefilter("ignore")
        out["ad_p"] = round(float(stats.anderson_ksamp([x, y]).significance_level), 4)
    return out


def compact_blue_mask(t: Table) -> np.ndarray:
    """Post hoc: compact (low CI_50_30) sources with a stellar-like blue F277W-F356W colour."""
    with np.errstate(invalid="ignore"):
        return (as_float(t, "ref_log_ci_50_30") < COMPACT_LOG_CI) & (
            as_float(t, "color_f277w_f356w") < BLUE_F277W_F356W
        )


def finite(obj):
    """JSON-safe copy: non-finite floats become None."""
    if isinstance(obj, dict):
        return {k: finite(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [finite(v) for v in obj]
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    return obj


def verdict(d_values: list[float]) -> str:
    if not d_values:
        return "untested"
    if max(d_values) > D_STRONG:
        return "differ strongly"
    if max(d_values) <= D_COMPARABLE:
        return "comparable"
    return "intermediate"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run-dir", type=Path, required=True)
    ap.add_argument("--science", default="smacs0723_nircam")
    ap.add_argument("--control", default="ceers_t021")
    ap.add_argument("--out", type=Path, default=Path("results/science_control_features"))
    args = ap.parse_args(argv)

    sci = ranked_features(args.run_dir, args.science)
    ctl = ranked_features(args.run_dir, args.control)
    cols = shared_features(sci, ctl)
    sci_mag, ctl_mag = as_float(sci, "ref_mag"), as_float(ctl, "ref_mag")

    per_feature: dict[str, dict] = {}
    for col in cols:
        entry: dict = {"all": compare(values(sci, col), values(ctl, col)), "bins": {}}
        for lo, hi in MAG_BINS:
            key = f"{lo:g}-{hi:g}"
            ms = (sci_mag > lo) & (sci_mag <= hi)
            mc = (ctl_mag > lo) & (ctl_mag <= hi)
            entry["bins"][key] = compare(values(sci, col, ms), values(ctl, col, mc))
        binned = [b["ks_d"] for b in entry["bins"].values() if b["ks_d"] is not None]
        entry["max_binned_ks_d"] = max(binned) if binned else None
        entry["verdict"] = verdict(binned)
        per_feature[col] = entry

    compact_blue: dict[str, dict[str, dict[str, int]]] = {}
    for name, t, mag in (("science", sci, sci_mag), ("control", ctl, ctl_mag)):
        flags = compact_blue_mask(t)
        compact_blue[name] = {
            f"{lo:g}-{hi:g}": {
                "n": int(((mag > lo) & (mag <= hi)).sum()),
                "compact_blue": int((flags & (mag > lo) & (mag <= hi)).sum()),
            }
            for lo, hi in POSTHOC_BINS
        }

    summary = {
        "provenance": "derived",
        "source": f"run {args.run_dir.name}: {args.science} vs {args.control} ranked strata",
        "science_only_features": [c for c in shared_features(sci, sci) if c not in cols],
        "control_only_features": [c for c in shared_features(ctl, ctl) if c not in cols],
        "n_ranked": {"science": len(sci), "control": len(ctl)},
        "assumptions": {
            "mag_bins": [f"{lo:g}-{hi:g}" for lo, hi in MAG_BINS],
            "d_comparable": D_COMPARABLE,
            "d_strong": D_STRONG,
            "min_n": MIN_N,
        },
        "features": per_feature,
        "posthoc_compact_blue": {
            "definition": f"log10 CI_50_30 < {COMPACT_LOG_CI} and F277W-F356W < {BLUE_F277W_F356W}",
            "counts": compact_blue,
        },
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(
        json.dumps(finite(summary), indent=1, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    plot(sci, ctl, cols, sci_mag, ctl_mag, args.out / "feature_distributions.png")
    for col, e in per_feature.items():
        a = e["all"]
        print(
            f"{col:28s} sci {a['science_median']!s:>8} ctl {a['control_median']!s:>8} "
            f"D_all {a['ks_d']!s:>6} max_bin_D {e['max_binned_ks_d']!s:>6} {e['verdict']}"
        )
    return 0


def plot(sci, ctl, cols, sci_mag, ctl_mag, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(cols)
    ncol = 4
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 2.6 * nrow), squeeze=False)
    faint = (24.0, 26.0)  # the bin where both fields are well populated
    for ax, col in zip(axes.flat, cols, strict=False):
        x_all, y_all = values(sci, col), values(ctl, col)
        if not (x_all.size and y_all.size):
            ax.set_title(col, fontsize=8)
            continue
        lo, hi = np.percentile(np.concatenate([x_all, y_all]), [1, 99])
        bins = np.linspace(lo, hi, 30) if hi > lo else 10
        for mask_s, mask_c, ls, lab in (
            (None, None, "-", "all"),
            (
                (sci_mag > faint[0]) & (sci_mag <= faint[1]),
                (ctl_mag > faint[0]) & (ctl_mag <= faint[1]),
                ":",
                "24-26",
            ),
        ):
            for v, colour, name in (
                (values(sci, col, mask_s), "C3", "SMACS"),
                (values(ctl, col, mask_c), "C0", "CEERS"),
            ):
                if v.size:
                    ax.hist(
                        v,
                        bins=bins,
                        density=True,
                        histtype="step",
                        color=colour,
                        ls=ls,
                        label=f"{name} {lab}",
                    )
        ax.set_title(col, fontsize=8)
        ax.tick_params(labelsize=6)
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    axes.flat[0].legend(fontsize=5)
    fig.tight_layout()
    fig.savefig(path, dpi=80)


if __name__ == "__main__":
    raise SystemExit(main())
