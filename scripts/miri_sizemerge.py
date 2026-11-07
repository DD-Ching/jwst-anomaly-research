"""Size-aware MIRI band association on SMACS 0723 (TASKS Then 7, D-014 addendum 5).

Plan: ``results/miri_confirm/PLAN_sizemerge.md`` (rule frozen before measuring). Reads the four
MAST L3 MIRI ``_cat.ecsv`` of program 2736 o002 and writes ``summary_sizemerge.json``.

Provenance: catalogs are ``observed``; merges, match counts and chance fractions are ``derived``;
r0 = 0.1", the alpha grid, the 8" shifts and the 0.10 chance ceiling are ASSUMPTIONs from the plan.
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
from astropy.table import Table

from jwst_anomaly import catalog

PREFIX = "jw02736-o002_t001_miri_"
BANDS = ("F770W", "F1000W", "F1500W", "F1800W")
R0 = 0.1  # ASSUMPTION (plan): current radius
ALPHAS = (0.0, 0.25, 0.5, 0.75, 1.0)  # ASSUMPTION (plan)
SHIFTS = ((8.0, 8.0), (8.0, -8.0), (-8.0, 8.0), (-8.0, -8.0))  # arcsec, ASSUMPTION (plan)
# Deviation (post hoc fix): each non-reference band gets its own direction (base shift rotated by
# 120 deg per band, same 11.3" length). A common shift keeps the shifted bands aligned with each
# other, so their real mutual pairs counted as chance (93 "chance" pairs at 0.1", not stored).
BAND_ROTATION_DEG = 120.0
MAX_CHANCE = 0.10  # ASSUMPTION (plan)
VETTED = {221: True, 239: True, 118: True, 240: True, 176: False}  # F770W label -> merge expected


def shifted(cat: Table, dx: float, dy: float) -> Table:
    """Copy of ``cat`` with positions offset by (dx, dy) arcsec on the sky (chance test)."""
    out = cat.copy()
    dec0 = np.deg2rad(float(np.nanmedian(out["dec"])))
    out["ra"] = np.asarray(out["ra"], float) + dx / 3600.0 / np.cos(dec0)
    out["dec"] = np.asarray(out["dec"], float) + dy / 3600.0
    return out


def n_pairs(merged: Table) -> int:
    """Number of cross-band associations: sum over rows of (n_bands - 1)."""
    return int(np.sum(np.asarray(merged["n_bands"]) - 1))


POST_HOC_RADII = (0.1, 0.15, 0.2, 0.3, 0.4, 0.6)  # arcsec; post hoc diagnostic, not a rule


def merge(cats: dict[str, Table], alpha: float, r0: float = R0) -> Table:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return catalog.merge_bands(cats, "F770W", radius_arcsec=r0, size_factor=alpha)


def chance_pairs(cats: dict[str, Table], alpha: float, r0: float = R0) -> float:
    """Mean pair count with every non-reference band shifted by each of ``SHIFTS``."""
    counts = []
    for dx, dy in SHIFTS:
        sh = {}
        for i, (b, c) in enumerate(cats.items()):
            if b == "F770W":
                sh[b] = c
                continue
            t = np.deg2rad(BAND_ROTATION_DEG * (i - 1))
            sh[b] = shifted(c, dx * np.cos(t) - dy * np.sin(t), dx * np.sin(t) + dy * np.cos(t))
        counts.append(n_pairs(merge(sh, alpha, r0)))
    return float(np.mean(counts))


def vetted_status(merged: Table) -> dict[str, int]:
    """n_bands of the merged row anchored on each vetted F770W label."""
    lab = np.asarray(merged["f770w_label"])
    return {f"f770w_{k}": int(merged["n_bands"][np.flatnonzero(lab == k)[0]]) for k in VETTED}


def chance_fraction(extra_real: float, extra_chance: float) -> float | None:
    """Chance share of the pairs added over alpha = 0.

    None (never passes) when nothing was added or the chance count fell below alpha = 0 (a net
    loss from one-to-one reassignment is not a measured low chance rate).
    """
    if extra_real <= 0 or extra_chance < 0:
        return None
    return extra_chance / extra_real


def adopt_alpha(alphas, fractions, max_chance: float = MAX_CHANCE) -> float:
    """Frozen rule: the largest alpha > 0 with chance fraction <= ``max_chance``, else 0."""
    ok = [
        al
        for al, f in zip(alphas, fractions, strict=True)
        if al > 0 and f is not None and f <= max_chance
    ]
    return max(ok, default=0.0)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("data_dir", type=Path, help="directory with the four MIRI _cat.ecsv files")
    p.add_argument("--out", type=Path, default=Path("results/miri_confirm/summary_sizemerge.json"))
    a = p.parse_args(argv)
    cats = {
        b: catalog.load_pipeline_catalog(a.data_dir / f"{PREFIX}{b.lower()}_cat.ecsv")
        for b in BANDS
    }
    meta = {
        b: {
            "n": len(c),
            "jwst_version": c.meta["jwst_version"],
            "pixel_scale_arcsec": round(c.meta["pixel_scale_arcsec"], 4),
            "sha256": c.meta["input_sha256"],
        }
        for b, c in cats.items()
    }
    real = {al: merge(cats, al) for al in ALPHAS}
    chance = {al: chance_pairs(cats, al) for al in ALPHAS}
    post_hoc = [
        {
            "radius_arcsec": r,
            "pairs_real": n_pairs(merge(cats, 0.0, r)),
            "pairs_chance_mean": chance_pairs(cats, 0.0, r),
        }
        for r in POST_HOC_RADII
    ]
    rows = []
    fracs = []
    for al in ALPHAS:
        extra_real = n_pairs(real[al]) - n_pairs(real[0.0])
        extra_chance = chance[al] - chance[0.0]
        frac = chance_fraction(extra_real, extra_chance)
        fracs.append(frac)  # unrounded: the rule decides on these, the JSON shows them rounded
        passes = adopt_alpha([al], [frac]) == al > 0
        n_multi = int(np.sum(np.asarray(real[al]["n_bands"]) >= 2))
        rows.append(
            {
                "alpha": al,
                "n_rows": len(real[al]),
                "n_multiband_rows": n_multi,
                "pairs_real": n_pairs(real[al]),
                "pairs_chance_mean": chance[al],
                "extra_real": extra_real,
                "extra_chance": extra_chance,
                "chance_fraction": None if frac is None else round(frac, 4),
                "passes": passes,
                "vetted_n_bands": vetted_status(real[al]),
            }
        )
    adopted = adopt_alpha(ALPHAS, fracs)
    vet = rows[ALPHAS.index(adopted)]["vetted_n_bands"]
    hits = sum(1 for k, exp in VETTED.items() if exp and vet[f"f770w_{k}"] >= 2)
    out = {
        "plan": "results/miri_confirm/PLAN_sizemerge.md",
        "inputs": meta,
        "r0_arcsec": R0,
        "shifts_arcsec": SHIFTS,
        "band_rotation_deg": BAND_ROTATION_DEG,
        "max_chance_fraction": MAX_CHANCE,
        "grid": rows,
        "adopted_alpha": adopted,
        "vetted_expected_merges_hit": f"{hits} / {sum(VETTED.values())}",
        "vetted_176_single": vet["f770w_176"] == 1,
        "post_hoc_fixed_radius": post_hoc,
    }
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    for r in rows:
        print(
            {
                k: r[k]
                for k in (
                    "alpha",
                    "pairs_real",
                    "extra_real",
                    "extra_chance",
                    "chance_fraction",
                    "passes",
                )
            },
            r["vetted_n_bands"],
        )
    for r in post_hoc:
        print(r)
    print("adopted", adopted, "vetted hits", out["vetted_expected_merges_hit"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
