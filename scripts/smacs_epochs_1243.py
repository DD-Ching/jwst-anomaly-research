"""Two-epoch F444W check of SMACS ``f356w_1243`` (2736 o001 -> VENUS 6882 o057).

Plan: ``results/smacs_jwst3/epochs_1243/PLAN.md`` (frozen in the claim commit). Reuses the D-017
tool (``compare_epochs``) and its peer null (``peer_null``, D-017 addenda 1 and 4). ``decide``
encodes the plan's outcome rules.

    python scripts/smacs_epochs_1243.py <scratch dir> --cat-dir <dir with the CATs> \\
        --out results/smacs_jwst3/epochs_1243
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

from jwst_anomaly import catalog, cutouts

sys.path.insert(0, str(Path(__file__).resolve().parent))
from epoch_compare import compare_epochs, peer_null  # noqa: E402
from smacs_epochs_1816_384 import (  # noqa: E402
    BASELINE_YR,
    NULL_DMAG,
    NULL_RADIUS_ARCSEC,
    NULL_SIZE_FRAC,
)

NAME = "f356w_1243"
TARGET = (110.682947, -73.487722)  # record position (DJA / run catalog)
CATS = {
    "e1": "jw02736-o001_t001_nircam_clear-{b}_cat.ecsv",
    "e2": "jw06882-o057_t057_nircam_clear-{b}_cat.ecsv",
}
S3 = {"e1": "jw02736/L3/t/o001/", "e2": "jw06882/L3/t/o057/"}
BAND = "f444w"
# Frozen plan thresholds (ASSUMPTIONs).
MIN_PEERS = 10
VAR_SIG = 5.0
CONST_SIG = 2.0
FADE_BOUND_SIG = 3.0
SHEET_ARCSEC = 4.0
SHEET_PANELS = (
    ("2736 F356W", "e1", "f356w"),
    ("2736 F444W", "e1", "f444w"),
    ("6882 F444W", "e2", "f444w"),
    ("6882 F410M", "e2", "f410m"),
)


def _finite_pos(x: Any) -> bool:
    return x is not None and math.isfinite(x) and x > 0


def decide(st: dict[str, Any] | None) -> dict[str, str]:
    """Plan outcome rules on the F444W statistics ``st`` (None if the target is not measured).

    Variability: |dm_target| >= VAR_SIG x sigma_peer -> "variable"; <= CONST_SIG x sigma_peer ->
    "constant"; else "undetermined"; "not measured" without a target dm. Motion: residual > peer
    p99 -> "moves"; <= peer p95 -> "no motion detected"; else "undetermined". Fewer than
    MIN_PEERS peers -> both "undetermined"; a zero or missing peer sigma / percentile ->
    "undetermined" for that axis.
    """
    if st is None:
        return {"variability": "not measured", "motion": "not measured"}
    if st["n_peers"] < MIN_PEERS:
        return {"variability": "undetermined", "motion": "undetermined"}
    sig = st["peer_dmag_sigma"]
    if st["target_dmag"] is None:
        var = "not measured"
    elif not _finite_pos(sig):
        var = "undetermined"
    else:
        dm = abs(st["target_dmag"])
        if dm >= VAR_SIG * sig:
            var = "variable"
        elif dm <= CONST_SIG * sig:
            var = "constant"
        else:
            var = "undetermined"
    r, p95, p99 = (
        st["target_residual_mas"],
        st["peer_residual_p95_mas"],
        st["peer_residual_p99_mas"],
    )
    if not (_finite_pos(p95) and _finite_pos(p99)):
        mot = "undetermined"
    elif r > p99:
        mot = "moves"
    elif r <= p95:
        mot = "no motion detected"
    else:
        mot = "undetermined"
    return {"variability": var, "motion": mot}


def stats(cat1: Table, cat2: Table, target: SkyCoord, cmp: dict[str, Any]) -> dict[str, Any]:
    """Target values plus the peer null, with the plan's derived bounds."""
    st = peer_null(
        cat1,
        cat2,
        target,
        cmp,
        radius_arcsec=NULL_RADIUS_ARCSEC,
        dmag=NULL_DMAG,
        size_frac=NULL_SIZE_FRAC,
    )
    sig, p95 = st["peer_dmag_sigma"], st["peer_residual_p95_mas"]
    st.update(
        {
            "target_dmag": cmp["target_dmag"],
            "target_residual_mas": cmp["target_residual_total_mas"],
            "target_residual_vec_mas": cmp["target_residual_mas"],
            "constant_threshold_mag": CONST_SIG * sig if sig is not None else None,
            "fade_bound_mag": FADE_BOUND_SIG * sig if sig is not None else None,
            "p95_motion_bound_mas_per_yr": p95 / BASELINE_YR if p95 is not None else None,
        }
    )
    return st


def sheet(path: Path, files: dict[str, str | None], st: dict[str, Any] | None) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.io import fits

    fig, axes = plt.subplots(1, len(files), figsize=(3.2 * len(files), 3.4))
    for ax, (lab, p) in zip(np.atleast_1d(axes), files.items(), strict=True):
        ax.set_xticks([])
        ax.set_yticks([])
        if not p:
            ax.set_title(f"{lab}: no cutout")
            continue
        d = fits.getdata(p).astype(float)
        lo, hi = np.nanpercentile(d, [5, 99.5])
        ax.imshow(d, origin="lower", cmap="gray_r", vmin=lo, vmax=hi)
        ax.plot(d.shape[1] / 2 - 0.5, d.shape[0] / 2 - 0.5, "r+", ms=14, mew=0.8)
        ax.set_title(lab, fontsize=9)
    note = ""
    if st is not None and st["target_dmag"] is not None and st["peer_dmag_sigma"] is not None:
        note = f"dm {st['target_dmag']:+.3f}, sigma_peer {st['peer_dmag_sigma']:.3f}"
    fig.suptitle(f'{NAME} ({SHEET_ARCSEC:.0f}" cutouts, native orientation) {note}', fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("work", type=Path, help="scratch directory for cutout FITS files")
    ap.add_argument("--cat-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    target = SkyCoord(*TARGET, unit="deg")
    cat1 = catalog.load_pipeline_catalog(args.cat_dir / CATS["e1"].format(b=BAND))
    cat2 = catalog.load_pipeline_catalog(args.cat_dir / CATS["e2"].format(b=BAND))
    out: dict[str, Any] = {
        "plan": "results/smacs_jwst3/epochs_1243/PLAN.md",
        "target": list(TARGET),
        "band": BAND,
        "epochs": {"e1": "2736 o001 (2022-06)", "e2": "6882 o057 (2026-06-05)"},
        "files": {
            ep: {"file": CATS[ep].format(b=BAND), "sha256": c.meta.get("input_sha256")}
            for ep, c in (("e1", cat1), ("e2", cat2))
        },
    }
    st: dict[str, Any] | None = None
    try:
        cmp = compare_epochs(cat1, cat2, target)
    except ValueError as exc:
        out["not_measured"] = str(exc)
    else:
        st = stats(cat1, cat2, target, cmp)
        out["compare"] = cmp
        out["peers"] = st
        dm, pm = cmp["target_dmag"], st["peer_dmag_median"]
        out["post_hoc"] = {
            # The frozen rule centres dm on the bright-reference median; these check that choice.
            "dmag_over_reference_sigma": abs(dm) / cmp["reference_dmag_sigma"]
            if dm is not None and cmp["reference_dmag_sigma"] > 0
            else None,
            "dmag_minus_peer_median": dm - pm if dm is not None and pm is not None else None,
            "label": "derived (post hoc)",
        }
    out["decision"] = decide(st)
    out["assumptions"] = {
        "min_peers": MIN_PEERS,
        "var_sig": VAR_SIG,
        "const_sig": CONST_SIG,
        "fade_bound_sig": FADE_BOUND_SIG,
        "null_dmag": NULL_DMAG,
        "null_size_frac": NULL_SIZE_FRAC,
        "null_radius_arcsec": NULL_RADIUS_ARCSEC,
        "baseline_yr": BASELINE_YR,
    }
    tab = Table({"source_uid": [NAME], "ra": [TARGET[0]], "dec": [TARGET[1]]})

    def fetch(panel: tuple[str, str, str]) -> str | None:
        _, ep, b = panel
        i2d = CATS[ep].format(b=b).replace("_cat.ecsv", "_i2d.fits")
        got = cutouts.make_cutouts(
            f"s3://stpubdata/jwst/public/{S3[ep]}{i2d}",
            tab,
            size_arcsec=SHEET_ARCSEC,
            out_dir=args.work / f"{ep}_{b}",
            storage_options={"anon": True},
        )
        return got["path"][0] or None

    with ThreadPoolExecutor(max_workers=len(SHEET_PANELS)) as pool:  # independent S3 range reads
        paths = list(pool.map(fetch, SHEET_PANELS))
    files = {lab: p for (lab, _, _), p in zip(SHEET_PANELS, paths, strict=True)}
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")
    sheet(args.out / "contact_sheet.png", files, st)
    print(json.dumps({"decision": out["decision"], "post_hoc": out.get("post_hoc")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
