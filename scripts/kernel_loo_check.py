"""Kernel-level check of ``arc_diff.kernel_fit_stamps`` (shared) against per-stamp ``kernel_fit``.

Runs ``arc_diff.band_diff`` for one band up to ``kernel_match``, captures its exact inputs, then
computes every stamp's leave-one-out kernel both ways on the same arrays (results/transients/
kernel_loo/PLAN.md §3 test 1). Reports max|ΔK| / max|K| per stamp, kept-row counts are not
exposed by ``kernel_fit``, so only kernels are compared; plus the condition number of the
weighted stamp design and the time of each method.

    python scripts/kernel_loo_check.py --band F150W --e1 <i2d> --e2 <i2d> --ra R --dec D \
        [--max-stamps N] --out check.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import arc_diff as ad  # noqa: E402


class _Captured(BaseException):
    pass


def capture_inputs(job: dict) -> dict:
    box = {}

    def grab(e1, e2, v1, v2, blur, half, stamps=False, **_):
        box.update(e1=e1, e2=e2, v1=v1, v2=v2, blur=blur, half=half)
        raise _Captured

    orig = ad.kernel_match
    ad.kernel_match = grab
    try:
        ad.band_diff(job)
    except _Captured:
        pass
    finally:
        ad.kernel_match = orig
    return box


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--band", required=True)
    p.add_argument("--e1", required=True)
    p.add_argument("--e2", required=True)
    p.add_argument("--ra", type=float, required=True)
    p.add_argument("--dec", type=float, required=True)
    p.add_argument("--half-arcsec", type=float, default=12.0)
    p.add_argument("--kernel-px", type=int, default=4)
    p.add_argument("--max-stamps", type=int, default=0, help="0 = all")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args(argv)
    job = dict(
        band=a.band.upper(), e1=a.e1, e2=a.e2, ra=a.ra, dec=a.dec, half=a.half_arcsec,
        bg_arcsec=1.0, tile_arcsec=3.0, kernel_px=a.kernel_px, kernel_stamps=True,
        core_arcsec=0.8, core_pct=99.9, min_sigma=5.0, inject=0, inject_mags=[], seed=0,
        inj_sky=[], inj_pre=[], reg_tile_arcsec=6.0,
    )  # fmt: skip
    x = capture_inputs(job)
    half = x["half"]
    conv1 = x["blur"] >= 0
    ref, tgt = (x["e1"], x["e2"]) if conv1 else (x["e2"], x["e1"])
    var = x["v1"] + x["v2"]
    stamp = ad.compact_stamps(x["e1"], half) & ad.compact_stamps(x["e2"], half)
    lab, nlab = ndimage.label(stamp)
    t0 = time.perf_counter()
    tim: dict = {}
    shared = ad.kernel_fit_stamps(ref, tgt, var, half, stamp, lab, nlab, timing=tim)
    t_shared = time.perf_counter() - t0
    labels = list(range(1, nlab + 1))
    if a.max_stamps:
        labels = labels[: a.max_stamps]
    rel, skipped_both, mismatch_skip = [], 0, []
    t0 = time.perf_counter()
    for i in labels:
        try:
            k_old, _ = ad.kernel_fit(ref, tgt, var, half, mask=stamp & (lab != i))
        except ValueError:
            skipped_both += i not in shared
            if i in shared:
                mismatch_skip.append(i)
            continue
        if i not in shared:
            mismatch_skip.append(i)
            continue
        rel.append(float(np.max(np.abs(shared[i] - k_old)) / np.max(np.abs(k_old))))
    t_refit = time.perf_counter() - t0
    # condition number of the weighted stamp design (all stamps)
    n = 2 * half + 1
    good = np.isfinite(ref) & np.isfinite(tgt) & np.isfinite(var) & (var > 0)
    good[:half, :] = good[-half:, :] = good[:, :half] = good[:, -half:] = False
    fin = ndimage.minimum_filter(np.isfinite(ref).astype(np.uint8), size=n) > 0
    ys, xs = np.nonzero(good & fin & stamp)
    r0 = np.where(np.isfinite(ref), ref, 0.0)
    cols = [r0[ys - dy, xs - dx] for dy in range(-half, half + 1) for dx in range(-half, half + 1)]
    des = np.column_stack([*cols, np.ones(ys.size)]) / np.sqrt(var[ys, xs])[:, None]
    out = {
        "label": "derived; kernel-level comparison, same inputs (PLAN §3 test 1, tol 1e-6)",
        "band": a.band.upper(),
        "box": [a.ra, a.dec, a.half_arcsec],
        "n_stamps": int(nlab),
        "stamps_compared": len(rel),
        "skips_both": int(skipped_both),
        "skip_mismatch": mismatch_skip,
        "max_rel_dK": max(rel) if rel else None,
        "median_rel_dK": float(np.median(rel)) if rel else None,
        "n_over_tolerance": int(sum(r > 1e-6 for r in rel)),
        "cond_design": float(np.linalg.cond(des)),
        "rows": int(ys.size),
        "time_shared_all_s": round(t_shared, 2),
        "shared_build_s": tim.get("stamp_design_build"),
        "time_refit_compared_s": round(t_refit, 2),
    }
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
