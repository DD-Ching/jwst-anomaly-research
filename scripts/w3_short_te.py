"""W3 short-event recovery study (owner directive V2.1, priority 1): stage 1, a frozen sample.

Runs the production injection chain (``w3_moa.run_inject``, CHAIN_VERSION unchanged) for t_E = 3 d
only, on one field, with a seed that no limit uses, and writes a tracked table with a frozen
``split`` column (dev / validation, by crc32 of the row index and seed) fixed before any analysis.
Later stages tune classifiers on ``dev`` only and report recovery on ``validation`` once.

    python scripts/w3_short_te.py --field gb12 --per-cell 1000 --seed 3003 --procs 4
    python scripts/w3_short_te.py --field gb12 --seed 3003 --diagnose   # dev split only
    python scripts/w3_short_te.py --field gb12 --seed 3003 --revet --procs 4   # stage 2, dev only
    python scripts/w3_short_te.py --field gb12 --seed 3003 --revet --night-robust --procs 4
    python scripts/w3_short_te.py --field gb7 --real gb7-R-8-6-94052 --night-robust   # FP check
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import w3_moa as wm  # noqa: E402
from astropy.table import Table  # noqa: E402

OUT = ROOT / "results" / "w3_moa" / "short_te"


SPIKE_SIGMA = 3.0  # ASSUMPTION: spike / umbra epochs: injected flux beyond 3 median errors
JOB_KEYS = ("event_id", "tE", "rho", "u0", "t0", "Is")  # the inputs that define one injection


NIGHT_MARGIN_D = (
    1.0  # ASSUMPTION: nights within 1 d of the feature are left out of the σ_night estimate
)


def night_chi2(r, w, nt, s2: float) -> float:
    """χ² of residuals ``r`` (weights ``w``) when every exposure of a night shares one Gaussian
    offset of variance ``s2`` (per-night systematic error, marginalised): per night,
    Σwr² − s2 (Σwr)² / (1 + s2 Σw) (Sherman–Morrison). The log-determinant is the same for every
    model on the same light curve, so it cancels in a ΔBIC."""
    import numpy as np

    _, inv = np.unique(nt, return_inverse=True)
    a = np.bincount(inv, w)
    b = np.bincount(inv, w * r)
    c = np.bincount(inv, w * r * r)
    return float((c - s2 * b * b / (1.0 + s2 * a)).sum())


def night_variance(r, w, nt) -> float:
    """Maximum-likelihood σ_night² from per-night weighted mean residuals ``r`` (the light curve
    outside the feature): m_n ~ N(0, 1/Σw + s2). ``derived``; 0 when the nights scatter as white
    noise."""
    import numpy as np
    from scipy.optimize import minimize_scalar

    _, inv = np.unique(nt, return_inverse=True)
    a = np.bincount(inv, w)
    m = np.bincount(inv, w * r) / a
    v0 = 1.0 / a

    def nll(ls):
        v = v0 + 10.0**ls
        return float(np.sum(np.log(v) + m * m / v))

    hi = float(np.log10(max(np.var(m), 1e-12))) + 1.0
    ls = minimize_scalar(nll, bounds=(hi - 12.0, hi), method="bounded").x
    s2 = 10.0**ls
    return 0.0 if nll(ls) >= nll(hi - 12.0) - 1e-6 else float(s2)


def night_robust(lc, res: dict, ex: str, best_o: str, ecl: dict) -> dict:
    """Stage 2 step 2 (D-082): the exotic preference under a per-night systematic error instead of
    dropping nights. σ_night² is fitted to the best ordinary model's residuals outside the feature;
    ΔBIC = χ²_night(exotic) − χ²_night(alt) + the BIC penalty difference, against the best single
    lens and against the eclipse dip, at the white-noise best fits (no refit: an approximation);
    the larger (least favourable to the exotic model) is returned."""
    import numpy as np

    fe = wm.w3.model_flux(ex, lc, res[ex])
    fo = wm.w3.model_flux(best_o, lc, res[best_o])
    nt = wm.nights(lc.t)
    inside = np.abs(fe - fo) > 3.0 * np.median(lc.sf)
    if not inside.any():
        return {
            "nr_s_night": np.nan,
            "nr_dbic": np.nan,
            "nr_dbic_lens": np.nan,
            "nr_dbic_ecl": np.nan,
        }
    lo, hi = float(lc.t[inside].min()), float(lc.t[inside].max())
    out_ = (lc.t < lo - NIGHT_MARGIN_D) | (lc.t > hi + NIGHT_MARGIN_D)
    s2 = night_variance((lc.f - fo)[out_], lc.w[out_], nt[out_])
    fecl = wm.eclipse_flux(lc, ecl)
    pen_e = res[ex]["bic"] - res[ex]["chi2"]
    ce = night_chi2(lc.f - fe, lc.w, nt, s2)
    d_lens = (
        ce
        - night_chi2(lc.f - fo, lc.w, nt, s2)
        + pen_e
        - (res[best_o]["bic"] - res[best_o]["chi2"])
    )
    d_ecl = ce - night_chi2(lc.f - fecl, lc.w, nt, s2) + pen_e - (ecl["bic"] - ecl["chi2"])
    return {"nr_s_night": float(np.sqrt(s2)), "nr_dbic": float(max(d_lens, d_ecl)),
            "nr_dbic_lens": float(d_lens), "nr_dbic_ecl": float(d_ecl)}  # fmt: skip


def patch_night_robust() -> None:
    """Record ``night_robust`` next to the production jackknife (its result and the vetting are
    unchanged): the jackknife's returned dict carries the extra keys into ``vet_one``'s record."""
    orig = wm.jackknife_nights
    if getattr(orig, "_night_robust", False):
        return

    def both(lc, res, ex, best_o, ecl):
        out = orig(lc, res, ex, best_o, ecl)
        out.update(night_robust(lc, res, ex, best_o, ecl))
        return out

    both._night_robust = True
    wm.jackknife_nights = both


def split_of(row, seed: int) -> str:
    """Frozen assignment (ASSUMPTION: 50/50) from a SHA-256 of the injection's own inputs and the
    seed. It does not depend on the table's row order (pool completion order) or on any outcome."""
    key = f"{seed}|" + "|".join(repr(float(row[k])) if k != "event_id" else str(row[k])
                                for k in JOB_KEYS)  # fmt: skip
    return "validation" if hashlib.sha256(key.encode()).digest()[0] & 1 else "dev"


def sampling_diagnostics(tab: Table, field: str) -> Table:
    """Per flagged injection (``derived``): how well the injected signal is sampled. Spike nights
    have an epoch where the injected flux rises above 3 median errors; umbra nights one where it
    falls below −3 median errors (the jackknife's own feature threshold)."""
    import numpy as np

    wm.set_field(field)
    pre = wm.read_prescreen()
    flagged = tab[np.asarray(tab["flagged"], bool)]
    arrays = wm.load_arrays(list(dict.fromkeys(flagged["event_id"])), aux=False, pre=pre)
    rows = []
    for r in flagged:
        t, f, sf = arrays[r["event_id"]][:3]
        fs = float(wm.moa.mag_to_counts(float(r["Is"]), wm.moa.parse_event_id(r["event_id"])[1]))
        prm = {k: float(r[k]) for k in ("tE", "rho", "u0", "t0")}
        sig = wm.injected_signal(t, "W3", fs, prm)
        thr = SPIKE_SIGMA * float(np.median(sf))
        nt = wm.nights(t)
        up, down = sig > thr, sig < -thr
        tc = prm["t0"]
        rows.append(
            {
                "event_id": r["event_id"],
                "split": r["split"],
                "rho": prm["rho"],
                "u0": prm["u0"],
                "Is": float(r["Is"]),
                "failed_test": str(r["failed_test"] or "")
                if not np.ma.is_masked(r["failed_test"])
                else "",
                "recovered": bool(r["recovered"]),
                "spike_nights": int(np.unique(nt[up]).size),
                "spike_before": bool((up & (t < tc)).any()),
                "spike_after": bool((up & (t > tc)).any()),
                "umbra_nights": int(np.unique(nt[down]).size),
                "feature_nights": int(np.unique(nt[up | down]).size),
            }
        )
    out = Table(rows=rows)
    out.meta.update(
        provenance="derived",
        source=f"scripts/w3_short_te.py --diagnose: injected-signal sampling of flagged {field} "
        "t_E = 3 d injections (simulated signals on real carriers)",
    )
    return out


def _revet_worker(job):
    """Re-run the production vetting on one injected light curve (``_inject_worker``'s inputs) and
    keep what the jackknife itself used: its feature nights, the nights it dropped and its ΔBIC."""
    import numpy as np

    r, t, f, sf, ax = job
    eid = r["event_id"]
    prm = {k: float(r[k]) for k in ("tE", "rho", "u0", "t0")}
    sig = wm.injected_signal(
        t, "W3", float(wm.moa.mag_to_counts(r["Is"], wm.moa.parse_event_id(eid)[1])), prm
    )
    f2 = f + sig
    scan = wm.deficit_scan(t, f2, sf)
    d = {"event_id": eid, "ra": r["ra"], "dec": r["dec"], "ref_mag": r["Is"]}
    o = wm.vet_one((d, t, f2, sf, ax, {k: scan[k] for k in wm.SCAN_KEYS}, [], False, None))
    failed = next(((nm, note) for nm, ok, note in o["tests"] if not ok), ("", ""))
    jn = o.get("jackknife_nights", {})
    nt = wm.nights(t)
    thr = SPIKE_SIGMA * float(np.median(sf))
    spike = set(np.unique(nt[sig > thr]).tolist())
    dropped = [int(x) for x in jn.get("dropped", [])]
    return {
        "event_id": eid,
        "rho": prm["rho"],
        "u0": prm["u0"],
        "Is": float(r["Is"]),
        "stage1_failed": "" if np.ma.is_masked(r["failed_test"]) else str(r["failed_test"]),
        "failed_test": failed[0],
        "failed_note": failed[1][:200],
        "survives": bool(o["survives"]),
        "exotic": str(o.get("exotic", "")),
        "dbic_all": float(o.get("dbic_all", np.nan)),
        "jk_ran": bool(jn),
        "jk_feature_nights": int(jn.get("n_nights", -1)),
        "jk_n_dropped": int(jn.get("n_dropped", -1)),
        "jk_dbic": float(jn.get("dbic", np.nan)),
        "jk_dropped": ",".join(map(str, dropped)),
        "spike_nights": len(spike),
        "dropped_spike_nights": sum(x in spike for x in dropped),
        "spikes_left": len(spike - set(dropped)),
        **{k: float(jn.get(k, np.nan)) for k in NR_KEYS},
    }


NR_KEYS = ("nr_s_night", "nr_dbic", "nr_dbic_lens", "nr_dbic_ecl")


def vet_real(field: str, eid: str) -> dict:
    """The production vetting of one real flag (``run_vet``'s job, without neighbours, binary-lens
    fit or catalogue queries: the record stops where the chain stops)."""
    import numpy as np

    wm.use_moa_fit_params()
    wm.set_field(field)
    pre = wm.read_prescreen()
    wm.set_field_context(pre)
    from astropy.table import vstack

    chunks = sorted((ROOT / "results" / "w3_moa").glob(f"fits_{field}_chunk*.ecsv.gz"))
    fits = vstack([Table.read(c, format="ascii.ecsv") for c in chunks])
    r = fits[np.asarray(fits["event_id"], str) == eid][0]
    ev = wm.moa.MoaField(field).events()
    pos = ev[np.asarray(ev["event_id"], str) == eid][0]
    t, f, sf, ax = wm.load_arrays([eid], aux=True, pre=pre)[eid]
    rp, _ = wm.gaia_rp_mags([float(pos["ra"])], [float(pos["dec"])])
    d = {"event_id": eid, "ra": float(pos["ra"]), "dec": float(pos["dec"]),
         "dophot_mag": float(pos["dophot_magnitude"]), "ref_mag": float(rp[0])}  # fmt: skip
    scan = {k: float(r[k]) for k in wm.SCAN_KEYS}
    return wm.vet_one((d, t, f, sf, ax, scan, [], False, wm.res_from_row(r)))


def _filled(col):
    """String column with ECSV's masked empty strings back as ""."""
    import numpy as np

    return np.asarray(col.filled("") if hasattr(col, "filled") else col, str)


def revet_dev(tab: Table, field: str, procs: int) -> Table:
    """Stage 2, step 1 (D-082 "Revisit if"): re-vet the flagged dev injections with the production
    chain and record the jackknife's own feature nights, dropped nights and ΔBIC (``derived``)."""
    import numpy as np

    wm.use_moa_fit_params()
    wm.set_field(field)
    pre = wm.read_prescreen()
    wm.set_field_context(pre)
    sel = tab[np.asarray(tab["flagged"], bool) & (np.asarray(tab["split"], str) == "dev")]
    arrays = wm.load_arrays(list(dict.fromkeys(sel["event_id"])), aux=True, pre=pre)
    jobs = [(dict(zip(sel.colnames, r, strict=True)), *arrays[r["event_id"]]) for r in sel]
    with wm._pool(procs) as pool:
        rows = list(pool.imap(_revet_worker, jobs, 1))
    out = Table(rows=rows)
    out.meta.update(
        provenance="derived",
        source=f"scripts/w3_short_te.py --revet: production vetting re-run on flagged dev {field} "
        "t_E = 3 d injections (simulated signals on real carriers)",
    )
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--field", default="gb12")
    ap.add_argument("--per-cell", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=3003)
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--diagnose", action="store_true", help="sampling diagnostics of the sample")
    ap.add_argument("--revet", action="store_true", help="stage 2: re-vet flagged dev injections")
    ap.add_argument("--night-robust", action="store_true", help="stage 2 step 2: record ΔBIC "
                    "under a per-night systematic error next to the jackknife")  # fmt: skip
    ap.add_argument("--real", default="", help="vet one real flag of --field (FP check)")
    ap.add_argument("--resplit", action="store_true", help="rewrite the split column (job key)")
    a = ap.parse_args(argv)
    if a.night_robust:
        patch_night_robust()
    if a.real:
        import json

        o = vet_real(a.field, a.real)
        print(json.dumps({k: o[k] for k in ("event_id", "tests", "survives", "jackknife_nights")
                          if k in o}, default=str, indent=1))  # fmt: skip
        return 0
    if a.resplit:
        path = OUT / f"injections_{a.field}_te3_seed{a.seed}.ecsv.gz"
        tab = Table.read(path, format="ascii.ecsv")
        tab["split"] = [split_of(r, a.seed) for r in tab]
        wm.w3.write_ecsv_gz(tab, path)
        print(
            f"re-split {path}: {sum(tab['split'] == 'dev')} dev, "
            f"{sum(tab['split'] != 'dev')} validation"
        )
        return 0
    if a.revet:
        tab = Table.read(
            OUT / f"injections_{a.field}_te3_seed{a.seed}.ecsv.gz", format="ascii.ecsv"
        )
        d = revet_dev(tab, a.field, a.procs)
        nr = "_nightrobust" if a.night_robust else ""
        out = OUT / f"revet_dev_{a.field}_te3_seed{a.seed}{nr}.ecsv"
        d.write(out, overwrite=True)
        same = int(sum(_filled(d["failed_test"]) == _filled(d["stage1_failed"])))
        print(f"wrote {out}: {len(d)} dev flags; failing test reproduced for {same}")
        return 0
    if a.diagnose:
        tab = Table.read(
            OUT / f"injections_{a.field}_te3_seed{a.seed}.ecsv.gz", format="ascii.ecsv"
        )
        d = sampling_diagnostics(tab, a.field)
        out = OUT / f"diagnostics_{a.field}_te3_seed{a.seed}.ecsv"
        d.write(out, overwrite=True)
        print(f"wrote {out}: {len(d)} flagged injections")
        return 0
    wm.set_field(a.field)
    prod = wm.out_dir() / f"injections_{a.field}.ecsv"
    keep = prod.with_suffix(".production.ecsv")
    had_prod = prod.exists()
    if keep.exists():
        raise SystemExit(f"{keep} exists (an interrupted run?): restore it to {prod} first")
    if had_prod:
        shutil.copy2(prod, keep)  # run_inject writes to the production path
    te_grid = wm.INJ_TE
    wm.INJ_TE = (3.0,)
    try:
        path = wm.run_inject(a.procs, a.per_cell, 0, a.seed, False, "lf")
        tab = Table.read(path)
    finally:
        wm.INJ_TE = te_grid
        if had_prod:
            shutil.move(keep, prod)
        else:
            prod.unlink(missing_ok=True)  # never leave the study table as production injections
    tab["split"] = [split_of(r, a.seed) for r in tab]
    tab.meta["study"] = "W3 short-event recovery, stage 1 (t_E = 3 d; frozen dev/validation split)"
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"injections_{a.field}_te3_seed{a.seed}.ecsv.gz"
    wm.w3.write_ecsv_gz(tab, out)
    print(f"wrote {out}: {len(tab)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
