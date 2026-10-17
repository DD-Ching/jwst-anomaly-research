"""Bounded comparison: ``arc_diff.py`` with its kernel matching replaced by ordinary SFFT.

Development tool for ``results/transients/sfft_trial/PLAN.md`` (owner, 2026-10-17); not part of
the search chain. ``arc_diff.kernel_match`` is replaced by ``sfft_kernel_match``, which runs
``sfft.CustomizedPacket.Customized_Packet.CP``: ordinary SFFT, polynomial spatial variation,
Numpy CPU backend. It is *not* the B-spline JWST method of Hu & Wang 2024, which is GPU-only in
sfft 1.7.3. Every other step of ``arc_diff`` runs unchanged. Settings are fixed in PLAN §2:
REF = the side ``arc_diff`` convolves, ``ForceConv='REF'``, ``GKerHW = kernel_px``,
``KerPolyOrder=2``, ``BGPolyOrder=0``, ``ConstPhotRatio=True``, no source masking; both masked
images are 0 on the union NaN mask grown by the kernel footprint. Variance: var_ref propagated
through an effective whole-box 9x9 kernel fitted to SFFT's own matched image, plus var_sci
(a mean kernel, so it is biased low where SFFT's local kernel departs from the mean; PLAN
addendum 1). ``significance`` rescales S by its robust std in both tools.

    python scripts/sfft_trial.py <arc_diff.py arguments ...>
    python scripts/sfft_trial.py regions --maps M.npz --hits H.ecsv [--spike RA DEC] --out R.json
    python scripts/sfft_trial.py summarize --runs DIR [--arc-runs DIR] --out S.json
    python scripts/sfft_trial.py diagnostics --run DIR [--run DIR ...] --out D.json

Requires ``pip install sfft==1.7.3`` (not a project dependency).
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import arc_diff as ad  # noqa: E402

SFFT_SETTINGS = {"KerPolyOrder": 2, "BGPolyOrder": 0, "ConstPhotRatio": True}


def sfft_kernel_match(e1, e2, v1, v2, blur, half, stamps=False, stamp_fit=None, timing=None):
    """Drop-in for ``arc_diff.kernel_match``: matched e1, e2, var, kernel, scale and epoch."""
    from sfft.CustomizedPacket import Customized_Packet

    t0 = time.perf_counter()
    conv1 = blur >= 0
    ref, tgt = (e1, e2) if conv1 else (e2, e1)
    vref, vtgt = (v1, v2) if conv1 else (v2, v1)
    # Masked images: both set to 0 on the union of NaNs, grown by the kernel half-width, as
    # SFFT's own packets do; no source masking (PLAN §2; union mask: post-review correction).
    foot = np.ones((2 * half + 1, 2 * half + 1), bool)
    bad = ndimage.binary_dilation(~(np.isfinite(ref) & np.isfinite(tgt)), structure=foot)
    with tempfile.TemporaryDirectory() as d:
        paths = {}
        for name, arr in (("ref", ref), ("sci", tgt)):
            paths[name] = os.path.join(d, name + ".fits")
            paths["m" + name] = os.path.join(d, "m" + name + ".fits")
            fits.writeto(paths[name], arr.astype(float), overwrite=True)
            fits.writeto(paths["m" + name], np.where(bad, 0.0, arr.astype(float)), overwrite=True)
        diff_path = os.path.join(d, "diff.fits")
        Customized_Packet.CP(
            FITS_REF=paths["ref"],
            FITS_SCI=paths["sci"],
            FITS_mREF=paths["mref"],
            FITS_mSCI=paths["msci"],
            ForceConv="REF",
            GKerHW=half,
            FITS_DIFF=diff_path,
            BACKEND_4SUBTRACT="Numpy",
            NUM_CPU_THREADS_4SUBTRACT=1,
            VERBOSE_LEVEL=0,
            **SFFT_SETTINGS,
        )
        diff = fits.getdata(diff_path).astype(float)
    m = tgt - diff  # SFFT DIFF = SCI - K * REF, so K * REF = SCI - DIFF
    m[~np.isfinite(ref)] = np.nan
    # Effective whole-box kernel of SFFT's own matching (m = K * ref + B is noise-free in ref,
    # so a 9x9 least-squares fit recovers K's mean shape; its constant term absorbs B). It
    # propagates the variance as arc_diff does, gives the kernel sum without the background and
    # convolves post-matching injections (post-review correction).
    kern, _ = ad.kernel_fit(ref, m, vref + vtgt, half)
    vfill = np.where(np.isfinite(vref), vref, 0.0)
    var = ndimage.convolve(vfill, kern**2, mode="constant") + vtgt
    var[~np.isfinite(vref)] = np.nan
    ksum = float(kern.sum())
    if timing is not None:
        timing["sfft_s"] = round(time.perf_counter() - t0, 3)
    if conv1:
        return m, e2, var, kern, ksum, 1
    return e1, m, var, kern, 1.0 / ksum, 2


def regions(argv) -> int:
    import argparse

    from astropy.table import Table
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales

    p = argparse.ArgumentParser()
    p.add_argument("regions")
    p.add_argument("--maps", type=Path, required=True)
    p.add_argument("--hits", type=Path, required=True)
    p.add_argument("--spike", type=float, nargs=2)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args(argv)
    z = np.load(a.maps, allow_pickle=True)
    h = Table.read(a.hits)
    out = {"label": "derived; S statistics by region (PLAN §3.3)", "hits": len(h)}
    out["unflagged"] = int((~np.asarray(h["near_bright_core"], bool)).sum()) if len(h) else 0
    for b in sorted({k.split("_")[0] for k in z.keys()}):
        s, e1 = z[b + "_s"], z[b + "_e1img"]
        w = WCS(fits.Header.fromstring(str(z[b + "_wcs"])))
        scale = float(np.mean(proj_plane_pixel_scales(w))) * 3600
        fin = np.isfinite(e1)
        thr = np.nanpercentile(e1, 99.9)
        rc = int(np.ceil(0.8 / scale))
        disk = np.hypot(*np.mgrid[-rc : rc + 1, -rc : rc + 1]) <= rc
        cores = ndimage.binary_dilation(fin & (np.nan_to_num(e1) > thr), structure=disk)
        stampm = ad.compact_stamps(e1, 4)
        regs = {"all": np.isfinite(s), "bright_cores": cores, "compact_stamps": stampm}
        if a.spike:
            x, y = w.all_world2pix(a.spike[0], a.spike[1], 0)
            yy, xx = np.mgrid[0 : s.shape[0], 0 : s.shape[1]]
            r = np.hypot(xx - float(x), yy - float(y)) * scale
            regs["spike_annulus"] = (r >= 1.5) & (r <= 4.0)
        res = {}
        for name, m in regs.items():
            v = s[m & np.isfinite(s)]
            res[name] = {
                "n_px": int(v.size),
                "robust_std": ad._robust_std(v) if v.size else None,
                "max_abs": float(np.max(np.abs(v))) if v.size else None,
                "frac_abs_gt5": float(np.mean(np.abs(v) > 5)) if v.size else None,
            }
        out[b] = res
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out)[:400])
    return 0


def _stage_rows(run: Path) -> dict:
    """fade_calib's stage rows (per-band S included, so 'mixed' signs are judged as in eval)."""
    import csv

    import fade_calib as fc

    rows = fc.stage_rows(run, run / "morph_all.json")
    with open(run / "injected_pre.csv", newline="") as f:
        names = [r["name"] for r in csv.DictReader(f)]
    return {(run.name.split("_", 1)[1], n): r for n, r in zip(names, rows, strict=True)}


def summarize(argv) -> int:
    """PLAN §4 decision inputs from the stored runs (``runs/<tool>_<name>/``)."""
    import argparse

    from astropy.table import Table

    p = argparse.ArgumentParser()
    p.add_argument("summarize")
    p.add_argument(
        "--runs", type=Path, required=True, help="directory with sfft_* (and arc_*) runs"
    )
    p.add_argument("--arc-runs", type=Path, help="directory with arc_* runs (default: --runs)")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args(argv)
    names = ["b1_s30", "b1_s31", "b2_s40"]
    base = {"arc": a.arc_runs or a.runs, "sfft": a.runs}
    tot, per = {}, {}
    for t in ("arc", "sfft"):
        rows = {}
        for n in names:
            rows.update(_stage_rows(base[t] / f"{t}_{n}"))
        per[t] = rows
        fades = [r for r in rows.values() if r["cls"] in ("partial", "vanish")]
        tot[t] = {
            "unchanged_false_hits": sum(r["hit"] for r in rows.values() if r["cls"] == "unchanged"),
            "unchanged_n": sum(r["cls"] == "unchanged" for r in rows.values()),
            "fade_review_sign": sum(r["review_sign"] for r in fades),
            "fade_n": len(fades),
            "appear_review_sign": sum(
                r["review_sign"] for r in rows.values() if r["cls"] == "appear"
            ),
            "appear_n": sum(r["cls"] == "appear" for r in rows.values()),
        }
        scr = 0
        for b in ("b1", "b2"):
            d = base[t] / f"{t}_{b}_screen"
            h = Table.read(d / "hits.ecsv")
            v = {
                x["id"]: x["verdict"]
                for x in json.loads((d / "morph_all.json").read_text())["targets"]
            }
            scr += sum((not r["near_bright_core"]) and v.get(r["hit_id"]) == "pass" for r in h)
        tot[t]["screen_review_entries"] = scr
        rss = {
            n: json.loads((base[t] / f"{t}_{n}" / "rss.json").read_text())
            for n in [*names, "b1_screen", "b2_screen"]
        }
        tot[t]["wall_s_total"] = round(sum(r["wall_s"] for r in rss.values()), 1)
        tot[t]["peak_rss_mb"] = {n: r["peak_rss_mb"] for n, r in rss.items()}
    diffs = []
    for k, ra_ in per["arc"].items():
        rs = per["sfft"][k]
        if ra_["review_sign"] != rs["review_sign"] or ra_["hit"] != rs["hit"]:
            diffs.append(
                {
                    "run": k[0],
                    "name": k[1],
                    "cls": ra_["cls"],
                    "mag": ra_["mag"],
                    "arc": {"hit": ra_["hit"], "review_sign": ra_["review_sign"]},
                    "sfft": {"hit": rs["hit"], "review_sign": rs["review_sign"]},
                }
            )
    ar, sf = tot["arc"], tot["sfft"]
    mem_ratio = [round(sf["peak_rss_mb"][n] / ar["peak_rss_mb"][n], 2) for n in ar["peak_rss_mb"]]
    crit = {
        "a_false_hits": sf["unchanged_false_hits"] <= ar["unchanged_false_hits"],
        "b_fades": sf["fade_review_sign"] >= ar["fade_review_sign"]
        and sf["appear_review_sign"] >= ar["appear_review_sign"] - 2,
        "c_screen_review": sf["screen_review_entries"] <= ar["screen_review_entries"],
        "d_wall": sf["wall_s_total"] <= 2 * ar["wall_s_total"],
    }
    if all(crit.values()):
        verdict = "try SFFT in the next batch"
    elif not crit["a_false_hits"] or not crit["b_fades"]:
        verdict = "keep arc_diff"
    else:
        verdict = "evidence insufficient"
    out = {
        "label": "derived; PLAN §4 rule (ordinary SFFT, not the B-spline JWST method)",
        "totals": tot,
        "criteria": crit,
        "verdict": verdict,
        "peak_rss_ratio_sfft_over_arc": [min(mem_ratio), max(mem_ratio)],
        "per_source_differences": diffs,
    }
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({"totals": tot, "criteria": crit, "verdict": verdict, "n_diff": len(diffs)}))
    return 0


def diagnostics(argv) -> int:
    """Median per-band S of injected sources per class x magnitude, per run."""
    import argparse
    import csv
    import statistics

    p = argparse.ArgumentParser()
    p.add_argument("diagnostics")
    p.add_argument("--run", type=Path, action="append", required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args(argv)
    out = {
        "label": "derived; median S at injected positions (injected_pre.csv), per run",
        "runs": {},
    }
    for run in a.run:
        with open(run / "injected_pre.csv", newline="") as f:
            rows = list(csv.DictReader(f))
        bands = [k for k in rows[0] if k.isupper() and k.startswith("F")]
        res = {}
        for c in sorted({r["cls"] for r in rows}):
            for m in sorted({r["mag"] for r in rows}):
                sel = [r for r in rows if r["cls"] == c and r["mag"] == m]
                if sel:
                    res[f"{c}_{m}"] = {
                        b: round(
                            statistics.median(float(r[b]) for r in sel if r[b] not in ("", "None")),
                            1,
                        )
                        for b in bands
                    }
        out["runs"][str(run)] = res
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    return 0


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "regions":
        return regions(argv)
    if argv and argv[0] == "summarize":
        return summarize(argv)
    if argv and argv[0] == "diagnostics":
        return diagnostics(argv)
    ad.kernel_match = sfft_kernel_match
    return ad.main(argv)


if __name__ == "__main__":
    sys.exit(main())
