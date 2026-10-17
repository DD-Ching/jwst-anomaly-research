"""Bounded comparison: ``arc_diff.py`` with its kernel matching replaced by ordinary SFFT.

Development tool for ``results/transients/sfft_trial/PLAN.md`` (owner, 2026-10-17); not part of
the search chain. ``arc_diff.kernel_match`` is replaced by ``sfft_kernel_match``, which runs
``sfft.CustomizedPacket.Customized_Packet.CP``: ordinary SFFT, polynomial spatial variation,
Numpy CPU backend. It is *not* the B-spline JWST method of Hu & Wang 2024, which is GPU-only in
sfft 1.7.3. Every other step of ``arc_diff`` runs unchanged. Settings are fixed in PLAN §2:
REF = the side ``arc_diff`` convolves, ``ForceConv='REF'``, ``GKerHW = kernel_px``,
``KerPolyOrder=2``, ``BGPolyOrder=0``, ``ConstPhotRatio=True``, no source masking (NaN -> 0).
Variance is var_ref + var_sci (no kernel propagation); ``significance`` rescales S by its robust
std in both tools.

    python scripts/sfft_trial.py <arc_diff.py arguments ...>
    python scripts/sfft_trial.py regions --maps M.npz --hits H.ecsv [--spike RA DEC] --out R.json

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
    with tempfile.TemporaryDirectory() as d:
        paths = {}
        for name, arr in (("ref", ref), ("sci", tgt)):
            paths[name] = os.path.join(d, name + ".fits")
            paths["m" + name] = os.path.join(d, "m" + name + ".fits")
            fits.writeto(paths[name], arr.astype(float), overwrite=True)
            fits.writeto(paths["m" + name], np.nan_to_num(arr.astype(float)), overwrite=True)
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
    var = np.where(np.isfinite(vref), vref, np.nan) + vtgt
    ksum = ad.flux_scale(ref, m)
    kern = np.zeros((2 * half + 1, 2 * half + 1))
    kern[half, half] = ksum  # summary only: SFFT's kernel varies over the box
    if timing is not None:
        timing["sfft_s"] = round(time.perf_counter() - t0, 3)
    if conv1:
        return m, e2, var, kern, ksum, 1
    return e1, m, var, kern, 1.0 / ksum, 2


def regions(argv) -> int:
    import argparse

    from astropy.table import Table
    from astropy.wcs import WCS

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
        scale = float(np.mean(np.abs(np.diag(w.pixel_scale_matrix)))) * 3600
        fin = np.isfinite(e1)
        thr = np.nanpercentile(e1, 99.9)
        rc = int(np.ceil(0.8 / scale))
        from scipy import ndimage

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
    import csv

    v = {t["id"]: t["verdict"] for t in json.loads((run / "morph_all.json").read_text())["targets"]}
    out = {}
    with open(run / "injected_pre.csv", newline="") as f:
        for r in csv.DictReader(f):
            hit = r["hit"] or None
            flag = r["near_bright_core"] == "True"
            review = bool(hit) and not flag and v.get(hit) == "pass"
            out[(run.name.split("_", 1)[1], r["name"])] = {
                "cls": r["cls"],
                "mag": r["mag"],
                "hit": bool(hit),
                "review": review,
                "sign": r["hit_sign"],
            }
    return out


def summarize(argv) -> int:
    """PLAN §4 decision inputs from the stored runs (``runs/<tool>_<name>/``)."""
    import argparse

    import fade_calib as fc
    from astropy.table import Table

    p = argparse.ArgumentParser()
    p.add_argument("summarize")
    p.add_argument("--runs", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args(argv)
    names = ["b1_s30", "b1_s31", "b2_s40"]
    tot, per = {}, {}
    for t in ("arc", "sfft"):
        rows = {}
        for n in names:
            rows.update(_stage_rows(a.runs / f"{t}_{n}"))
        per[t] = rows
        fades = [r for r in rows.values() if r["cls"] in ("partial", "vanish")]
        tot[t] = {
            "unchanged_false_hits": sum(r["hit"] for r in rows.values() if r["cls"] == "unchanged"),
            "unchanged_n": sum(r["cls"] == "unchanged" for r in rows.values()),
            "fade_review_sign": sum(
                r["review"] and fc.sign_ok(r["cls"], r["sign"], {}) for r in fades
            ),
            "fade_n": len(fades),
            "appear_review_sign": sum(
                r["review"] and fc.sign_ok(r["cls"], r["sign"], {})
                for r in rows.values()
                if r["cls"] == "appear"
            ),
            "appear_n": sum(r["cls"] == "appear" for r in rows.values()),
        }
        scr = 0
        for b in ("b1", "b2"):
            d = a.runs / f"{t}_{b}_screen"
            h = Table.read(d / "hits.ecsv")
            v = {
                x["id"]: x["verdict"]
                for x in json.loads((d / "morph_all.json").read_text())["targets"]
            }
            scr += sum((not r["near_bright_core"]) and v.get(r["hit_id"]) == "pass" for r in h)
        tot[t]["screen_review_entries"] = scr
        walls = [
            json.loads((a.runs / f"{t}_{n}" / "rss.json").read_text())["wall_s"]
            for n in [*names, "b1_screen", "b2_screen"]
        ]
        tot[t]["wall_s_total"] = round(sum(walls), 1)
    diffs = []
    for k, ra_ in per["arc"].items():
        rs = per["sfft"][k]
        if ra_["review"] != rs["review"] or ra_["hit"] != rs["hit"]:
            diffs.append(
                {
                    "run": k[0],
                    "name": k[1],
                    "cls": ra_["cls"],
                    "mag": ra_["mag"],
                    "arc": {"hit": ra_["hit"], "review": ra_["review"]},
                    "sfft": {"hit": rs["hit"], "review": rs["review"]},
                }
            )
    ar, sf = tot["arc"], tot["sfft"]
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
        "per_source_differences": diffs,
    }
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({"totals": tot, "criteria": crit, "verdict": verdict, "n_diff": len(diffs)}))
    return 0


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "regions":
        return regions(argv)
    if argv and argv[0] == "summarize":
        return summarize(argv)
    ad.kernel_match = sfft_kernel_match
    return ad.main(argv)


if __name__ == "__main__":
    sys.exit(main())
