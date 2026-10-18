"""Same-source pairing feasibility test (task 1: propose pairs; plan in
``results/same_source/a2744_feasibility/PLAN.md``).

Given a catalogue, every eligible source is a query and every other eligible source at 2-60" a
candidate. Candidates are ranked by the SED match of ``orphan_pairs.match_table`` (free flux scale =
achromatic magnification, 3 % floor, photo-z overlap); nothing else is allowed to differ. Published
multiple images are read only to score the ranking (recall, false-match burden) and to build
matched decoys (similar brightness, colour and size, different source). An optional morphology
step scores pairs with ``parity_test.compare`` (affine map of either handedness).

This proposes pairs; it does not test them against a lens model (task 2: ``lensmodel.find_images``
and ``lens_consistency.py images``). Every number is ``derived``; thresholds are ASSUMPTIONs in the
plan.

    python scripts/same_source_pairs.py --field abell2744 --model abell2744-bergamini23 \\
        --out results/same_source/a2744_feasibility [--morph-image s3://...f200w_i2d.fits]
"""

from __future__ import annotations

import argparse
import json
import os
import time
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import lens_consistency as lc
import numpy as np
import orphan_pairs as op
from scipy.spatial import cKDTree

SEP = (2.0, 60.0)  # arcsec, candidate separations (ASSUMPTION; plan)
MATCH_RADIUS = 0.3  # arcsec, published image -> catalogue row
DECOY_IMAGE_CLEAR = 1.0  # arcsec, decoys this far from every published image
DECOY_TOL = (0.5, 0.15, 0.3)  # |dm| mag, |d colour| dex, |d radius| fraction
TOP_K = (1, 5, 10)


def holdout(system: str) -> bool:
    return zlib.crc32(str(system).encode()) % 2 == 0


def build_sample(cat, images, elig, x, y, ix, iy):
    """Catalogue row (or -1) of each published image: nearest row within MATCH_RADIUS that is
    eligible and matched by no other image."""
    tree = cKDTree(np.column_stack([x, y]))
    d, row = tree.query(np.column_stack([ix, iy]))
    row = np.where(d <= MATCH_RADIUS, row, -1)
    ok = row >= 0
    ok &= np.where(ok, elig[np.clip(row, 0, None)], False)
    vals, counts = np.unique(row[row >= 0], return_counts=True)
    dup = set(vals[counts > 1].tolist())
    ok &= ~np.isin(row, list(dup))
    return np.where(ok, row, -1), {
        "images": len(images),
        "within_0.3": int((d <= MATCH_RADIUS).sum()),
        "eligible_unique": int(ok.sum()),
        "shared_rows": len(dup),
    }


def retrieval(f, e, z16, z84, x, y, elig_idx, queries):
    """``{query_row: (cand_rows, match, p_sed, chi2_nu, sep)}`` over eligible candidates."""
    pts = np.column_stack([x[elig_idx], y[elig_idx]])
    tree = cKDTree(pts)
    out = {}
    for q in queries:
        near = np.asarray(tree.query_ball_point([x[q], y[q]], SEP[1]), int)
        cand = elig_idx[near]
        sep = np.hypot(x[cand] - x[q], y[cand] - y[q])
        keep = (sep >= SEP[0]) & (cand != q)
        cand, sep = cand[keep], sep[keep]
        t = op.match_table(np.full(len(cand), q), cand, f, e, z16, z84)
        out[int(q)] = (
            cand,
            np.asarray(t["match"]),
            np.asarray(t["p_sed"]),
            np.asarray(t["chi2_nu"]),
            sep,
        )
    return out


def rank_of(cand, match, p, chi2_nu, targets):
    order = np.lexsort((chi2_nu, -p, ~match))  # match first, then p desc, then chi2_nu asc
    ranked = cand[order]
    pos = [int(np.flatnonzero(ranked == t)[0]) + 1 for t in targets if t in ranked]
    return min(pos) if pos else None


def pick_decoy(a, b, x, y, mag, col, rad, elig_idx, image_xy):
    sep = np.hypot(x[elig_idx] - x[a], y[elig_idx] - y[a])
    dimg = cKDTree(image_xy).query(np.column_stack([x[elig_idx], y[elig_idx]]))[0]
    dm, dc = mag[elig_idx] - mag[b], col[elig_idx] - col[b]
    dr = rad[elig_idx] / rad[b] - 1.0
    ok = (sep >= SEP[0]) & (sep <= SEP[1]) & (dimg > DECOY_IMAGE_CLEAR) & (elig_idx != a)
    ok &= (np.abs(dm) <= DECOY_TOL[0]) & (np.abs(dc) <= DECOY_TOL[1]) & (np.abs(dr) <= DECOY_TOL[2])
    if not ok.any():
        return None
    cost = (dm / DECOY_TOL[0]) ** 2 + (dc / DECOY_TOL[1]) ** 2 + (dr / DECOY_TOL[2]) ** 2
    return int(elig_idx[ok][np.argmin(np.where(ok, cost, np.inf)[ok])])


def _morph_job(args):
    import parity_test as pt

    tmpl_parts, b_parts = args
    tmpl = pt.Template(*tmpl_parts)
    bx, by, b, bvar = b_parts
    res = pt.compare(tmpl, bx, by, b, bvar, np.eye(2))
    best = min(res["pred"]["chi2"], res["wrong"]["chi2"])
    return {"chi2_nu_best": best / max(b.size, 1), "n_pix": int(b.size)}


def morph_scores(image_uri, ra, dec, pairs):
    """Best-family affine chi^2 / n_pix of B fitted with A's stamp, per (a, b) row pair."""
    import parity_test as pt

    want = sorted({k for p in pairs for k in p})
    tab = {"ra": ra[want], "dec": dec[want]}
    stamps = pt.read_stamps(image_uri, tab, list(range(len(want))), pt.HALO_FILTER_ARCSEC, True)
    pos = {k: i for i, k in enumerate(want)}
    jobs, keys = [], []
    for a, b in pairs:
        sa, sb = stamps.get(pos[a]), stamps.get(pos[b])
        if sa is None or sb is None:
            continue
        sci_a = np.where(sa[6], 0.0, sa[0])
        jobs.append(((sci_a, sa[1], sa[2], sa[3]), pt._fit_pixels(*sb[:4], sb[6])))
        keys.append((a, b))
    with ProcessPoolExecutor(os.cpu_count() or 1) as ex:
        res = list(ex.map(_morph_job, jobs))
    return dict(zip(keys, res, strict=True))


def auc(pos, neg):
    """P(positive score < negative score) + 0.5 ties (lower = more alike)."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if not len(pos) or not len(neg):
        return None
    less = (pos[:, None] < neg[None, :]).mean()
    tie = (pos[:, None] == neg[None, :]).mean()
    return float(less + 0.5 * tie)


def evaluate(rows, half):
    sel = [r for r in rows if r["holdout"] == half]
    ranks = [r["rank"] for r in sel]
    n = len(sel)
    rec = {
        f"recall@{k}": (sum(1 for x in ranks if x is not None and x <= k) / n if n else None)
        for k in TOP_K
    }
    burden = [r["false_matches"] for r in sel]
    return {
        "queries": n,
        **rec,
        "median_rank": float(np.median([x if x is not None else np.inf for x in ranks]))
        if n
        else None,
        "median_false_matches": float(np.median(burden)) if n else None,
        "median_candidates": float(np.median([r["n_cand"] for r in sel])) if n else None,
    }


def verdict(half: dict) -> dict:
    """The plan's pass conditions on one half: P1 recall@10 >= 0.5; P2 median false matches <= 10
    and decoy match rate <= 0.2 and positive match rate >= 0.6 (missing values fail)."""
    r, c = half["retrieval"], half["control"]
    p1 = (r.get("recall@10") or 0.0) >= 0.5
    burden = r.get("median_false_matches")
    dec, pos = c.get("decoy_match_rate"), c.get("positive_match_rate")
    p2 = (
        burden is not None
        and burden <= 10
        and dec is not None
        and dec <= 0.2
        and (pos or 0.0) >= 0.6
    )
    return {"P1": p1, "P2": p2, "pass": p1 and p2}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--field", default="abell2744")
    ap.add_argument("--model", default="abell2744-bergamini23")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--morph-image", default=None, help="i2d (S3 URI) for the morphology add-on")
    ap.add_argument("--morph-half", choices=["dev", "holdout"], default="dev")
    args = ap.parse_args(argv)
    t0 = time.time()
    cat = op.as_standard(op.load_catalogue(args.field))
    images, _ = lc.image_list(args.model, {}, None)
    lc.shift_images(args.model, images)
    elig = op.select_sources(cat, op.FIELDS[args.field]["z_cluster"])
    ra, dec = np.asarray(cat["ra"], float), np.asarray(cat["dec"], float)
    ra0, dec0 = float(np.nanmedian(ra)), float(np.nanmedian(dec))
    x, y = op.tangent_xy(ra, dec, ra0, dec0)
    ix, iy = op.tangent_xy(images["ra"], images["dec"], ra0, dec0)
    row, sample = build_sample(cat, images, elig, x, y, ix, iy)
    f, e = np.asarray(cat["flux"], float), np.asarray(cat["err"], float)
    z16, z84 = np.asarray(cat["z16"], float), np.asarray(cat["z84"], float)
    fs, _ = op.flux_matrix(cat, op.snr_bands_of(cat))
    with np.errstate(invalid="ignore", divide="ignore"):
        mag = -2.5 * np.log10(fs.sum(1)) + 31.4  # nJy -> AB
    col = op.source_colours(cat)
    rad = op.aperture_radius(cat)
    systems = np.asarray(images["system"]).astype(str)
    elig_idx = np.flatnonzero(elig)
    in_s = np.flatnonzero(row >= 0)
    counter = {}  # query row -> rows of the same system (2-60")
    n_close = 0
    pos_pairs = []
    for k in in_s:
        same = [j for j in in_s if j != k and systems[j] == systems[k]]
        tgt = []
        for j in same:
            s = np.hypot(x[row[j]] - x[row[k]], y[row[j]] - y[row[k]])
            if SEP[0] <= s <= SEP[1]:
                tgt.append(int(row[j]))
                if j > k:
                    pos_pairs.append((int(row[k]), int(row[j]), systems[k]))
            elif s < SEP[0] and j > k:
                n_close += 1
        if tgt:
            counter[int(row[k])] = (tgt, systems[k])
    other_images = {int(row[k]): systems[k] for k in in_s}
    res = retrieval(f, e, z16, z84, x, y, elig_idx, list(counter))
    rows = []
    for q, (tgt, sysid) in counter.items():
        cand, match, p, cnu, _ = res[q]
        is_tgt = np.isin(cand, tgt)
        is_img = np.array([c in other_images for c in cand]) & ~is_tgt
        rows.append(
            {
                "query_row": q,
                "system": sysid,
                "holdout": holdout(sysid),
                "n_cand": len(cand),
                "rank": rank_of(cand, match, p, cnu, tgt),
                "target_match": bool(match[is_tgt].any()),
                "false_matches": int((match & ~is_tgt & ~is_img).sum()),
                "other_image_matches": int((match & is_img).sum()),
            }
        )
    image_xy = np.column_stack([ix, iy])
    control = []
    for a, b, sysid in pos_pairs:
        c = pick_decoy(a, b, x, y, mag, col, rad, elig_idx, image_xy)
        tp = op.match_table(np.array([a]), np.array([b]), f, e, z16, z84)
        entry = {
            "a": a,
            "b": b,
            "system": sysid,
            "holdout": holdout(sysid),
            "pos_match": bool(tp["match"][0]),
            "pos_p": float(tp["p_sed"][0]),
            "decoy": c,
        }
        if c is not None:
            td = op.match_table(np.array([a]), np.array([c]), f, e, z16, z84)
            entry.update(decoy_match=bool(td["match"][0]), decoy_p=float(td["p_sed"][0]))
        control.append(entry)
    t_ret = time.time() - t0

    def ctrl(half):
        sel = [c for c in control if c["holdout"] == half and c["decoy"] is not None]
        return {
            "pairs": sum(1 for c in control if c["holdout"] == half),
            "with_decoy": len(sel),
            "positive_match_rate": float(np.mean([c["pos_match"] for c in sel])) if sel else None,
            "decoy_match_rate": float(np.mean([c["decoy_match"] for c in sel])) if sel else None,
        }

    summary = {
        "plan": "PLAN.md",
        "catalogue_rows": len(cat),
        "eligible_rows": int(elig.sum()),
        "sample": sample,
        "positive_pairs_2_60": len(pos_pairs),
        "positive_pairs_below_2": n_close,
        "dev": {"retrieval": evaluate(rows, False), "control": ctrl(False)},
        "holdout": {"retrieval": evaluate(rows, True), "control": ctrl(True)},
        "cost_s": {"retrieval_and_control": round(t_ret, 1)},
    }
    for half in ("dev", "holdout"):
        summary[half]["verdict"] = verdict(summary[half])
    if args.morph_image:
        t1 = time.time()
        want = args.morph_half == "holdout"
        sel = [c for c in control if c["holdout"] == want and c["decoy"] is not None]
        pairs = [(c["a"], c["b"]) for c in sel] + [(c["a"], c["decoy"]) for c in sel]
        sc = morph_scores(args.morph_image, ra, dec, pairs)
        for c in sel:
            c["morph_pos"] = sc.get((c["a"], c["b"]), {}).get("chi2_nu_best")
            c["morph_decoy"] = sc.get((c["a"], c["decoy"]), {}).get("chi2_nu_best")
        ok = [c for c in sel if c["morph_pos"] is not None and c["morph_decoy"] is not None]
        summary[args.morph_half]["morph"] = {
            "pairs_scored": len(ok),
            "auc": auc([c["morph_pos"] for c in ok], [c["morph_decoy"] for c in ok]),
            "paired_win_rate": float(np.mean([c["morph_pos"] < c["morph_decoy"] for c in ok]))
            if ok
            else None,
        }
        summary["cost_s"][f"morph_{args.morph_half}"] = round(time.time() - t1, 1)
    args.out.mkdir(parents=True, exist_ok=True)
    tag = f"_morph_{args.morph_half}" if args.morph_image else ""
    (args.out / f"queries{tag}.json").write_text(json.dumps(rows, indent=1) + "\n")
    (args.out / f"control{tag}.json").write_text(json.dumps(control, indent=1, default=str) + "\n")
    (args.out / f"summary{tag}.json").write_text(json.dumps(op._finite(summary), indent=1) + "\n")
    print(json.dumps(op._finite(summary), indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
