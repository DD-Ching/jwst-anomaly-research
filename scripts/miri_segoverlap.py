"""Segment-overlap MIRI band association on SMACS 0723 (TASKS Then 7, D-014 addendum 6).

Plan: ``results/miri_confirm/PLAN_segoverlap.md`` (rule frozen before measuring). Reads the four
MAST L3 MIRI ``_segm.fits`` and ``_cat.ecsv`` of program 2736 o002 and writes
``summary_segoverlap.json``.

Provenance: segmentation maps and catalogs are ``observed`` (pipeline products); resampled maps,
overlaps, pair counts and chance fractions are ``derived``; the f grid, the 0.1" radius, the
shifts and the 0.10 chance ceiling are ASSUMPTIONs from the plan.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from jwst_anomaly import catalog

PREFIX = "jw02736-o002_t001_miri_"
REF = "F770W"
OTHERS = ("F1000W", "F1500W", "F1800W")
F_GRID = (0.1, 0.25, 0.5, 0.75)  # ASSUMPTION (plan)
R0 = 0.1  # arcsec, ASSUMPTION (plan): current merge radius
SHIFTS = ((8.0, 8.0), (8.0, -8.0), (-8.0, 8.0), (-8.0, -8.0))  # arcsec, as miri_sizemerge.py
BAND_ROTATION_DEG = 120.0  # per non-reference band, as miri_sizemerge.py
MAX_CHANCE = 0.10  # ASSUMPTION (plan)
VETTED = {221: True, 239: True, 118: True, 240: True, 176: False}  # F770W label -> merge expected


def band_shift(i: int, dx: float, dy: float) -> tuple[float, float]:
    """Shift (arcsec) of the i-th non-reference band: base vector rotated by 120 deg per band."""
    t = np.deg2rad(BAND_ROTATION_DEG * i)
    return dx * np.cos(t) - dy * np.sin(t), dx * np.sin(t) + dy * np.cos(t)


def ref_sky(wcs_ref: WCS, shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    """RA, Dec (deg) of every reference-grid pixel centre."""
    yy, xx = np.indices(shape)
    return wcs_ref.pixel_to_world_values(xx, yy)


def resample_labels(
    segm: np.ndarray,
    wcs_b: WCS,
    ra: np.ndarray,
    dec: np.ndarray,
    shift: tuple[float, float] = (0.0, 0.0),
) -> np.ndarray:
    """Nearest-neighbour band-b labels at the reference pixels' sky positions.

    ``shift`` (arcsec, east and north) moves the band-b map on the sky, so the reference pixel at
    sky position s reads band b at s - shift. Outside the band-b map the label is 0.
    """
    cosd = np.cos(np.deg2rad(np.nanmedian(dec)))
    x, y = wcs_b.world_to_pixel_values(ra - shift[0] / 3600.0 / cosd, dec - shift[1] / 3600.0)
    xi = np.rint(np.nan_to_num(x, nan=-1)).astype(np.int64)
    yi = np.rint(np.nan_to_num(y, nan=-1)).astype(np.int64)
    ok = (xi >= 0) & (yi >= 0) & (xi < segm.shape[1]) & (yi < segm.shape[0])
    out = np.zeros(ra.shape, dtype=np.int64)
    out[ok] = segm[yi[ok], xi[ok]]
    return out


def overlap_pairs(ref_seg: np.ndarray, b_seg: np.ndarray, area_b: np.ndarray | None = None):
    """(i, j, o_ij, n_ij) for every pair of labels sharing pixels; o = n_ij / min(A_i, A_j).

    Areas are pixel counts on the common (reference) grid, or ``area_b[j]`` for band b when given
    (native-grid areas, unaffected by truncation at the reference grid's edge; post hoc).
    """
    a = ref_seg.ravel()
    b = b_seg.ravel()
    both = (a > 0) & (b > 0)
    if not both.any():
        return np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0), np.zeros(0, np.int64)
    pairs, n = np.unique(np.stack([a[both], b[both]]), axis=1, return_counts=True)
    area_a = np.bincount(a, minlength=int(a.max()) + 1)
    if area_b is None:
        area_b = np.bincount(b, minlength=int(b.max()) + 1)
    i, j = pairs
    o = n / np.minimum(area_a[i], area_b[j])
    return i, j, o, n


def greedy_one_to_one(
    i: np.ndarray,
    j: np.ndarray,
    score: np.ndarray,
    threshold: float,
    tiebreak: np.ndarray | None = None,
):
    """Accept pairs with score >= threshold, one-to-one, greedy by decreasing score.

    Ties are broken by decreasing ``tiebreak`` when given (shared pixels; post hoc), then by
    (i, j), so the result is deterministic.
    """
    keep = score >= threshold
    tb = np.zeros(len(i)) if tiebreak is None else np.asarray(tiebreak, float)
    i, j, score, tb = i[keep], j[keep], score[keep], tb[keep]
    order = np.lexsort((j, i, -tb, -score))
    used_i: set[int] = set()
    used_j: set[int] = set()
    out = []
    for k in order:
        a, b = int(i[k]), int(j[k])
        if a in used_i or b in used_j:
            continue
        used_i.add(a)
        used_j.add(b)
        out.append((a, b))
    return out


def chance_fraction(extra_real: float, extra_chance: float) -> float | None:
    """Chance share of the extra pairs; None (never passes) when there are no real extras."""
    if extra_real <= 0:
        return None
    return extra_chance / extra_real


def adopt_f(fs, fractions, max_chance: float = MAX_CHANCE) -> float | None:
    """Frozen rule: the smallest f whose chance fraction <= ``max_chance``; None if none passes."""
    ok = [f for f, c in zip(fs, fractions, strict=True) if c is not None and c <= max_chance]
    return min(ok, default=None)


def sep_arcsec(ra1, dec1, ra2, dec2) -> np.ndarray:
    """Small-angle separation (arcsec)."""
    cosd = np.cos(np.deg2rad(dec1))
    return 3600.0 * np.hypot((ra1 - ra2) * cosd, dec1 - dec2)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("data_dir", type=Path, help="directory with the four MIRI _segm.fits/_cat.ecsv")
    p.add_argument("--out", type=Path, default=Path("results/miri_confirm/summary_segoverlap.json"))
    a = p.parse_args(argv)
    seg, wcs, pos, meta = {}, {}, {}, {}
    for b in (REF, *OTHERS):
        fp = a.data_dir / f"{PREFIX}{b.lower()}_segm.fits"
        with fits.open(fp) as h:
            seg[b] = np.asarray(h["SCI"].data, dtype=np.int64)
            wcs[b] = WCS(h["SCI"].header)
        cat = catalog.load_pipeline_catalog(a.data_dir / f"{PREFIX}{b.lower()}_cat.ecsv")
        pos[b] = {
            int(lab): (float(r), float(d))
            for lab, r, d in zip(cat["label"], cat["ra"], cat["dec"], strict=True)
        }
        meta[b] = {
            "n_catalog": len(cat),
            "jwst_version": cat.meta["jwst_version"],
            "cat_sha256": cat.meta["input_sha256"],
            "segm_sha256": _sha(fp),
        }
    ra, dec = ref_sky(wcs[REF], seg[REF].shape)

    def extra_count(pairs, b, shift=(0.0, 0.0)) -> tuple[int, int]:
        """(pairs with centroids > R0 apart, pairs within R0) after shifting band b."""
        cosd = np.cos(np.deg2rad(np.nanmedian(dec)))
        n_far = n_near = 0
        for i, j in pairs:
            if i not in pos[REF] or j not in pos[b]:
                n_far += 1  # no catalog row: cannot be a 0.1" pair
                continue
            r1, d1 = pos[REF][i]
            r2, d2 = pos[b][j]
            s = sep_arcsec(r1, d1, r2 + shift[0] / 3600.0 / cosd, d2 + shift[1] / 3600.0)
            if s > R0:
                n_far += 1
            else:
                n_near += 1
        return n_far, n_near

    real_ov, chance_ov, real_nat, chance_nat = {}, {}, {}, {}
    for k, b in enumerate(OTHERS):
        real_ov[b] = overlap_pairs(seg[REF], resample_labels(seg[b], wcs[b], ra, dec))
        chance_ov[b] = []
        area_nat = np.bincount(seg[b].ravel())
        real_nat[b] = overlap_pairs(seg[REF], resample_labels(seg[b], wcs[b], ra, dec), area_nat)
        chance_nat[b] = []
        for dx, dy in SHIFTS:
            sh = band_shift(k, dx, dy)
            moved = resample_labels(seg[b], wcs[b], ra, dec, sh)
            chance_ov[b].append((sh, overlap_pairs(seg[REF], moved)))
            chance_nat[b].append((sh, overlap_pairs(seg[REF], moved, area_nat)))

    # Validation on known cases: every 0.1" centroid pair (not one-to-one) should overlap.
    validation = {}
    for b in OTHERS:
        i, j, o, _ = real_ov[b]
        best = {}
        for a_, b_, o_ in zip(i, j, o, strict=True):
            best[(int(a_), int(b_))] = float(o_)
        near = []
        for li, (r1, d1) in pos[REF].items():
            for lj, (r2, d2) in pos[b].items():
                if abs(d1 - d2) * 3600 <= R0 and sep_arcsec(r1, d1, r2, d2) <= R0:
                    near.append(best.get((li, lj), 0.0))
        near = np.asarray(near)
        validation[b] = {
            "n_pairs_within_r0": int(near.size),
            "overlap_quantiles_10_50_90": [
                round(float(q), 3) for q in np.quantile(near, [0.1, 0.5, 0.9])
            ]
            if near.size
            else None,
            "n_overlap_ge_f": {str(f): int(np.sum(near >= f)) for f in F_GRID},
        }

    rows, fracs = [], []
    for f in F_GRID:
        real_far = real_near = 0
        chance_far = 0.0
        vet_bands = {k: 0 for k in VETTED}
        per_band = {}
        for b in OTHERS:
            acc = greedy_one_to_one(*real_ov[b][:3], f)
            far, near = extra_count(acc, b)
            real_far += far
            real_near += near
            for li, _ in acc:
                if li in vet_bands:
                    vet_bands[li] += 1
            cf = [extra_count(greedy_one_to_one(*ov[:3], f), b, sh)[0] for sh, ov in chance_ov[b]]
            chance_far += float(np.mean(cf))
            per_band[b] = {
                "accepted": len(acc),
                "extra_real": far,
                "extra_chance_mean": float(np.mean(cf)),
            }
        frac = chance_fraction(real_far, chance_far)
        fracs.append(frac)
        rows.append(
            {
                "f": f,
                "accepted_within_r0": real_near,
                "extra_real": real_far,
                "extra_chance_mean": chance_far,
                "chance_fraction": None if frac is None else round(frac, 4),
                "passes": frac is not None and frac <= MAX_CHANCE,
                "per_band": per_band,
                "vetted_bands_gained": {f"f770w_{k}": v for k, v in vet_bands.items()},
            }
        )
    adopted = adopt_f(F_GRID, fracs)

    # Post hoc (code review): the frozen counting above (a) breaks overlap ties by label, (b) lets
    # chance pairs use labels that the real run spends on 0.1" pairs, (c) measures band-b areas
    # after truncation at the F770W grid edge and (d) counts labels without a catalog row as
    # extra. This variant fixes all four: native band-b areas, ties by shared pixels, pairs only
    # between catalogued labels, and the labels of real 0.1" pairs removed from the chance run.
    corrected, cfracs = [], []
    for f in F_GRID:
        real_far = chance_far = 0.0
        vet_bands = {k: 0 for k in VETTED}
        for b in OTHERS:

            def accept(ov, drop_i=frozenset(), drop_j=frozenset(), b=b, f=f):
                i, j, o, n = ov
                keep = np.array(
                    [
                        int(x) in pos[REF]
                        and int(y) in pos[b]
                        and x not in drop_i
                        and y not in drop_j
                        for x, y in zip(i, j, strict=True)
                    ],
                    dtype=bool,
                )
                return greedy_one_to_one(i[keep], j[keep], o[keep], f, n[keep])

            acc = accept(real_nat[b])
            near_pairs = [(x, y) for x, y in acc if extra_count([(x, y)], b)[1]]
            real_far += len(acc) - len(near_pairs)
            for li, _ in acc:
                if li in vet_bands:
                    vet_bands[li] += 1
            di = frozenset(x for x, _ in near_pairs)
            dj = frozenset(y for _, y in near_pairs)
            cf = [extra_count(accept(ov, di, dj), b, sh)[0] for sh, ov in chance_nat[b]]
            chance_far += float(np.mean(cf))
        frac = chance_fraction(real_far, chance_far)
        cfracs.append(frac)
        corrected.append(
            {
                "f": f,
                "extra_real": int(real_far),
                "extra_chance_mean": chance_far,
                "chance_fraction": None if frac is None else round(frac, 4),
                "passes": frac is not None and frac <= MAX_CHANCE,
                "vetted_bands_gained": {f"f770w_{k}": v for k, v in vet_bands.items()},
            }
        )
    if adopted is None:
        hits = None
    else:
        vet = rows[F_GRID.index(adopted)]["vetted_bands_gained"]
        hits = sum(1 for k, exp in VETTED.items() if exp and vet[f"f770w_{k}"] >= 1)
    out = {
        "plan": "results/miri_confirm/PLAN_segoverlap.md",
        "inputs": meta,
        "ref_band": REF,
        "r0_arcsec": R0,
        "shifts_arcsec": SHIFTS,
        "band_rotation_deg": BAND_ROTATION_DEG,
        "max_chance_fraction": MAX_CHANCE,
        "validation_r0_pairs": validation,
        "validation_n_pairs_within_r0_total": sum(
            v["n_pairs_within_r0"] for v in validation.values()
        ),
        "grid": rows,
        "adopted_f": adopted,
        "post_hoc_corrected_grid": corrected,
        "post_hoc_corrected_adopted_f": adopt_f(F_GRID, cfracs),
        "vetted_expected_merges_hit": None if hits is None else f"{hits} / {sum(VETTED.values())}",
        "vetted_176_single": None
        if adopted is None
        else rows[F_GRID.index(adopted)]["vetted_bands_gained"]["f770w_176"] == 0,
    }
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(validation))
    for r in rows:
        print(
            {
                k: r[k]
                for k in (
                    "f",
                    "accepted_within_r0",
                    "extra_real",
                    "extra_chance_mean",
                    "chance_fraction",
                )
            },
            r["vetted_bands_gained"],
        )
    for r in corrected:
        print("post hoc", r)
    print("adopted", adopted, "vetted hits", out["vetted_expected_merges_hit"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
