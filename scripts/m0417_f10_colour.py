"""MACS0417 family 10: catalogue colour test of CANUCS 1115938 against 10.1 / 10.2.

Lensing is achromatic, so a third image of family 10 must match 10.1 / 10.2 in colour up to a flux
ratio. Pair chi^2 (flux ratio free) on CANUCS DR1 0.3" colour fluxes with a 5 % error floor; the
test is validated first on the first two images of every other Jauzac+19 family.
Frozen plan: ``results/lens_consistency/macs0417_f10_colour/PLAN.md``.

    python scripts/m0417_f10_colour.py --cat <CANUCS DR1 photometry-cat.fits.gz> [--posthoc]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import m0417_frame as mf
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from scipy import optimize, stats

from jwst_anomaly.provenance import file_sha256

OUT = Path("results/lens_consistency/macs0417_f10_colour")
FRAME = Path("results/lens_consistency/macs0417_frame")
THIRD = 1115938  # CANUCS DR1 SOURCE id (``../macs0417_f10_zscan``)
FAMILY = "10"
# Frozen rules (PLAN.md; ASSUMPTIONs)
BANDS = ("F090W", "F115W", "F150W", "F200W", "F277W", "F356W", "F410M", "F444W")
MATCH_ARCSEC = 0.3
ERR_FLOOR = 0.05
P_DIFFER = 0.01
MIN_BANDS = 5
SNR_KNOWN = 10.0
GATE_MIN_PAIRS = 8
GATE_MIN_FRAC = 0.8


def pair_chi2(fa, ea, fb, eb) -> tuple[float, int, float]:
    """min over r of sum (fa - r fb)^2 / (sa^2 + r^2 sb^2) with the error floor; returns
    (chi2, dof, r). Bands with a non-finite flux or error, or an error <= 0, are dropped."""
    fa, ea, fb, eb = (np.asarray(v, float) for v in (fa, ea, fb, eb))
    ok = np.isfinite(fa) & np.isfinite(fb) & np.isfinite(ea) & np.isfinite(eb) & (ea > 0) & (eb > 0)
    fa, ea, fb, eb = fa[ok], ea[ok], fb[ok], eb[ok]
    sa2 = ea**2 + (ERR_FLOOR * fa) ** 2
    sb2 = eb**2 + (ERR_FLOOR * fb) ** 2
    if fa.size == 0:
        return float("nan"), -1, float("nan")

    def chi(r):
        return float(np.sum((fa - r * fb) ** 2 / (sa2 + r * r * sb2)))

    r0 = float(np.sum(fa * fb / sa2) / np.sum(fb * fb / sa2)) if np.any(fb) else 1.0
    lo, hi = sorted((r0 / 10, r0 * 10)) if r0 != 0 else (-1.0, 1.0)
    res = optimize.minimize_scalar(chi, bounds=(lo, hi), method="bounded")
    return float(res.fun), int(fa.size - 1), float(res.x)


def classify(chi2: float, dof: int) -> str:
    """Pair class: "untestable" below MIN_BANDS bands, else "differ" (p < P_DIFFER) or
    "consistent"."""
    if dof + 1 < MIN_BANDS or not np.isfinite(chi2):
        return "untestable"
    return "differ" if stats.chi2.sf(chi2, dof) < P_DIFFER else "consistent"


def gate(classes: list[str]) -> dict:
    testable = [c for c in classes if c != "untestable"]
    ok = sum(c == "consistent" for c in testable)
    frac = ok / len(testable) if testable else None
    return {
        "n_testable": len(testable),
        "n_consistent": ok,
        "fraction": frac,
        "passes": len(testable) >= GATE_MIN_PAIRS and frac is not None and frac >= GATE_MIN_FRAC,
    }


def decide(gate_passes: bool, target_classes: list[str]) -> str:
    if not gate_passes:
        return "not validated"
    if "differ" in target_classes:
        return "H2 favoured: colours differ"
    if target_classes and all(c == "consistent" for c in target_classes):
        return "H1 not excluded"
    return "untestable"


def fluxes(row) -> tuple[np.ndarray, np.ndarray]:
    f = np.array([float(row[f"FLUX_COLOR03_TOTAL_{b}"]) for b in BANDS])
    e = np.array([float(row[f"FLUXERR_COLOR03_TOTAL_{b}"]) for b in BANDS])
    return f, e


def counterparts(images: Table, cat: Table) -> dict[str, int]:
    """Image id -> CANUCS row index of the nearest source within MATCH_ARCSEC."""
    ci = SkyCoord(np.asarray(images["ra"], float), np.asarray(images["dec"], float), unit="deg")
    cc = SkyCoord(np.asarray(cat["RA"], float), np.asarray(cat["DEC"], float), unit="deg")
    idx, sep, _ = ci.match_to_catalog_sky(cc)
    return {
        str(i): int(k)
        for i, k, s in zip(images["image_id"], idx, sep.arcsec, strict=True)
        if s <= MATCH_ARCSEC
    }


def test_pair(cat: Table, ka: int, kb: int, label: str) -> dict:
    fa, ea = fluxes(cat[ka])
    fb, eb = fluxes(cat[kb])
    chi2, dof, r = pair_chi2(fa, ea, fb, eb)
    cls = classify(chi2, dof)
    return {
        "pair": label,
        "sources": [int(cat["SOURCE"][ka]), int(cat["SOURCE"][kb])],
        "chi2": None if not np.isfinite(chi2) else round(chi2, 3),
        "dof": dof,
        "p": None if dof < 1 or not np.isfinite(chi2) else float(stats.chi2.sf(chi2, dof)),
        "ratio": None if not np.isfinite(r) else round(r, 4),
        "class": cls,
    }


def specificity(
    cat: Table, refs: dict[str, int], flux_lo: float, flux_hi: float, exclude: int
) -> dict:
    """Post hoc (not a frozen rule): of the catalogue sources with F200W S/N >= SNR_KNOWN and
    F200W flux in [flux_lo, flux_hi], other than the references and ``exclude`` (the target),
    the fraction "consistent" with every reference image."""
    f = np.asarray(cat["FLUX_COLOR03_TOTAL_F200W"], float)
    e = np.asarray(cat["FLUXERR_COLOR03_TOTAL_F200W"], float)
    with np.errstate(divide="ignore", invalid="ignore"):
        sel = np.flatnonzero((f / e >= SNR_KNOWN) & (f >= flux_lo) & (f <= flux_hi))
    sel = [k for k in sel if k not in refs.values() and k != exclude]
    both = 0
    for k in sel:
        if all(test_pair(cat, r, k, "")["class"] == "consistent" for r in refs.values()):
            both += 1
    return {
        "label": "post hoc",
        "flux_range_njy": [flux_lo, flux_hi],  # catalogue fluxes are in nJy
        "n_sources": len(sel),
        "n_consistent_with_all_refs": both,
        "fraction": both / len(sel) if sel else None,
    }


def main(argv: list[str] | None = None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cat", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--posthoc", action="store_true", help="also write posthoc.json (specificity)")
    args = ap.parse_args(argv)
    cat = Table.read(args.cat)
    off = json.loads((FRAME / "summary.json").read_text())["measurement"]["offset"]
    images = mf.shifted_loader(off["dra"], off["ddec"])("macs0417")
    cp = counterparts(images, cat)
    sys_of = dict(zip(map(str, images["image_id"]), map(str, images["system"]), strict=True))
    known = []
    for s in sorted(set(sys_of.values()) - {FAMILY}):
        ids = sorted(i for i in cp if sys_of[i] == s)
        if len(ids) < 2:
            continue
        ka, kb = cp[ids[0]], cp[ids[1]]
        if ka == kb:
            continue  # both images matched one source: no pair to test
        snr = [
            float(cat["FLUX_COLOR03_TOTAL_F200W"][k] / cat["FLUXERR_COLOR03_TOTAL_F200W"][k])
            for k in (ka, kb)
        ]
        if min(snr) < SNR_KNOWN:
            continue
        known.append(test_pair(cat, ka, kb, f"{ids[0]}-{ids[1]}"))
    g = gate([k["class"] for k in known])
    k3 = int(np.flatnonzero(np.asarray(cat["SOURCE"]) == THIRD)[0])
    targets = [test_pair(cat, cp[i], k3, f"{i}-{THIRD}") for i in ("10.1", "10.2") if i in cp]
    ref = (
        test_pair(cat, cp["10.1"], cp["10.2"], "10.1-10.2") if {"10.1", "10.2"} <= set(cp) else None
    )
    out = {
        "catalogue_sha256": file_sha256(args.cat),
        "rules": {
            "bands": list(BANDS),
            "match_arcsec": MATCH_ARCSEC,
            "err_floor": ERR_FLOOR,
            "p_differ": P_DIFFER,
            "min_bands": MIN_BANDS,
            "snr_known": SNR_KNOWN,
            "gate": [GATE_MIN_PAIRS, GATE_MIN_FRAC],
        },
        "n_counterparts": len(cp),
        "known_pairs": known,
        "gate": g,
        "targets": targets,
        "reference_10_1_10_2": ref,
        "decision": decide(g["passes"], [t["class"] for t in targets]),
        "provenance": "observed fluxes (CANUCS DR1); derived chi2",
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "known_pairs"}, indent=1))
    if args.posthoc:
        f3 = float(cat["FLUX_COLOR03_TOTAL_F200W"][k3])
        refs = {i: cp[i] for i in ("10.1", "10.2") if i in cp}
        ph = specificity(cat, refs, f3 / 2, f3 * 2, k3)
        (args.out / "posthoc.json").write_text(json.dumps(ph, indent=1, allow_nan=False) + "\n")
        print(json.dumps(ph, indent=1))
    return out


if __name__ == "__main__":
    main()
