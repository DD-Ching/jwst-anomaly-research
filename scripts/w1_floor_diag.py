"""What sets the W1 shear-screen B-mode floor? Per-source contributions at the |S_x| maximum.

Diagnostic for TASKS "JWST focus" 1 (W1): at a field's largest |S_x| centre, list the sources that
carry S_x (their share of the sum), the aperture's mean ellipticity (a coherent PSF-anisotropy
residual shows as a mean e that survives at many centres) and whether the top contributors are
collinear chains (streak segments: major axes aligned with each other and with the line joining
them). All outputs are ``derived``; the chain tolerances are ASSUMPTIONs.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import inject_shear as ish  # noqa: E402

es = ish.es
CHAIN_PA_DEG = 10.0  # ASSUMPTION: neighbour's major axis within this of the source's
CHAIN_LINE_DEG = 10.0  # ASSUMPTION: direction to the neighbour within this of the source's axis
CHAIN_RADIUS = 20.0  # arcsec, ASSUMPTION
NEIGHBOUR_DMAG = 2.0  # ASSUMPTION: a neighbour at least this much brighter can bias a shape
NEIGHBOUR_K = 5.0  # ASSUMPTION: within this many of the neighbour's semimajor sigmas


def bright_neighbour(x, y, mag, a_arcsec, dmag=NEIGHBOUR_DMAG, k=NEIGHBOUR_K) -> np.ndarray:
    """True where a source lies within ``k`` semimajor sigmas of a neighbour >= ``dmag`` brighter
    (positions in arcsec; the neighbour's size sets the reach of its light)."""
    from scipy.spatial import cKDTree

    x, y, mag, a = (np.asarray(v, float) for v in (x, y, mag, a_arcsec))
    out = np.zeros(len(x), bool)
    ok = np.isfinite(mag) & np.isfinite(a) & np.isfinite(x) & np.isfinite(y)
    tree = cKDTree(np.c_[x[ok], y[ok]])
    idx_ok = np.flatnonzero(ok)
    for b in idx_ok:
        for j in tree.query_ball_point([x[b], y[b]], k * a[b]):
            i = idx_ok[j]
            if i != b and mag[i] - mag[b] >= dmag:
                out[i] = True
    return out


def chain_neighbours(x, y, pa, radius=CHAIN_RADIUS, pa_tol=CHAIN_PA_DEG, line_tol=CHAIN_LINE_DEG):
    """Number of other sources within ``radius`` whose major axis and connecting line both lie
    within the tolerances of each source's major axis (model frame: x West, y North; PA E of N)."""
    from scipy.spatial import cKDTree

    x, y, pa = (np.asarray(v, float) for v in (x, y, pa))
    hits = cKDTree(np.c_[x, y]).query_ball_point(np.c_[x, y], radius)
    out = np.zeros(len(x), int)
    for i, h in enumerate(hits):
        h = np.asarray([j for j in h if j != i], int)
        if not len(h):
            continue
        line = np.mod(np.degrees(np.arctan2(-(x[h] - x[i]), y[h] - y[i])), 180.0)
        ok = (es.lensmodel.axis_offset_deg(pa[h], np.full(len(h), pa[i])) <= pa_tol) & (
            es.lensmodel.axis_offset_deg(line, np.full(len(h), pa[i])) <= line_tol
        )
        out[i] = int(ok.sum())
    return out


def diagnose(name: str, out: Path, top: int, dmag: float, k_nb: float) -> dict:
    spec = ish.field_spec(name)
    (out / name).mkdir(parents=True, exist_ok=True)
    model, shapes, extra, cat, _ = ish.ir.load_field(name, out / name, spec)
    sargs = SimpleNamespace(**ish.SHEAR_DEFAULTS)
    sargs.max_radius = spec["max_radius"]
    psf = es.psf_sigma_px(shapes)
    res = es.shear_screen(model, shapes, sargs, psf, extra)
    scale = ish.ir.pixel_scale_arcsec(cat)
    xa, ya = model.to_frame(np.asarray(shapes["ra"]), np.asarray(shapes["dec"]))
    nb = bright_neighbour(
        xa, ya, shapes["mag"], np.asarray(shapes["semimajor_px"], float) * scale, dmag, k_nb
    )
    # exploratory: the same screen without bright-neighbour rows (no injections, no adoption)
    alt = es.shear_screen(model, shapes[~nb], sargs, psf, extra)
    sx, s = res["s_cross"], res["s"]
    k = int(np.nanargmax(np.abs(sx)))
    cx, cy = float(res["gx"].ravel()[k]), float(res["gy"].ravel()[k])
    ap = res["ap"]
    e = res["e"][res["use"]]
    x, y = res["xy"]
    row = ap.w.getrow(k)
    cols, w = row.indices, row.data
    proj = w * e[cols]
    w2 = ap.w2[k, cols].toarray().ravel()
    noise = float(np.sqrt(w2 @ (np.abs(e[cols]) ** 2) / 2.0))
    contrib = -proj.imag / noise  # per-source share of S_x at this centre
    total = float(contrib.sum())
    if not np.isclose(total, float(sx.ravel()[k])):
        raise RuntimeError(f"per-source S_x sum {total} != map value {sx.ravel()[k]}")
    order = np.argsort(-np.sign(total) * contrib)
    pa = np.mod(0.5 * np.degrees(np.angle(e)), 180.0)
    chains = chain_neighbours(x, y, pa)
    q = np.abs(w)
    mean_e = complex(np.sum(q * e[cols]) / np.sum(q))
    idx_use = np.flatnonzero(res["use"])
    rows = []
    for j in order[:top]:
        c = int(cols[j])
        r = shapes[idx_use[c]]
        rows.append(
            {
                "row": int(idx_use[c]),
                "ra": float(r["ra"]),
                "dec": float(r["dec"]),
                "x": float(x[c]),
                "y": float(y[c]),
                "mag": float(r["mag"]),
                "snr": float(r["snr"]),
                "abs_e": float(abs(e[c])),
                "pa_corr_deg": float(pa[c]),
                "contrib": float(contrib[j]),
                "chain_neighbours": int(chains[c]),
                "bright_neighbour": bool(nb[idx_use[c]]),
            }
        )
    nb_c = nb[idx_use[cols]]
    srt = np.sort(np.sign(total) * contrib)[::-1]
    return {
        "field": name,
        "provenance": "derived",
        "centre_xy": [cx, cy],
        "s_cross": float(sx.ravel()[k]),
        "s": float(s.ravel()[k]),
        "n_sources": int(len(cols)),
        "top1_share": float(srt[0] / abs(total)),
        "top3_share": float(srt[:3].sum() / abs(total)),
        "top5_share": float(srt[:5].sum() / abs(total)),
        "aperture_mean_e": [mean_e.real, mean_e.imag],
        "aperture_mean_abs_e": float(abs(mean_e)),
        "field_mean_abs_e": float(np.mean(np.abs(e))),
        "field_frac_chain_ge2": float(np.mean(chains >= 2)),
        "aperture_frac_chain_ge2": float(np.mean(chains[cols] >= 2)),
        "share_from_bright_neighbour_rows": float(contrib[nb_c].sum() / total),
        "aperture_frac_bright_neighbour": float(nb_c.mean()),
        "field_frac_bright_neighbour_used": float(nb[idx_use].mean()),
        "without_bright_neighbour": {
            "n_used": int(alt["use"].sum()),
            "n_used_all": int(res["use"].sum()),
            "b_max": float(np.nanmax(np.abs(alt["s_cross"]))),
            "b_max_all": float(np.nanmax(np.abs(sx))),
            "s_max": float(np.nanmax(alt["s"])),
            "s_max_all": float(np.nanmax(s)),
            "responsivity": float(alt["counts"]["responsivity"]),
            "responsivity_all": float(res["counts"]["responsivity"]),
            "note": "exploratory screen only: no rotation-null p, no injections, not adopted",
        },
        "neighbour_rule": {"dmag": dmag, "k_semimajor_sigma": k_nb},
        "top": rows,
        "chain_tolerances": {
            "pa_deg": CHAIN_PA_DEG,
            "line_deg": CHAIN_LINE_DEG,
            "radius_arcsec": CHAIN_RADIUS,
        },
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fields", nargs="+", default=list(ish.DEFAULT_FIELDS))
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--dmag", type=float, nargs="+", default=[NEIGHBOUR_DMAG])
    ap.add_argument("--k", type=float, nargs="+", default=[NEIGHBOUR_K])
    a = ap.parse_args(argv)
    res = [diagnose(f, a.work, a.top, dm, k) for f in a.fields for dm in a.dmag for k in a.k]
    a.out.write_text(json.dumps(es.json_safe(res), indent=1) + "\n")
    for r in res:
        print(
            r["field"],
            r["neighbour_rule"],
            [round(v, 1) for v in r["centre_xy"]],
            f"Sx {r['s_cross']:.2f} n {r['n_sources']} top1 {r['top1_share']:.2f} "
            f"top3 {r['top3_share']:.2f} top5 {r['top5_share']:.2f} "
            f"<e> {r['aperture_mean_abs_e']:.3f} field<|e|> {r['field_mean_abs_e']:.3f} "
            f"chain>=2 ap {r['aperture_frac_chain_ge2']:.3f} field {r['field_frac_chain_ge2']:.3f}",
        )
        print(
            "   ",
            r["without_bright_neighbour"],
            "share",
            round(r["share_from_bright_neighbour_rows"], 2),
        )
        for t in r["top"][:6]:
            print("   ", {k: round(v, 3) if isinstance(v, float) else v for k, v in t.items()})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
