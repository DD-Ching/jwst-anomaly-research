"""Post hoc clump checks for Abell S1063 ``arcs`` anti rows (vetting step 3; not in PLAN_1557.md).

1. ``third``: 1557's catalogue segment holds a third, PSF-sized component along its axis that the
   two-component fit of ``arcs_vet_deblend.py`` does not model. Fit it (round Gaussian at the
   brightest 1″-filtered F200W pixel 0.25–0.6″ from the knot pair), then fit fixed-position
   fluxes of all three in F090W–F200W and compare flux ratios (``ratio_consistency``).
2. ``clump_anti_fraction``: the fraction of randomly oriented point-pair separations that a lens
   with reduced shear ``g`` (κ < 1) maps ≥ 60° from the shear axis (model_prediction; ASSUMPTIONs:
   isotropic separations, point-like clumps, i.e. the largest possible intrinsic elongation).

    python scripts/arcs_vet_clumps.py --cut-dir <dir with the 1557 cutouts of arcs_vet_deblend.py>

Output: ``results/arcs/abells1063-cats/vetting/clumps_1557.json`` (derived).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.optimize import least_squares

from jwst_anomaly import paths

_spec = importlib.util.spec_from_file_location(
    "arcs_vet_deblend", Path(__file__).resolve().parent / "arcs_vet_deblend.py"
)
vd = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vd)


def clump_anti_fraction(g: float, anti_deg: float = 60.0, n: int = 200001) -> float:
    """Fraction of isotropic source separations whose image direction is ≥ ``anti_deg`` from the
    shear axis, for reduced shear ``g`` along x (image = A⁻¹ source, A ∝ diag(1 − g, 1 + g))."""
    phi = np.linspace(-np.pi / 2, np.pi / 2, n)
    x, y = np.cos(phi) / (1 - g), np.sin(phi) / (1 + g)
    return float(np.mean(np.degrees(np.arctan2(np.abs(y), np.abs(x))) >= anti_deg))


def third(cut: dict, m2: list) -> dict:
    sci, sub, err, wcs, px = vd.load(cut["f200w"])
    yy, xx = np.indices(sub.shape, dtype=float)
    cx, cy = 0.5 * (m2[0] + m2[4]), 0.5 * (m2[1] + m2[5])
    r = np.hypot(xx - cx, yy - cy) * px
    sm = gaussian_filter(np.nan_to_num(sub), 1.0)
    j, i = np.unravel_index(np.argmax(np.where((r > 0.25) & (r < 0.6), sm, -np.inf)), sm.shape)
    m = (np.hypot(xx - i, yy - j) <= 4) & np.isfinite(sub) & np.isfinite(err) & (err > 0)
    p3 = least_squares(
        lambda p: (vd._gauss(xx[m], yy[m], *p) - sub[m]) / err[m],
        [i, j, float(np.clip(sm[j, i] * 10, 1e-6, 999.0)), 1.0],
        bounds=([i - 3, j - 3, 0, 0.5], [i + 3, j + 3, 1e3, 4]),
    ).x
    knots = [(m2[0], m2[1], m2[3]), (m2[4], m2[5], m2[7]), (p3[0], p3[1], p3[3])]
    sky = [wcs.pixel_to_world(k[0], k[1]) for k in knots]
    sig = [k[2] * px for k in knots]
    sed = {}
    for b in vd.SW:
        _, s_b, e_b, w_b, px_b = vd.load(cut[b])
        kn = []
        for c, sa in zip(sky, sig, strict=True):
            kx, ky = w_b.world_to_pixel(c)
            kn.append((float(kx), float(ky), sa / px_b))
        mx, my = np.mean([k[0] for k in kn]), np.mean([k[1] for k in kn])
        f, e = vd.forced_fluxes(s_b, e_b, kn, 0.6 / px_b, mx, my)
        k = vd.to_ujy(px_b)
        sed[b] = [
            [round(float(a * k), 4), round(float(c * k), 4)] for a, c in zip(f, e, strict=True)
        ]

    def cmp(a, b):
        fa, ea = zip(*[sed[x][a] for x in vd.SW], strict=True)
        fb, eb = zip(*[sed[x][b] for x in vd.SW], strict=True)
        p, n = vd.ratio_consistency(fa, ea, fb, eb)
        return {"p": p, "n_bands": n, "class": vd.sed_class(p)}

    return {
        "third_sigma_arcsec": round(sig[2], 3),
        "third_sep_from_knot1_knot2_arcsec": [
            round(float(sky[2].separation(s).arcsec), 3) for s in sky[:2]
        ],
        "pa_knot1_to_third": round(float(sky[0].position_angle(sky[2]).deg % 180), 1),
        "sed_ujy_knot1_knot2_third": sed,
        "third_vs_knot1": cmp(2, 0),
        "third_vs_knot2": cmp(2, 1),
        "knot1_vs_knot2_with_third": cmp(0, 1),
    }


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--cut-dir", type=Path, required=True)
    p.add_argument(
        "--vet-dir",
        type=Path,
        default=paths.repo_root() / "results" / "arcs" / "abells1063-cats" / "vetting",
    )
    a = p.parse_args(argv)
    cut = {b: str(next((a.cut_dir / b).rglob("*.fits"))) for b in vd.SW}
    m2 = json.loads((a.vet_dir / "deblend_1557.json").read_text())["1557"]["observed"]["m2"]
    out = {
        "label": "post hoc (not in PLAN_1557.md); fits derived, fractions model_prediction",
        "1557_third_component": third(cut, m2),
        "clump_anti_fraction": {
            "assumptions": (
                "isotropic point-pair separations; reduced shear along the predicted axis; κ < 1"
            ),
            "g": {str(g): round(clump_anti_fraction(g), 4) for g in (0.43, 0.45, 1.5, 1.6, 1.9)},
        },
    }
    (a.vet_dir / "clumps_1557.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(out, indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
