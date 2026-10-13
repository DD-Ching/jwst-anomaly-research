"""Forward-modelled arc curvature: the curvature side each published model predicts for an arc.

``arcs_curvature.py`` compares an arc's measured curvature with the mass direction -alpha, a
circular-lens proxy that fails near saddle points (D-024 addendum 5). Here the model itself
predicts the shape: a round Gaussian source (``--sigma-src`` arcsec, ASSUMPTION) placed at the
source-plane position of the arc's catalogue centre (its photo-z, else ``--z-default``) is
ray-traced onto a ``--size`` arcsec image-plane grid. The predicted image is cut at
``--frac`` of its peak (connected component at the centre) and fitted with the same parabola
(:func:`arcs_curvature.fit_curvature`). The predicted and
measured curvature sides, taken through their signs against -alpha (``cos_mass_obs``,
``cos_mass_model``): equal signs mean the arc bends the way the model says.

Only sources whose stored ``curvature.ecsv`` row is curved (``agree`` or ``reverse``) are compared.
``model_sig`` is |c_model| / sigma_c(observed): the predicted bend in units of the measurement
error. A predicted bend below ``--min-sig`` sigma is not testable (ASSUMPTION). Outputs
``results/arcs_curvature/<model>/forward.ecsv`` and ``forward.json`` (model_prediction against
derived).

    python scripts/arcs_curvature_model.py --model elgordo-caminha23
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy import ndimage

from jwst_anomaly import paths, schema

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ac = _load("arcs_curvature")
lc = ac.lc

ASSUMPTIONS = ("sigma_src", "size", "step", "frac", "z_default", "min_sig")


def predicted_fit(model, ra: float, dec: float, z_s: float, args) -> dict:
    """Parabola fit to the model image of a round source at the source position of (ra, dec)."""
    x0, y0 = (float(v) for v in model.to_frame(ra, dec))
    h = args.size / 2.0
    g = np.arange(-h, h + args.step / 2, args.step)
    gx, gy = np.meshgrid(x0 + g, y0 + g)
    bx, by = model.source_points(gx, gy, z_s)
    b0x, b0y = model.source_points(np.array([x0]), np.array([y0]), z_s)
    img = np.exp(-((bx - b0x[0]) ** 2 + (by - b0y[0]) ** 2) / (2 * args.sigma_src**2))
    lab, _ = ndimage.label(img >= args.frac * img.max())
    c = len(g) // 2
    if lab[c, c] == 0:
        return {}
    m = lab == lab[c, c]
    east, north = -(gx[m] - x0), gy[m] - y0  # model frame x points West
    return ac.fit_curvature(east, north, img[m])


def _same(t) -> np.ndarray:
    """True where the predicted and measured curvature sides fall on the same side of -alpha."""
    return np.sign(np.asarray(t["cos_mass_obs"], float)) == np.sign(
        np.asarray(t["cos_mass_model"], float)
    )


def run(args) -> dict:
    obs = Table.read(args.results / args.model / "curvature.ecsv")
    arcs = Table.read(args.arcs_dir / args.model / "arcs.ecsv")
    zmap = {int(r["label"]): float(r["z_phot"]) for r in arcs}
    obs = obs[np.isin(obs["curvature_class"], ["agree", "reverse"])]
    model, _, _ = lc.load_model(args.model)
    lc.apply_frame_offset(args.model, model)
    rows = []
    for r in obs:
        z = zmap.get(int(r["label"]), np.nan)
        z = z if np.isfinite(z) and z > model.z_lens + 0.1 else args.z_default
        fit = predicted_fit(model, float(r["ra"]), float(r["dec"]), z, args)
        # both curvature sides lie along nearly parallel minor axes, so compare their signs
        # against one reference, the mass direction -alpha (east = +alpha_x, north = -alpha_y)
        d = model.deflection([r["ra"]], [r["dec"]], [z])
        me, mn = float(d["alpha_x"][0]), -float(d["alpha_y"][0])
        norm = np.hypot(me, mn)
        me, mn = me / norm, mn / norm
        c_err = abs(float(r["c"])) / float(r["c_sig"])
        rows.append(
            {
                "label": int(r["label"]),
                "z_used": z,
                "curvature_class": str(r["curvature_class"]),
                "testable": bool(r["testable"]),
                "c_obs": float(r["c"]),
                "c_model": fit.get("c", np.nan),
                "length_model": fit.get("length", np.nan),
                "model_sig": abs(fit["c"]) / c_err if "c" in fit else np.nan,
                "cos_mass_obs": float(r["cos_mass"]),
                "cos_mass_model": fit["side_e"] * me + fit["side_n"] * mn if "c" in fit else np.nan,
            }
        )
    t = Table(rows=rows) if rows else Table()
    t.meta.update(
        provenance=schema.Provenance.MODEL_PREDICTION.value,
        source=f"{args.model} ray-traced round source vs results/arcs_curvature/{args.model}",
        assumptions={k: getattr(args, k) for k in ASSUMPTIONS},
    )
    test = t[np.asarray(t["model_sig"], float) >= args.min_sig] if rows else []
    summary = {
        "model": args.model,
        "assumptions": {k: getattr(args, k) for k in ASSUMPTIONS},
        "n_curved": len(t),
        "n_model_testable": len(test),
        "consistent": int(np.sum(_same(test))) if len(test) else 0,
        "inconsistent": [
            {
                "label": int(r["label"]),
                "model_sig": round(float(r["model_sig"]), 1),
                "cos_mass_obs": round(float(r["cos_mass_obs"]), 2),
                "cos_mass_model": round(float(r["cos_mass_model"]), 2),
            }
            for r, ok in zip(test, _same(test), strict=True)
            if not ok
        ],
    }
    out = args.results / args.model
    t.write(out / "forward.ecsv", overwrite=True)
    (out / "forward.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--model", required=True, choices=sorted(ac.IMAGES))
    p.add_argument("--arcs-dir", type=Path, default=paths.repo_root() / "results" / "arcs")
    p.add_argument("--results", type=Path, default=paths.repo_root() / "results" / "arcs_curvature")
    p.add_argument("--sigma-src", type=float, default=0.1, help="source Gaussian sigma (arcsec)")
    p.add_argument("--size", type=float, default=8.0, help="image-plane grid side (arcsec)")
    p.add_argument("--step", type=float, default=0.03, help="grid step (arcsec)")
    p.add_argument("--frac", type=float, default=0.1, help="cut as a fraction of the peak")
    p.add_argument("--z-default", type=float, default=2.0)
    p.add_argument("--min-sig", type=float, default=3.0)
    args = p.parse_args(argv)
    print(json.dumps(run(args), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
