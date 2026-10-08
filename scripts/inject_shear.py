"""Injection-recovery of W1 negative-mass lenses through the ``shear`` screen (D-050).

The same simulated lenses as ``inject_radial.py`` (D-049: ``paint_lens``, lens of mass |M| at the
cluster redshift, the field's own background rows lensed, blends, detection floor), measured with
the catalogue aperture-mass screen of ``exotic_screens.py shear`` instead of ``radial``. The
efficiencies decide whether the shear screen replaces ``radial`` for W1 (D-050 "Adoption").

Recovery (ASSUMPTIONs): the largest S of the screen's grid centres within ``--recover-tol`` of the
injected centre has ``p_random`` < 0.05, i.e. fewer than 5 % of the rotation-null draws reach that
value anywhere in the field, *and* exceeds the real field's largest B-mode excursion |S_×| (D-053:
the real E and B maps have the same heavier-than-rotation tails, so shape systematics, not the
rotation null, set the floor; the real field is null under the same rule). Each batch of 10
trials draws its own ``--n-random`` rotation null of the *real* field (the injected rows change
the field maximum only near the lens; ASSUMPTION). Rows that ``paint_lens`` keeps use the real
field's corrected ellipticities. Painted images get theirs from their lensed moments, minus R
times the cluster model's reduced shear at their position, and keep only R of the lens-induced
change, because the catalogue's isophotal moments respond to shear by R ~ 0.45 (D-053).

Everything injected is ``simulated``; efficiencies and limits are ``derived``.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from astropy.table import Table

from jwst_anomaly import paths, schema

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ir = _load("inject_radial")
es = ir.es

SHEAR_DEFAULTS = es.shear_defaults()
# fields with photo-z whose inputs are public without tarballs (DJA v7.4/v7.0 tarballs 404 on
# 2026-10-08, so SMACS and El Gordo have no photo-z here and stay out of the headline)
DEFAULT_FIELDS = ("abell2744", "macs0416", "macs1149", "abell370")


def injected_ellipticity(e_src, raw_img, raw_src, r: float) -> np.ndarray:
    """Corrected ε of a painted image: its source's corrected ε plus R times the lens-induced
    change of the *measured* moments (raw image minus raw source). The catalogue's moments respond
    to shear by R (D-053), the painted moments by 1; taking the change between raw moments keeps
    the cluster shear out of it. |ε| is capped at 0.99. ``ShearInjector`` paints only resolved
    sources; for an unresolved one (NaN) this returns R times the raw image ε."""
    e_src, raw_src = np.asarray(e_src, complex), np.asarray(raw_src, complex)
    resolved = np.isfinite(raw_src) & np.isfinite(e_src)
    with np.errstate(invalid="ignore"):
        lensed = e_src + r * (np.asarray(raw_img, complex) - raw_src)
    out = np.where(resolved, lensed, r * np.asarray(raw_img, complex))
    mod = np.abs(out)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(mod > 0.99, out * 0.99 / mod, out)


class ShearInjector:
    """The shear screen on one field: real ellipticities cached, recovery checked near the lens."""

    def __init__(self, model, shapes: Table, args, psf_sigma: float, extra: Table | None = None):
        shapes = shapes.copy(copy_data=False)  # the caller's table gets no `_src` column
        self.model, self.shapes, self.args, self.psf = model, shapes, args, psf_sigma
        self.base = es.shear_screen(model, shapes, args, psf_sigma, extra)
        if not np.isfinite(self.base["s"]).any():
            raise ValueError("no grid centre has enough shear sources")
        self.e = self.base["e"]
        self.e_all = self.base["e_all"]
        self.eps = self.base["eps"]
        self.r = self.base["counts"]["responsivity"]
        xs, ys = model.to_frame(shapes["ra"], shapes["dec"])
        self.xs, self.ys = np.asarray(xs), np.asarray(ys)
        self.use = self.base["use"]
        self.ap = self.base["ap"]
        self.null = self.base["null_max"]
        self.b_max = float(np.nanmax(np.abs(self.base["s_cross"])))
        # spike segments are never painted: the real field vetoes them too
        # and sources where the cluster shear cannot be removed (κ >= 1, |g| or |R g| >= 1) are
        # not painted
        # only resolved sources are lensed: an unresolved one has no measured shape to paint from
        self.background = (
            es.lensable_mask(shapes, model.z_lens)
            & ~self.base["spike"]
            & self.base["correctable"]
            & np.isfinite(self.eps)
        )
        shapes["_src"] = np.arange(len(shapes))  # carried into painted rows by paint_lens
        self.next_label = int(np.max(shapes["label"])) + 1

    def draw_null(self, rng) -> None:
        self.null = self.ap.null_max(self.e[self.use], self.args.n_random, rng)

    def p_random(self, s: float) -> float:
        return float(np.mean(self.null >= s)) if np.isfinite(s) else 1.0

    def trial(self, x_l, y_l, theta_row, pixel_scale) -> dict:
        keep, img, info = ir.paint_lens(
            self.shapes,
            self.xs,
            self.ys,
            self.background,
            theta_row,
            self.model,
            x_l,
            y_l,
            self.psf,
            self.next_label,
            pixel_scale,
        )
        a = self.args
        near_l = np.hypot(self.xs - x_l, self.ys - y_l) <= a.recover_tol + a.aperture_arcsec + 1.0
        k = keep & self.use & near_l
        x, y, e = self.xs[k], self.ys[k], self.e[k]
        n_img_used = 0
        if len(img):
            e_img, c = es.shear_sources(
                self.model, img, self.psf, a.min_snr, a.max_g, self.r, spike_veto=False
            )
            src = np.asarray(img["_src"], int)
            e_img = np.where(
                np.isfinite(e_img),
                injected_ellipticity(self.e_all[src], c["eps"], self.eps[src], self.r),
                np.nan + 0j,
            )
            ok = np.isfinite(e_img)
            n_img_used = int(ok.sum())
            xi, yi = self.model.to_frame(img["ra"], img["dec"])
            x = np.concatenate([x, np.atleast_1d(xi)[ok]])
            y = np.concatenate([y, np.atleast_1d(yi)[ok]])
            e = np.concatenate([e, e_img[ok]])
        gx, gy = self.base["gx"], self.base["gy"]
        near = np.hypot(gx - x_l, gy - y_l) <= a.recover_tol
        cx, cy = gx[near], gy[near]
        ap = es.ApertureMass(
            x, y, cx, cy, a.aperture_arcsec, a.filter, a.r_min_arcsec, a.min_sources
        )
        s, _ = ap.snr(e)
        best = float(np.nanmax(s)) if np.isfinite(s).any() else float("nan")
        p = self.p_random(best)
        return info | {
            "x": float(x_l),
            "y": float(y_l),
            "n_image_used": n_img_used,
            "s_best": best,
            "p_random": p,
            "recovered": bool(p < ir.P_RECOVER and best > self.b_max),
        }


def run_field(name: str, args) -> dict:
    t0 = time.time()
    out = args.out / name
    out.mkdir(parents=True, exist_ok=True)
    model, shapes, extra, cat, photoz = ir.load_field(name, out)
    sargs = SimpleNamespace(**SHEAR_DEFAULTS)
    for k in ("aperture_arcsec", "filter", "r_min_arcsec"):
        setattr(sargs, k, getattr(args, k))
    sargs.max_radius = ir.FIELDS[name]["max_radius"]
    sargs.recover_tol = args.recover_tol
    psf = es.psf_sigma_px(shapes)
    inj = ShearInjector(model, shapes, sargs, psf, extra)
    t_base = time.time() - t0
    scale = ir.pixel_scale_arcsec(cat)
    fx, fy, area = ir.screened_footprint(inj.xs, inj.ys, sargs.max_radius)
    border = ir.footprint_border(inj.xs, inj.ys, sargs.max_radius)
    z_src = ir.source_redshift(shapes)
    s_map = inj.base["s"]
    s_max = float(np.nanmax(s_map))
    iy, ix = np.unravel_index(np.nanargmax(s_map), s_map.shape)
    base_null = inj.base["null_max"]
    seeds = np.random.SeedSequence(args.seed).spawn(len(args.mass))
    trials, eff = [], {}
    for mass, seed in zip(args.mass, seeds, strict=True):
        rng = np.random.default_rng(seed)
        theta_row = np.where(inj.background, ir.theta_e_arcsec(mass, model.z_lens, z_src), np.nan)
        theta_2 = float(ir.theta_e_arcsec(mass, model.z_lens, ir.Z_SOURCE)[0])
        sub = []
        for k in range(args.n_inject):
            if k % ir.BATCH == 0:
                inj.draw_null(rng)
            i = rng.integers(0, len(fx))
            jx, jy = rng.uniform(-ir.FOOTPRINT_STEP / 2, ir.FOOTPRINT_STEP / 2, 2)
            r = inj.trial(fx[i] + jx, fy[i] + jy, theta_row, scale)
            r["mass_msun"] = mass
            sub.append(r)
        trials.extend(sub)
        n_rec = int(sum(r["recovered"] for r in sub))
        p = n_rec / args.n_inject
        eff[f"{mass:.0e}"] = {
            "mass_msun": mass,
            "theta_e_zs2_arcsec": theta_2,
            "n_inject": args.n_inject,
            "n_recovered": n_rec,
            "efficiency": p,
            "efficiency_err": float(np.sqrt(p * (1 - p) / args.n_inject)),
            "mean_lensed": float(np.mean([r["n_lensed"] for r in sub])),
            "mean_image_used": float(np.mean([r["n_image_used"] for r in sub])),
            "median_s_best": float(np.nanmedian([r["s_best"] for r in sub])),
        }
        print(f"{name} {mass:.0e}: {n_rec}/{args.n_inject}", flush=True)
    table = Table(rows=trials)
    table.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"inject_shear.py {name}: simulated W1 lenses in {cat}",
    )
    table.write(out / "trials.ecsv", overwrite=True)
    finite = np.isfinite(s_map)
    ra_pk, dec_pk = model.to_sky(inj.base["gx"][iy, ix], inj.base["gy"][iy, ix])
    summary = {
        "field": name,
        "model": ir.FIELDS[name]["model"],
        "z_lens": model.z_lens,
        "catalog": str(cat),
        "photoz": str(photoz) if photoz else None,
        "screen": inj.base["counts"]
        | {
            "psf_sigma_px": psf,
            "n_centres_valid": int(finite.sum()),
            "s_max": s_max,
            "s_max_xy": [float(inj.base["gx"][iy, ix]), float(inj.base["gy"][iy, ix])],
            "s_max_radec": [float(ra_pk), float(dec_pk)],
            "p_random_max": float(np.mean(base_null >= s_max)),
            "b_max": inj.b_max,
            "e_exceeds_b": bool(s_max > inj.b_max),
            "passes_eb_rule": bool(
                np.mean(base_null >= s_max) < ir.P_RECOVER and s_max > inj.b_max
            ),
            **es.cross_p_values(inj.base["s_cross"], base_null),
            "null_max_p50_p95": [float(v) for v in np.percentile(base_null, [50, 95])],
        },
        "screened_area_arcsec2": area,
        "screened_area_deg2": area / 3600.0**2,
        "footprint_border": border,
        "efficiency": eff,
        "wall_time_s": {"base_screen": round(t_base, 1), "total": round(time.time() - t0, 1)},
        "assumptions": {
            "recover_tol_arcsec": args.recover_tol,
            "p_recover": ir.P_RECOVER,
            "trials_per_null": ir.BATCH,
            "seed": args.seed,
            "n_inject": args.n_inject,
            "spike_stars": len(extra) if extra is not None else None,
            "screen": SHEAR_DEFAULTS
            | {k: getattr(sargs, k) for k in ("aperture_arcsec", "filter", "r_min_arcsec")}
            | {"max_radius": sargs.max_radius},
        },
    }
    (out / "summary.json").write_text(json.dumps(es.json_safe(summary), indent=1))
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fields", nargs="+", default=list(DEFAULT_FIELDS), choices=sorted(ir.FIELDS))
    ap.add_argument("--mass", nargs="+", type=float, default=list(ir.MASSES), help="|M| in M_sun")
    ap.add_argument("--n-inject", type=int, default=200, help="lenses per field and mass")
    ap.add_argument("--recover-tol", type=float, default=ir.RECOVER_TOL)
    ap.add_argument("--aperture-arcsec", type=float, default=SHEAR_DEFAULTS["aperture_arcsec"])
    ap.add_argument(
        "--filter", default=SHEAR_DEFAULTS["filter"], choices=("pointmass", "schirmer", "tophat")
    )
    ap.add_argument("--r-min-arcsec", type=float, default=SHEAR_DEFAULTS["r_min_arcsec"])
    ap.add_argument("--seed", type=int, default=50)
    ap.add_argument("--out", type=Path, default=paths.outputs_dir() / "inject_shear")
    ap.add_argument(
        "--combine-only",
        action="store_true",
        help="combine existing <out>/<field>/summary.json files (fields run as separate processes)",
    )
    args = ap.parse_args(argv)
    keys = [f"{m:.0e}" for m in args.mass]
    if len(set(keys)) != len(keys):
        ap.error(f"--mass values must differ at one significant figure (keys {keys})")
    if args.combine_only:
        summaries = [json.loads((args.out / f / "summary.json").read_text()) for f in args.fields]
    else:
        summaries = [run_field(f, args) for f in args.fields]
    combined = ir.combine(summaries, args.mass)

    def setting(s):  # grid half-width and star table are per field; everything else must agree
        a = json.loads(json.dumps(s["assumptions"]))
        a["screen"].pop("max_radius", None)
        a.pop("spike_stars", None)
        return a

    if len({json.dumps(setting(s), sort_keys=True) for s in summaries}) != 1:
        ap.error("the fields were run with different screen settings; re-run them alike")
    combined["screen"] = setting(summaries[0])
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "combined.json").write_text(json.dumps(es.json_safe(combined), indent=1))
    print(json.dumps(combined, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
