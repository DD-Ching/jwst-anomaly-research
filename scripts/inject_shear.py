"""Injection-recovery of W1 negative-mass lenses through the ``shear`` screen (D-050).

The simulated lenses of ``inject_radial.py`` (D-049: ``paint_lens``, lens of mass |M| at the
cluster redshift, the field's own background rows lensed, blends, detection floor), measured with
the catalogue aperture-mass screen of ``exotic_screens.py shear`` instead of ``radial``. One
difference (ASSUMPTION, conservative): only rows whose shape can be corrected are lensed, i.e.
resolved after PSF deconvolution, not spike segments, and with κ < 1, |g| < 1, |R g| < 1; the S/N,
``max_g`` and spike-segment cuts apply to the painted images, since magnification and position
change all three. The
radial injections also lens unresolved rows, so the efficiency ratio mixes screen power with
this (docs/exotic_limits.md caveats). The efficiencies decide whether the shear screen replaces
``radial`` for W1 (D-050).

Recovery (ASSUMPTIONs): the largest S of the screen's grid centres within ``--recover-tol`` of the
injected centre has ``p_random`` < 0.05, i.e. fewer than 5 % of the rotation-null draws reach that
value anywhere in the field, *and* exceeds the real field's largest B-mode excursion |S_×| (D-053:
the real E and B maps have the same heavier-than-rotation tails, so shape systematics, not the
rotation null, set the floor; the real field is null under the same rule). Each batch of 10 trials
draws its own rotation null (``SHEAR_DEFAULTS["n_random"]`` draws, not a CLI option) of the *real*
field (the injected rows change the field maximum only near the lens; ASSUMPTION). Rows that
``paint_lens`` keeps use the real field's corrected ellipticities. A painted image gets its source's
corrected ε plus R times the change of the raw (measured) moments between image and source, because
the catalogue's isophotal moments respond to shear by R ~ 0.45 (D-053); the sum is linear and |ε|
capped at 0.99 (ASSUMPTION; how often the cap binds is not measured).

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
from astropy.table import Table, vstack

from jwst_anomaly import paths, schema

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ir = _load("inject_radial")
es = ir.es
wsm = _load("w1_streak_mask")

SHEAR_DEFAULTS = es.shear_defaults()
# fields with public photo-z (SMACS and El Gordo added 2026-10-17: their DJA tarballs answer again)
DEFAULT_FIELDS = (
    "abell2744",
    "macs0416",
    "macs1149",
    "abell370",
    "abells1063",
    "smacs0723",
    "elgordo",
    "macs0417",
    "whl0137",
    "rxj2129",
    "macs0647",
)  # macs1423 is screened but not combined by default: results/w1_shear_more/RESULT.md
# Screen options added after summaries were first committed, with the value those runs used
# (before the option existed). A summary without the key is read with this value.
LEGACY_SCREEN_OPTIONS = {"psf_chi": [0.0, 0.0], "min_coverage": 0.0}
# Shear-only inputs that differ from the radial screen's FIELDS entry. The shear screen needs
# photo-z (cluster shear at each row's z), so Abell S1063 uses the DJA v7.5 zout here
# (docs/fields/abells1063.md), with the Gaia spike veto (D-043).
SHEAR_FIELD_OVERRIDES = {
    "abells1063": {
        "photoz": (
            "eazy",
            "cache/dja/photoz-abells1063-v7.5/abells1063-grizli-v7.5-fix.eazypy.zout.fits",
        ),
        "spike_stars": True,
        "centre": (342.1832, -44.5309),
    },
    # photo-z already in FIELDS (DJA tarballs reachable 2026-10-17); Gaia veto at the model centre
    "smacs0723": {"spike_stars": True, "centre": (110.826989, -73.454723)},
    "elgordo": {"spike_stars": True, "centre": (15.740693, -49.271992)},
    # Shear-only fields (2026-10-17; results/w1_shear_more/PLAN.md): JWST-era map models that
    # passed the `arcs` alignment control, the photo-z of those runs, the 100" grid of the other
    # CANUCS fields (ASSUMPTION) and the Gaia veto at the model reference point.
    "macs0417": {
        "model": "macs0417-canucs-maps",
        "catalog": ir._MAST.format("jw01208-o006_t003_nircam_clear-f200w"),
        "photoz": ("canucs", ir._CANUCS.format("macs0417")),
        "max_radius": 100.0,
        "spike_stars": True,
        "centre": (64.394582, -11.908964),
    },
    "macs1423": {
        "model": "macs1423-canucs-maps",
        "catalog": ir._MAST.format("jw01208-o010_t005_nircam_clear-f200w"),
        "photoz": ("canucs", ir._CANUCS.format("macs1423")),
        "max_radius": 100.0,
        "spike_stars": True,
        "centre": (215.949485, 24.078455),
    },
    "whl0137": {
        "model": "whl0137-relics-lenstool",
        "catalog": ir._MAST.format("jw02282-o010_t001_nircam_clear-f200w"),
        "photoz": (
            "eazy",
            "cache/dja/photoz-sunrise-v7.5/sunrise-grizli-v7.5-fix.eazypy.zout.fits",
        ),
        "max_radius": 100.0,
        "spike_stars": True,
        "centre": (24.354, -8.4563),
    },
    # 2026-10-17 (results/w1_shear_rxj_m0647/PLAN.md): the primary models and DJA v7.0 photo-z of
    # the `arcs` runs; MACS0647 uses epoch 1 (program 1433 o010)
    "rxj2129": {
        "model": "rxj2129-lenstool-maps",
        "catalog": ir._MAST.format("jw02767-o002_t001_nircam_clear-f200w"),
        "photoz": (
            "eazy",
            "cache/dja/photoz-rxj2129-v7.0/rxj2129-grizli-v7.0-fix.eazypy.zout.fits",
        ),
        "max_radius": 100.0,
        "spike_stars": True,
        "centre": (322.416503, 0.089236),
    },
    "macs0647": {
        "model": "macs0647-clash-nfw",
        "catalog": ir._MAST.format("jw01433-o010_t001_nircam_clear-f200w"),
        "photoz": (
            "eazy",
            "cache/dja/photoz-macs0647-v7.0/macs0647-grizli-v7.0-fix.eazypy.zout.fits",
        ),
        "max_radius": 100.0,
        "spike_stars": True,
        "centre": (101.959458, 70.248611),
    },
    # 2026-10-17 (results/w1_shear_macs0717/PLAN.md): DJA has no MACS0717 photo-z; HST-based
    # ASTRODEEP-FF redshifts (scripts/astrodeep_zout.py) over the HFF WFC3 footprint
    "macs0717": {
        "photoz": ("eazy", "cache/external/astrodeep/m0717clu_zout.fits"),
        "max_radius": 100.0,
        "spike_stars": True,
        "centre": (109.3864, 37.7550),  # VENUS o029 pointing (configs/macs0717.yaml)
    },
}
# 2026-10-18 (results/w1_shear_ncf/PLAN.md): the CANUCS NIRCam flanking fields, no lens model
# (``inject_radial.NullLensModel`` at the catalogue's median position), F150W (no F200W there),
# CANUCS NCF photo-z, R fixed by gate R (``w1_ncf_gate_r.py``; the NCF have no shear to fit it)
NCF_FIELDS = {
    "macs0416-ncf": ("o023_t002", "macs0416", (64.14500, -24.09645)),
    "macs1149-ncf": ("o029_t004", "macs1149", (177.40342, 22.29696)),
    "abell370-ncf": ("o020_t001", "a370", (40.06161, -1.62260)),
    "macs0417-ncf": ("o026_t003", "macs0417", (64.36952, -11.82049)),
    "macs1423-ncf": ("o032_t005", "macs1423", (215.87389, 24.15656)),
}
GATE_R = _DIR.parent / "results" / "w1_shear_ncf" / "gate_r.json"
_NCF_CANUCS = (
    "cache/external/canucs/hlsp_canucs_jwst-hst_multi_{0}-ncf_multi_v1_photometry-cat.fits.gz"
)
for _name, (_obs, _c, _centre) in NCF_FIELDS.items():
    SHEAR_FIELD_OVERRIDES[_name] = {
        "model": None,
        "catalog": None,  # set below by ncf_catalog
        "photoz": ("canucs", _NCF_CANUCS.format(_c)),
        "max_radius": 200.0,
        "spike_stars": True,
        "centre": _centre,
        "responsivity": "ncf",
    }


def ncf_catalog(obs: str) -> str:
    """Data-root path of a program 1208 F150W pipeline catalogue (NCF screens and gate R)."""
    return ir._MAST.format(f"jw01208-{obs}_nircam_clear-f150w")


for _name, (_obs, _c, _centre) in NCF_FIELDS.items():
    SHEAR_FIELD_OVERRIDES[_name]["catalog"] = ncf_catalog(_obs)
# 2026-10-18 (results/w1_shear_goodsn/PLAN.md): six JADES GOODS-N F200W pointings (program
# 1181) with small mutual overlap, no lens model, DJA v7.3 photo-z (tarball member, D-051),
# R fixed at GOODSN_R (the smallest F200W cluster-fit R apart from MACS0647; ASSUMPTION)
GOODSN_FIELDS = {
    "goodsn-o001": ("o001_t001", (189.18115, 62.27479)),
    "goodsn-o003": ("o003_t007", (189.04836, 62.26091)),
    "goodsn-o004": ("o004_t008", (189.28111, 62.20377)),
    "goodsn-o005": ("o005_t008", (189.21071, 62.18853)),
    "goodsn-o009": ("o009_t009", (189.42173, 62.25673)),
    "goodsn-o010": ("o010_t009", (189.37152, 62.30012)),
}
GOODSN_R = 0.414  # Abell 370 F200W fit (results/w1_volume/summaries/z0.3_abell370.json), rounded
_GOODSN_ZOUT = "cache/dja/photoz-goodsn-v7.3/gdn-grizli-v7.3-fix.eazypy.zout.fits"
for _name, (_obs, _centre) in GOODSN_FIELDS.items():
    SHEAR_FIELD_OVERRIDES[_name] = {
        "model": None,
        "catalog": ir._MAST.format(f"jw01181-{_obs}_nircam_clear-f200w"),
        "photoz": ("eazy", _GOODSN_ZOUT),
        "max_radius": 200.0,
        "spike_stars": True,
        "centre": _centre,
        "responsivity": GOODSN_R,
    }
FIELD_NAMES = sorted(set(ir.FIELDS) | set(SHEAR_FIELD_OVERRIDES))


def fixed_responsivity(name: str, spec: dict, lens_z, gate: Path = GATE_R) -> float | None:
    """R for a field without a lens model: gate R's ``adopted_r`` when
    ``spec["responsivity"] == "ncf"`` (an error if gate R has not passed), the number itself when
    it is a number (GOODS-N); an error if no ``--lens-z`` is given (``z_lens`` 0). None for every
    other field (R is fitted)."""
    rspec = spec.get("responsivity")
    if rspec is None:
        return None
    if lens_z is None:
        raise ValueError(f"{name} has no lens model; give --lens-z")
    if rspec != "ncf":
        return float(rspec)
    res = json.loads(Path(gate).read_text())["gate_r"] if Path(gate).exists() else {}
    if not res.get("passed"):
        raise ValueError(f"{name}: gate R has not passed ({gate}; w1_ncf_gate_r.py)")
    return float(res["adopted_r"])


def field_spec(name: str) -> dict:
    return {**ir.FIELDS.get(name, {}), **SHEAR_FIELD_OVERRIDES.get(name, {})}


def injected_ellipticity(e_src, raw_img, raw_src, r: float) -> np.ndarray:
    """Corrected ε of a painted image: its source's corrected ε plus R times the lens-induced
    change of the *measured* moments (raw image minus raw source). The catalogue's moments respond
    to shear by R (D-053), the painted moments by 1; taking the change between raw moments keeps
    the cluster shear out of it. |ε| is capped at 0.99. An unresolved source (NaN) gives NaN:
    ``ShearInjector`` paints only resolved sources."""
    e_src, raw_src = np.asarray(e_src, complex), np.asarray(raw_src, complex)
    out = e_src + r * (np.asarray(raw_img, complex) - raw_src)
    mod = np.abs(out)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(mod > 0.99, out * 0.99 / mod, out)


class ShearInjector:
    """The shear screen on one field: real ellipticities cached, recovery checked near the lens."""

    def __init__(
        self,
        model,
        shapes: Table,
        args,
        psf_sigma: float,
        extra: Table | None = None,
        streak=None,
    ):
        shapes = shapes.copy(copy_data=False)  # the caller's table gets no `_src` column
        self.model, self.shapes, self.args, self.psf = model, shapes, args, psf_sigma
        # streak(ra, dec) -> bool: rows the field's streak mask vetoes (results/w1_streak_veto);
        # shear_screen blanks real rows on streaks, painted images on streaks are vetoed here
        self.streak = streak
        on = None
        if streak is not None:
            on = streak(np.asarray(shapes["ra"], float), np.asarray(shapes["dec"], float))
        self.n_streak = None if on is None else int(on.sum())
        self.base = es.shear_screen(model, shapes, args, psf_sigma, extra, on)
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
        # only rows whose shape can be corrected are lensed (ASSUMPTION, conservative; docstring)
        self.background = (
            es.lensable_mask(shapes, model.z_lens)
            & ~self.base["spike"]
            & self.base["correctable"]
            & np.isfinite(self.eps)
            & (True if on is None else ~on)  # streak rows are not galaxies to lens
        )
        shapes["_src"] = np.arange(len(shapes))  # carried into painted rows by paint_lens
        self.next_label = int(np.max(shapes["label"])) + 1
        # painted images face the real field's spike veto: the field's spike segments go in with
        # them so the spike axes are estimated as on the real field (ASSUMPTION: the vetoed real
        # segments dominate that estimate)
        self.extra = extra
        self.spike_ref = self._spike_cols(shapes[self.base["spike"]])

    @staticmethod
    def _spike_cols(t: Table) -> Table:
        return Table({k: np.asarray(t[k], float) for k in ("ra", "dec", "pa_obs")})

    def spike_veto(self, img: Table) -> np.ndarray:
        """Painted rows that ``spike_segments`` would veto on the real field."""
        both = vstack([self.spike_ref, self._spike_cols(img)])
        return es.spike_segments(both, self.shapes, self.extra)[len(self.spike_ref) :]

    def draw_null(self, rng) -> None:
        self.null = self.ap.null_max(self.e[self.use], self.args.n_random, rng)

    def p_random(self, s: float) -> float:
        return float(np.mean(self.null >= s)) if np.isfinite(s) else 1.0

    def lensed_sample(self, x_l, y_l, theta_row, pixel_scale):
        """Model-frame positions and corrected ε of the used rows near one lens painted at
        (``x_l``, ``y_l``): kept real rows within ``recover_tol + aperture + 1``″ plus the painted
        images that pass the screen's cuts. Returns ``x, y, e, info, n_image_used``."""
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
            psf_chi=getattr(self.args, "psf_chi", 0j),
        )
        a = self.args
        near_l = np.hypot(self.xs - x_l, self.ys - y_l) <= a.recover_tol + a.aperture_arcsec + 1.0
        k = keep & self.use & near_l
        x, y, e = self.xs[k], self.ys[k], self.e[k]
        n_img_used = 0
        if len(img):
            e_img, c = es.shear_sources(
                self.model,
                img,
                self.psf,
                a.min_snr,
                a.max_g,
                self.r,
                spike_veto=False,
                psf_chi=getattr(a, "psf_chi", 0j),
            )
            src = np.asarray(img["_src"], int)
            on_streak = (
                np.zeros(len(img), bool)
                if self.streak is None
                else self.streak(np.asarray(img["ra"], float), np.asarray(img["dec"], float))
            )
            e_img = np.where(
                np.isfinite(e_img) & ~self.spike_veto(img) & ~on_streak,
                injected_ellipticity(self.e_all[src], c["eps"], self.eps[src], self.r),
                np.nan + 0j,
            )
            ok = np.isfinite(e_img)
            n_img_used = int(ok.sum())
            xi, yi = self.model.to_frame(img["ra"], img["dec"])
            x = np.concatenate([x, np.atleast_1d(xi)[ok]])
            y = np.concatenate([y, np.atleast_1d(yi)[ok]])
            e = np.concatenate([e, e_img[ok]])
        return x, y, e, info, n_img_used

    def trial(self, x_l, y_l, theta_row, pixel_scale) -> dict:
        x, y, e, info, n_img_used = self.lensed_sample(x_l, y_l, theta_row, pixel_scale)
        a = self.args
        gx, gy = self.base["gx"], self.base["gy"]
        near = np.hypot(gx - x_l, gy - y_l) <= a.recover_tol
        cx, cy = gx[near], gy[near]
        valid = None if self.base["valid"] is None else self.base["valid"][near]
        ap = es.ApertureMass(
            x, y, cx, cy, a.aperture_arcsec, a.filter, a.r_min_arcsec, a.min_sources, valid
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
    spec = field_spec(name)
    model, shapes, extra, cat, photoz = ir.load_field(name, out, spec)
    streak = None
    if getattr(args, "streak_masks", None) is not None:
        streak = streak_veto(args.streak_masks, name)
    sargs = SimpleNamespace(**SHEAR_DEFAULTS)
    for k in ("aperture_arcsec", "filter", "r_min_arcsec", "min_coverage"):
        setattr(sargs, k, getattr(args, k))
    sargs.max_radius = spec["max_radius"]
    sargs.recover_tol = args.recover_tol
    sargs.psf_chi = complex(*getattr(args, "psf_chi", (0.0, 0.0)))
    r_fixed = fixed_responsivity(name, spec, getattr(args, "lens_z", None))
    if r_fixed is not None:  # no lens model: R fixed by gate R, the lens redshift explicit
        sargs.responsivity = r_fixed
    psf = es.psf_sigma_px(shapes)
    inj = ShearInjector(model, shapes, sargs, psf, extra, streak)
    t_base = time.time() - t0
    scale = ir.pixel_scale_arcsec(cat)
    fx, fy, area = ir.screened_footprint(inj.xs, inj.ys, sargs.max_radius)
    border = ir.footprint_border(inj.xs, inj.ys, sargs.max_radius)
    z_src = ir.source_redshift(shapes)
    # W1 lens redshift (results/w1_volume/PLAN.md): the cluster's unless --lens-z is given; then
    # only screen rows behind that redshift too are painted (single-plane ASSUMPTION: the cluster
    # model and its shear removal stay as they are)
    z_w1 = model.z_lens if getattr(args, "lens_z", None) is None else float(args.lens_z)
    if z_w1 != model.z_lens:
        inj.background = inj.background & es.lensable_mask(shapes, z_w1)
    s_map = inj.base["s"]
    s_max = float(np.nanmax(s_map))
    iy, ix = np.unravel_index(np.nanargmax(s_map), s_map.shape)
    base_null = inj.base["null_max"]
    seeds = np.random.SeedSequence(args.seed).spawn(len(args.mass))
    trials, eff = [], {}
    for mass, seed in zip(args.mass, seeds, strict=True):
        rng = np.random.default_rng(seed)
        theta_row = np.where(inj.background, ir.theta_e_arcsec(mass, z_w1, z_src), np.nan)
        theta_2 = float(ir.theta_e_arcsec(mass, z_w1, ir.Z_SOURCE)[0])
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
        "model": spec["model"],
        "z_lens": model.z_lens,
        "catalog": str(cat),
        "photoz": str(photoz) if photoz else None,
        "screen": inj.base["counts"]
        | {
            "psf_sigma_px": psf,
            "n_streak_vetoed": inj.n_streak,
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
            **({} if getattr(args, "lens_z", None) is None else {"w1_lens_z": z_w1}),
            **({} if r_fixed is None else {"fixed_responsivity": r_fixed}),
            **(
                {}
                if streak is None
                else {
                    "streak_veto": {"veto_arcsec": wsm.VETO_ARCSEC, "masks": "w1_streak_mask.py"},
                    "streak_mask": wsm.mask_identity(
                        args.streak_masks / f"{name}_streak_mask.fits.gz"
                    ),
                }
            ),
            "screen": SHEAR_DEFAULTS
            | {
                k: getattr(sargs, k)
                for k in ("aperture_arcsec", "filter", "r_min_arcsec", "min_coverage")
            }
            | {"max_radius": sargs.max_radius}
            | {"psf_chi": [sargs.psf_chi.real, sargs.psf_chi.imag]},
        },
    }
    (out / "summary.json").write_text(json.dumps(es.json_safe(summary), indent=1) + "\n")
    return summary


def streak_veto(mask_dir: Path, name: str):
    """``(ra, dec) -> bool`` for the field's streak mask (``w1_streak_mask.py``)."""
    path = mask_dir / f"{name}_streak_mask.fits.gz"
    if not path.exists():
        raise FileNotFoundError(f"no streak mask for {name}: {path}")
    mask, wcs = wsm.read_mask(path)
    return wsm.veto_function(mask, wcs)


def setting(summary: dict) -> dict:
    """A summary's screen settings that must agree across combined fields: the grid half-width and
    star table are per field; options missing from older summaries take their legacy value."""
    a = json.loads(json.dumps(summary["assumptions"]))
    a["screen"].pop("max_radius", None)
    a.pop("spike_stars", None)
    a.pop("streak_mask", None)  # per field
    a.pop("fixed_responsivity", None)  # per field kind (D-098)
    a["screen"] = {**LEGACY_SCREEN_OPTIONS, **a["screen"]}
    return a


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fields", nargs="+", default=list(DEFAULT_FIELDS), choices=FIELD_NAMES)
    ap.add_argument("--mass", nargs="+", type=float, default=list(ir.MASSES), help="|M| in M_sun")
    ap.add_argument("--n-inject", type=int, default=200, help="lenses per field and mass")
    ap.add_argument("--recover-tol", type=float, default=ir.RECOVER_TOL)
    ap.add_argument("--aperture-arcsec", type=float, default=SHEAR_DEFAULTS["aperture_arcsec"])
    ap.add_argument(
        "--filter", default=SHEAR_DEFAULTS["filter"], choices=("pointmass", "schirmer", "tophat")
    )
    ap.add_argument("--r-min-arcsec", type=float, default=SHEAR_DEFAULTS["r_min_arcsec"])
    ap.add_argument(
        "--min-coverage",
        type=float,
        default=SHEAR_DEFAULTS["min_coverage"],
        help="drop grid centres with less of their aperture inside the footprint (D-053 cut)",
    )
    ap.add_argument(
        "--psf-chi",
        type=float,
        nargs=2,
        default=(0.0, 0.0),
        metavar=("CHI1", "CHI2"),
        help="constant PSF distortion for the shape deconvolution (default isotropic; "
        "scripts/w1_psf_floor.py)",
    )
    ap.add_argument(
        "--streak-masks",
        type=Path,
        default=None,
        help="directory of <field>_streak_mask.fits.gz (scripts/w1_streak_mask.py): drop shape "
        "rows and painted images on streaks (results/w1_streak_veto/PLAN.md)",
    )
    ap.add_argument(
        "--lens-z",
        type=float,
        default=None,
        help="redshift of the injected W1 lenses (default: each field's cluster redshift; "
        "results/w1_volume/PLAN.md)",
    )
    ap.add_argument("--seed", type=int, default=50)
    ap.add_argument("--out", type=Path, default=paths.outputs_dir() / "inject_shear")
    ap.add_argument(
        "--combine-only",
        action="store_true",
        help="combine existing <out>/<field>/summary.json files (fields run as separate processes)",
    )
    args = ap.parse_args(argv)
    if tuple(args.psf_chi) != (0.0, 0.0) and len(args.fields) > 1:
        ap.error("--psf-chi is measured per field; run one field at a time")
    keys = [f"{m:.0e}" for m in args.mass]
    if len(set(keys)) != len(keys):
        ap.error(f"--mass values must differ at one significant figure (keys {keys})")
    if args.combine_only:
        summaries = [json.loads((args.out / f / "summary.json").read_text()) for f in args.fields]
    else:
        summaries = [run_field(f, args) for f in args.fields]
    combined = ir.combine(summaries, args.mass, {f: field_spec(f) for f in args.fields})

    if len({json.dumps(setting(s), sort_keys=True) for s in summaries}) != 1:
        ap.error("the fields were run with different screen settings; re-run them alike")
    combined["screen"] = setting(summaries[0])
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "combined.json").write_text(json.dumps(es.json_safe(combined), indent=1) + "\n")
    print(json.dumps(combined, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
