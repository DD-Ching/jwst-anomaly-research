"""Injection-recovery of dark-deflector image pairs through the orphan-pair screen (D-051).

Simulated lenses are painted into the real deep-field catalogues and the unchanged orphan-pair
rules (``orphan_pairs.select_sources``, ``match_table``, ``classify_pairs``) are run on the
result. Three lens types (exotic_sim, D-047):

- ``point``: dark point mass, n = 1, eps > 0 (sanity check; an ordinary compact dark lens);
- ``W2``: Ellis wormhole, n = 2, eps > 0, no visible deflector;
- ``W1``: negative mass, n = 1, eps < 0: two radially stretched images on the source's side,
  sources at beta < 2 theta_E removed (umbra), lensed sources at beta in [2, 4] theta_E.

Every injected image is ``simulated``; efficiencies and limits are ``derived`` from them;
masses and throat radii are a ``model_prediction`` of the lens formulae. Nothing here is a
detection. Every threshold is an ASSUMPTION (constants below; docs/exotic_limits.md "W2").

Painting (per lensed catalogue row; the row itself is removed):

- image positions and signed magnifications from ``exotic_sim.inject_images``;
- fluxes: ``|mu| x`` the row's (noisy) SED, plus fresh Gaussian noise of sqrt(max(1 - mu^2, 0))
  times the row's errors, so the scatter is max(|mu|, 1) x the row's errors (the row's own noise
  is scaled with it); errors = the row's errors (sky-limited);
- photo-z (z_low, z16, z84, z_best) and the catalogue mu are inherited from the row;
- size: Kron aperture radius ``sqrt(r_psf^2 + (r^2 - r_psf^2) s^2)``, with s the larger
  stretch 1/|lambda| and r_psf the 5th percentile of the selected sources' radii;
- deblending: two images closer than ``d_blend`` (the 2nd percentile of catalogue
  nearest-neighbour separations) become one row at the flux-weighted centroid with the summed
  flux; such a lens can never be recovered as a pair.

A lensed source is **recovered** when its two image rows both pass the source selection, are
0.3-3" apart, match in SED and photo-z, and are classified ``orphan``.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))

import orphan_pairs as op  # noqa: E402

from jwst_anomaly import exotic_sim, paths, schema  # noqa: E402

LENSES = {
    "point": {"n": 1.0, "sign": 1, "beta": (0.0, 3.0), "umbra": 0.0},
    "W2": {"n": 2.0, "sign": 1, "beta": (0.0, 3.0), "umbra": 0.0},
    "W1": {"n": 1.0, "sign": -1, "beta": (2.0, 4.0), "umbra": 2.0},
}
# beta_max 3 theta_E for attractive lenses: beyond it the inner image of a point mass carries
# < 1.3 % of the flux and cannot pass S/N >= 10 in 8 bands (ASSUMPTION; conservative).
THETA_E = (0.15, 0.3, 0.7, 1.5)  # arcsec (pair separations ~2 theta_E inside 0.3-3")
DEEP = ("macs0416-ncf", "macs1149-ncf", "abell370-ncf", "macs0417-ncf", "macs1423-ncf")
DEEP += ("goodsn-dja",)
N_LENS = 2000  # random-position lenses per field, type and theta_E (per-deflector efficiency)
N_SOURCE = 400  # lensed catalogue rows per field, type and theta_E (per-source efficiency)
BLEND_PERCENTILE = 2.0  # d_blend: this percentile of nearest-neighbour separations (ASSUMPTION)
BLEND_MAX = 0.5  # arcsec cap on d_blend (sparse catalogues; ASSUMPTION)
PSF_PERCENTILE = 5.0  # r_psf: this percentile of selected sources' aperture radii (ASSUMPTION)
LOCAL_PAD = 6.0  # arcsec around the lensed region evaluated per trial (> 0.5 x 3" + 0.3")
MAG_EDGES = (24.0, 25.0, 26.0, 27.0, 28.0)
POISSON_UL_95 = 2.996  # -ln(0.05): 95 % upper limit on a Poisson mean with zero events
Z_SOURCE = 2.0  # source redshift for theta_E -> mass / throat radius (ASSUMPTION; D-049)


# ----------------------------------------------------------------------------- physics


def _distances(z_l: float, z_s: float = Z_SOURCE):
    from astropy.cosmology import Planck18

    return (
        Planck18.angular_diameter_distance(z_l).si.value,
        Planck18.angular_diameter_distance(z_s).si.value,
        Planck18.angular_diameter_distance(z_l, z_s).si.value,
    )


def eps_bar_of(theta_e_arcsec, n: float, z_l: float, z_s: float = Z_SOURCE) -> np.ndarray:
    """Inverse of ``exotic_sim.einstein_radius``: |eps_bar| (m^n) for an angular theta_E."""
    d_l, d_s, d_ls = _distances(z_l, z_s)
    theta = np.deg2rad(np.asarray(theta_e_arcsec, float) / 3600.0)
    return theta ** (n + 1) * d_s * d_l**n / d_ls


def theta_e_to_mass(theta_e_arcsec, z_l: float, z_s: float = Z_SOURCE) -> np.ndarray:
    """|M| (M_sun) of an n = 1 lens: eps_bar = 4 G |M| / c^2 (model_prediction; Planck18)."""
    eps = eps_bar_of(theta_e_arcsec, 1.0, z_l, z_s)
    return eps * exotic_sim.C_SI**2 / (4.0 * exotic_sim.G_SI * exotic_sim.MSUN_KG)


def theta_e_to_throat_pc(theta_e_arcsec, z_l: float, z_s: float = Z_SOURCE) -> np.ndarray:
    """Ellis throat radius a (pc): eps_bar = pi a^2 / 4 (model_prediction; Planck18)."""
    eps = eps_bar_of(theta_e_arcsec, 2.0, z_l, z_s)
    return np.sqrt(4.0 * eps / np.pi) / exotic_sim.PC_M


def surface_density_limit(exposure_deg2: float, n_ul: float = POISSON_UL_95) -> float:
    """95 % upper limit (deg^-2): n_ul / sum(efficiency x area)."""
    return n_ul / exposure_deg2 if exposure_deg2 > 0 else float("inf")


def poisson_signal_ul(n_obs: int, background: float, cl: float = 0.95) -> float:
    """CLs upper limit s on a Poisson signal over a known background b (Read 2002): the s with
    P(N <= n_obs | s + b) / P(N <= n_obs | b) = 1 - cl. Unlike the classical limit it never
    collapses to 0 (excluding every density) when fewer events than b are observed."""
    from scipy.optimize import brentq
    from scipy.stats import poisson

    if not (np.isfinite(background) and background >= 0):
        raise ValueError(f"poisson_signal_ul: background must be finite and >= 0, not {background}")
    p_b = poisson.cdf(n_obs, background)
    f = lambda s: poisson.cdf(n_obs, s + background) / p_b - (1.0 - cl)  # noqa: E731
    hi = max(10.0, 3.0 * (n_obs + 10))
    return float(brentq(f, 0.0, hi))


# ----------------------------------------------------------------------------- field context


class Field:
    """A deep-field catalogue prepared for injections."""

    def __init__(self, name: str, std: Table | None = None, z_lens: float | None = None):
        self.name = name
        self.std = std if std is not None else op.load_catalogue(name)
        spec = op.FIELDS.get(name, {})
        self.z_lens = z_lens if z_lens is not None else op.z_lens_of(spec)
        gx, gy, area, (ra0, dec0) = op.searched_footprint(self.std)
        self.grid = np.column_stack([gx, gy])
        self.area_arcsec2 = area
        self.ra0, self.dec0 = ra0, dec0
        x, y = op.tangent_xy(self.std["ra"], self.std["dec"], ra0, dec0)
        good = np.isfinite(x) & np.isfinite(y)
        self.x, self.y = np.where(good, x, 1e9), np.where(good, y, 1e9)
        self.tree = cKDTree(np.column_stack([self.x, self.y]))
        nn = self.tree.query(np.column_stack([self.x[good], self.y[good]]), k=2)[0][:, 1]
        self.d_blend = min(float(np.percentile(nn, BLEND_PERCENTILE)), BLEND_MAX)
        sel = op.select_sources(self.std, self.z_lens)
        rad = np.asarray(self.std["ap_radius"], float)
        self.r_psf = float(np.nanpercentile(rad[sel], PSF_PERCENTILE)) if sel.any() else 0.1
        flux = np.asarray(self.std["flux"], float)
        self.lensable = (
            (np.asarray(self.std["z_best"], float) > self.z_lens) & np.isfinite(flux).any(1) & good
        )
        lw, _ = op.flux_matrix(self.std, op.snr_bands_of(self.std))
        with np.errstate(invalid="ignore", divide="ignore"):
            self.mag = 31.4 - 2.5 * np.log10(np.nanmean(lw, 1))  # AB, nJy; mean of S/N bands
        self.selected = sel
        self.n_selected = int(sel.sum())

    def to_radec(self, x, y):
        dec = self.dec0 + np.asarray(y) / 3600.0
        ra = self.ra0 + np.asarray(x) / 3600.0 / np.cos(np.deg2rad(self.dec0))
        return ra, dec


# ----------------------------------------------------------------------------- painting


def paint(field: Field, lx: float, ly: float, theta_e: float, ltype: str, rng, rows=None):
    """Lens catalogue rows around (lx, ly) arcsec: ``(removed_rows, images, lensed_rows)``.

    ``images`` is a standard-layout table with ``inj_src`` (catalogue row of the source) and
    ``inj_img`` (image index; 2 = merged pair). ``rows`` restricts lensing to those rows
    (default: every lensable row within beta_max theta_E)."""
    spec = LENSES[ltype]
    bmin, bmax = spec["beta"]
    near = np.asarray(field.tree.query_ball_point([lx, ly], bmax * theta_e), int)
    near = near[field.lensable[near]]
    if rows is not None:
        near = np.intersect1d(near, np.asarray(rows, int))
    dx, dy = field.x[near] - lx, field.y[near] - ly
    beta = np.hypot(dx, dy) / theta_e
    umbra = near[beta < spec["umbra"]]
    use = (beta >= bmin) & (beta <= bmax) & (beta >= spec["umbra"])
    src = near[use]
    removed = np.concatenate([umbra, src])
    if len(src) == 0:
        return removed, None, src
    imgs = exotic_sim.inject_images(
        dx[use], dy[use], theta_e, spec["n"], spec["sign"], source_id=src
    )
    out_rows, inj_src, inj_img = [], [], []
    std = field.std
    flux, err = np.asarray(std["flux"], float), np.asarray(std["err"], float)
    rad = np.asarray(std["ap_radius"], float)
    for s in src:
        im = imgs[np.asarray(imgs["source_id"]) == s]
        mus = np.abs(np.asarray(im["mu"], float))
        ix, iy = lx + np.asarray(im["dx"], float), ly + np.asarray(im["dy"], float)
        stretch = np.fmax(np.asarray(im["stretch_t"], float), np.asarray(im["stretch_r"], float))
        r0 = rad[s]
        r_img = np.sqrt(field.r_psf**2 + max(r0**2 - field.r_psf**2, 0.0) * stretch**2)
        sep = np.hypot(ix[1] - ix[0], iy[1] - iy[0]) if len(im) == 2 else np.inf
        if sep < field.d_blend:
            w = mus / mus.sum()  # merged: one row at the flux-weighted centroid
            ix, iy = np.array([np.dot(w, ix)]), np.array([np.dot(w, iy)])
            r_img = np.array([r_img.max() + 0.5 * sep])
            mus, idx = np.array([mus.sum()]), [2]
        else:
            idx = list(np.asarray(im["image"], int))
        for k in range(len(mus)):
            # the catalogue flux already carries its noise, scaled here by mu; add fresh noise only
            # up to the image's own (sky-limited) error, so the scatter is max(mu, 1) sigma
            extra = np.sqrt(max(1.0 - mus[k] ** 2, 0.0))
            f = mus[k] * flux[s] + extra * rng.normal(0.0, 1.0, flux.shape[1]) * np.nan_to_num(
                err[s]
            )
            out_rows.append((ix[k], iy[k], f, r_img[k]))
            inj_src.append(int(s))
            inj_img.append(int(idx[k]))
    if not out_rows:
        return removed, None, src
    base = std[np.asarray(inj_src, int)]
    base = Table(base, copy=True)
    ra, dec = field.to_radec([r[0] for r in out_rows], [r[1] for r in out_rows])
    base["ra"], base["dec"] = ra, dec
    f = np.vstack([r[2] for r in out_rows])
    f[~np.isfinite(np.asarray(base["err"], float))] = np.nan  # invalid bands stay invalid
    base["flux"] = f
    base["ap_radius"] = [r[3] for r in out_rows]
    for c in ("x_pix", "y_pix", "xmin", "xmax", "ymin", "ymax", "pa_deg"):
        base[c] = np.nan
    base["deblend"] = False
    base["src_id"] = -(np.asarray(inj_src, int) * 10 + np.asarray(inj_img, int) + 1)
    base["inj_src"], base["inj_img"] = inj_src, inj_img
    base.meta.update(
        provenance=schema.Provenance.SIMULATED.value,
        source=f"exotic_sim images painted into {base.meta.get('source', 'the catalogue')}",
    )
    return removed, base, src


def evaluate(field: Field, lx: float, ly: float, reach: float, removed, images) -> dict[int, dict]:
    """Run the orphan-pair rules on the painted neighbourhood; per lensed source, the outcome."""
    out: dict[int, dict] = {}
    if images is None:
        return out
    local = np.asarray(field.tree.query_ball_point([lx, ly], reach + LOCAL_PAD), int)
    local = np.setdiff1d(local, removed)
    real = field.std[local]
    real["inj_src"], real["inj_img"] = -1, -1
    cat = vstack([real, images], join_type="exact", metadata_conflicts="silent")
    cat.meta = dict(field.std.meta)
    inj_src = np.asarray(cat["inj_src"], int)
    for s in np.unique(inj_src[inj_src >= 0]):
        k = np.flatnonzero(inj_src == s)
        out[int(s)] = {"n_rows": len(k), "n_selected": 0, "matched": False, "cls": "", "sep": 0.0}
    sel = op.select_sources(cat, field.z_lens)
    for s in out:
        out[s]["n_selected"] = int(sel[inj_src == s].sum())
    rows = np.flatnonzero(sel)
    pair_src = [s for s, o in out.items() if o["n_selected"] >= 2]
    if not pair_src:
        return out
    x, y = op.tangent_xy(cat["ra"][rows], cat["dec"][rows], field.ra0, field.dec0)
    i, j, d = op.close_pairs(x, y, op.SEP_MIN, op.SEP_MAX)
    si, sj = inj_src[rows[i]], inj_src[rows[j]]
    same = (si >= 0) & (si == sj)
    if not same.any():
        return out
    sub = cat[rows]
    f, e = np.asarray(sub["flux"], float), np.asarray(sub["err"], float)
    z16, z84 = np.asarray(sub["z16"], float), np.asarray(sub["z84"], float)
    m = op.match_table(i[same], j[same], f, e, z16, z84)
    m["sep"] = d[same]
    m["ci"], m["cj"] = rows[m["i"]], rows[m["j"]]
    empty = Table({"image_id": [], "ra": [], "dec": []}, dtype=[str, float, float])
    m = op.classify_pairs(m, cat, empty, field.ra0, field.dec0)
    for r in m:
        s = int(inj_src[r["ci"]])
        o = out[s]
        o["sep"] = float(r["sep"])
        o["cls"] = str(r["pair_class"])
        o["matched"] = bool(r["match"])
        o["recovered"] = bool(r["match"]) and r["pair_class"] == "orphan"
    return out


# ----------------------------------------------------------------------------- experiments


def per_source(field: Field, ltype: str, theta_e: float, n: int, rng) -> Table:
    """Lens ``n`` random rows that the search selects unlensed (behind z_lens, S/N cuts), each at
    a random beta (uniform in area); per-source outcomes. Fainter rows that magnification lifts
    into the selection are counted by ``per_lens``."""
    bmin, bmax = LENSES[ltype]["beta"]
    bmin = max(bmin, LENSES[ltype]["umbra"])
    cand = np.flatnonzero(field.lensable & field.selected)
    rows = []
    for s in rng.choice(cand, size=min(n, len(cand)), replace=False):
        beta = np.sqrt(rng.uniform(bmin**2, bmax**2)) * theta_e
        phi = rng.uniform(0, 2 * np.pi)
        lx, ly = field.x[s] - beta * np.cos(phi), field.y[s] - beta * np.sin(phi)
        removed, images, _ = paint(field, lx, ly, theta_e, ltype, rng)
        res = evaluate(field, lx, ly, bmax * theta_e, removed, images).get(int(s), {})
        rows.append(
            (
                int(s),
                beta / theta_e,
                float(field.mag[s]),
                int(res.get("n_rows", 0)),
                int(res.get("n_selected", 0)),
                bool(res.get("matched", False)),
                str(res.get("cls", "")),
                bool(res.get("recovered", False)),
                float(res.get("sep", 0.0)),
            )
        )
    names = ("row", "beta", "mag", "n_images", "n_selected", "matched", "cls", "recovered", "sep")
    t = Table(rows=rows, names=names) if rows else Table(names=names)
    t.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"orphan-pair screen on {field.name} with injected {ltype} pairs",
        lens=ltype,
        theta_e=theta_e,
    )
    return t


def per_lens(field: Field, ltype: str, theta_e: float, n: int, rng) -> dict:
    """``n`` lenses at uniform random positions in the searched footprint: recovered fraction."""
    bmax = LENSES[ltype]["beta"][1]
    pick = rng.integers(0, len(field.grid), n)
    pos = field.grid[pick] + rng.uniform(-0.5, 0.5, (n, 2))
    n_rec = n_lensed = n_any = 0
    for lx, ly in pos:
        removed, images, src = paint(field, lx, ly, theta_e, ltype, rng)
        n_lensed += len(src)
        if images is None:
            continue
        n_any += 1
        res = evaluate(field, lx, ly, bmax * theta_e, removed, images)
        n_rec += any(o.get("recovered", False) for o in res.values())
    return {
        "n_lens": n,
        "n_recovered": int(n_rec),
        "efficiency": n_rec / n,
        "mean_lensed_sources": n_lensed / n,
        "frac_with_images": n_any / n,
    }


def mag_table(t: Table) -> list[dict]:
    """Recovered fraction per source-magnitude bin."""
    edges = (-np.inf, *MAG_EDGES, np.inf)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        s = (np.asarray(t["mag"]) >= lo) & (np.asarray(t["mag"]) < hi)
        n = int(s.sum())
        out.append(
            {
                "mag": [lo, hi],
                "n": n,
                "recovered": int(np.asarray(t["recovered"], bool)[s].sum()) if n else 0,
            }
        )
    return out


def run_field(name: str, args) -> dict:
    t0 = time.time()
    field = Field(name)
    rng = np.random.default_rng(args.seed + sum(map(ord, name)))
    res = {
        "field": name,
        "area_arcmin2": field.area_arcsec2 / 3600.0,
        "d_blend_arcsec": field.d_blend,
        "r_psf_arcsec": field.r_psf,
        "n_selected": field.n_selected,
        "z_lens": field.z_lens,
        "runs": {},
    }
    out_dir = args.out / name
    out_dir.mkdir(parents=True, exist_ok=True)
    for ltype in args.lenses:
        for te in args.theta_e:
            src = per_source(field, ltype, te, args.n_source, rng)
            src.write(out_dir / f"per_source_{ltype}_{te}.ecsv", overwrite=True)
            lens = per_lens(field, ltype, te, args.n_lens, rng)
            cls, counts = np.unique(np.asarray(src["cls"], str), return_counts=True)
            res["runs"][f"{ltype}:{te}"] = {
                "lens": ltype,
                "theta_e": te,
                "per_source": {
                    "n": len(src),
                    "recovered": int(np.sum(src["recovered"])),
                    "two_rows": int(np.sum(np.asarray(src["n_images"]) >= 2)),
                    "both_selected": int(np.sum(np.asarray(src["n_selected"]) >= 2)),
                    "matched": int(np.sum(src["matched"])),
                    "classes": {str(c): int(k) for c, k in zip(cls, counts, strict=True)},
                    "by_mag": mag_table(src),
                },
                "per_lens": lens,
            }
            print(
                f"{name} {ltype} theta_E={te}: per-source {int(np.sum(src['recovered']))}/"
                f"{len(src)}, per-lens {lens['n_recovered']}/{lens['n_lens']}",
                flush=True,
            )
    res["wall_s"] = time.time() - t0
    (out_dir / "injection_summary.json").write_text(json.dumps(op._finite(res), indent=1))
    return res


def combine(results: list[dict], orphans: dict[str, tuple[int, float]]) -> dict:
    """Per lens type and theta_E: summed exposure, both limits, and the mass scale."""
    out = {}
    n_obs = sum(v[0] for v in orphans.values())
    bkg = sum(v[1] for v in orphans.values())
    s95 = poisson_signal_ul(n_obs, bkg)
    keys = sorted(set.intersection(*(set(r["runs"]) for r in results)))  # run in every field
    for k in keys:
        ltype, te = k.split(":")
        te = float(te)
        exposure = sum(
            r["runs"][k]["per_lens"]["efficiency"] * r["area_arcmin2"] / 3600.0 for r in results
        )
        z_l = op.Z_LENS_REF
        row = {
            "lens": ltype,
            "theta_e": te,
            "exposure_deg2": exposure,
            "limit_no_candidate_deg2": surface_density_limit(exposure),
            "limit_background_deg2": s95 / exposure if exposure > 0 else float("inf"),
            "per_source_efficiency": sum(r["runs"][k]["per_source"]["recovered"] for r in results)
            / max(sum(r["runs"][k]["per_source"]["n"] for r in results), 1),
        }
        if LENSES[ltype]["n"] == 1:
            row["mass_msun"] = float(theta_e_to_mass(te, z_l))
        else:
            row["throat_pc"] = float(theta_e_to_throat_pc(te, z_l))
        out[k] = row
    return {
        "orphans_observed": n_obs,
        "orphans_expected_null_e": bkg,
        "signal_ul_95": s95,
        "z_lens": op.Z_LENS_REF,
        "z_source": Z_SOURCE,
        "limits": out,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    # deep fields only: combine converts theta_E to mass at Z_LENS_REF, which a cluster field
    # (its own lens redshift) would silently misreport
    ap.add_argument("--fields", nargs="+", default=list(DEEP), choices=sorted(DEEP))
    ap.add_argument("--lenses", nargs="+", default=list(LENSES), choices=list(LENSES))
    ap.add_argument("--theta-e", nargs="+", type=float, default=list(THETA_E))
    ap.add_argument("--n-lens", type=int, default=N_LENS)
    ap.add_argument("--n-source", type=int, default=N_SOURCE)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=Path, default=paths.outputs_dir() / "inject_pairs")
    ap.add_argument(
        "--orphans",
        type=Path,
        default=paths.outputs_dir() / "orphan_pairs",
        help="orphan_pairs.py output root: observed orphans and null (e) per field",
    )
    ap.add_argument(
        "--combine-only",
        action="store_true",
        help="skip the injections; combine <out>/<field>/injection_summary.json (fields run in "
        "parallel processes)",
    )
    args = ap.parse_args(argv)
    if args.combine_only:
        results = [
            json.loads((args.out / f / "injection_summary.json").read_text()) for f in args.fields
        ]
    else:
        results = [run_field(f, args) for f in args.fields]
    orphans = {}
    for f in args.fields:
        p = args.orphans / f / "summary.json"
        if p.exists():
            s = json.loads(p.read_text())
            expected = s["null_e_conditioned"]["expected_by_class"]["orphan"]
            if expected is None:  # NaN in the field run (no z-overlapping far pairs)
                raise SystemExit(f"error: {f}: null (e) orphan expectation undefined")
            orphans[f] = (int(s["classes"]["orphan"]), float(expected))
    missing = [f for f in args.fields if f not in orphans]
    if missing:  # the background-aware limit needs every field's observed and expected orphans
        raise SystemExit(f"error: no orphan_pairs summary for {missing}; run orphan_pairs.py first")
    summary = combine(results, orphans)
    summary["fields"] = results
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "limits.json").write_text(json.dumps(op._finite(summary), indent=1))
    print(json.dumps(op._finite({k: v for k, v in summary.items() if k != "fields"}), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
