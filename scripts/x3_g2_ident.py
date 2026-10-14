"""X3 G2 prerequisite: can Gaia astrometry separate W3 from an ordinary blended variable?

Unit file: docs/state/units/X3.md (G2 prerequisite). Analytic model plus a few representative cases;
everything here is `simulated` or `derived` (no observation is tested).

Photocentre of a source (flux F_s A(t), at the origin) blended with a constant star (F_b, at offset
d) when a lens deflects the source images:

    x_c(t) = [F_s A(t) (x_src + theta_E delta(t)) + F_b d] / F(t),   F(t) = F_s A(t) + F_b

  term 1 (blend ratio):  d F_b / F(t)                        -- depends on the observed flux only
  term 2 (lens images):  theta_E delta(t) F_s A(t) / F(t)     -- the image centroid shift, weighted
                                                                by the lensed share of the light
For the repulsive (negative-mass, n = 1) lens beta = theta + 1/theta (theta_E units), images exist
only for u >= 2, both on the source's side: theta = (u +- sqrt(u^2 - 4)) / 2, magnification
|1 / ((1 + theta^-2)(1 - theta^-2))|; delta is the flux-weighted image centroid minus the source.

Ordinary models given the same blend freedom:
  O1  one variable component and one constant one (either is the variable): x_c is affine in 1/F(t)
      along one fixed direction, so term 1 is reproduced exactly. Fitted basis (along-scan, AL):
      position (2), proper motion (2), parallax (1), blend vector times 1/F (2).
  O2  two variable components with free light curves: x_c = x_2 + n g(t), g = F_1 / F arbitrary in
      [0, 1]. With one AL measurement per transit this absorbs any AL signal unless g is
      constrained (colour, smoothness); treated analytically, not fitted.

Score: dchi2 = sum((AL residual after the weighted least-squares fit) / sigma)^2 for noise-free
data, i.e. the expected chi^2 gain of the W3 truth over the best O1 fit. Criterion (X3.md, before
computing): median >= 25 measurable, 9-25 marginal, < 9 not measurable.

ASSUMPTIONS (stated, not measured): scan angles uniform per transit; DR4 transit times = the DR3
times of one dev source repeated after 1,038 d up to 2020-01-20; per-transit AL sigma from dev
catalogue errors (max(ra_error, dec_error) * sqrt(n_transits / 2), n_transits =
astrometric_n_good_obs_al / 9); no microlens parallax; source and blend at M31's position for the
parallax factors; finite source rho = 0.01 by disc sampling; transits fainter than G 20.7 carry no
astrometry; AL sigma scales as F^-1/2 with the total flux.

  python scripts/x3_g2_ident.py        # writes results/x3/g2_ident.json
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


x3 = _load("x3_gaps")

DR3_END = 2457901.5  # 2017-05-28
DR4_END = 2458869.4  # 2020-01-20
DR3_SPAN_D = 1038.0
N_DRAWS = 25
DET_LIMIT_G = 20.7  # per-transit detection limit, the G1 ASSUMPTION (x3_g1.DET_LIMIT)
REF_SOURCE = "370857476684663040"  # dev source whose DR3 G transit times set the cadence
BASE = {
    "G": 17.0,
    "theta_e_mas": 3.0,
    "fs": 0.7,
    "u0": 0.5,
    "te_d": 100.0,
    "rho": 0.01,
    "d_mas": 200.0,
    "data": "DR4",
}
VARIATIONS = {
    "G": (15.0, 19.0),
    "theta_e_mas": (1.0, 10.0),
    "fs": (1.0, 0.5),
    "u0": (0.1, 1.5),
    "te_d": (30.0, 300.0),
    "d_mas": (20.0, 1000.0),
    "data": ("DR3",),
}


# ----------------------------------------------------------------------------- lens model
def images_neg(b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Image positions (theta_E units, along the source direction) and |magnifications|."""
    disc = b * b - 4.0
    s = np.sqrt(np.clip(disc, 0.0, None))
    th = np.stack([(b + s) / 2.0, (b - s) / 2.0])
    with np.errstate(divide="ignore", invalid="ignore"):
        mu = np.abs(1.0 / ((1.0 + th**-2) * (1.0 - th**-2)))
    mu = np.where(disc >= 0, np.nan_to_num(mu, posinf=0.0), 0.0)
    return th, mu


def lensed(ux: np.ndarray, uy: np.ndarray, rho: float, n: int = 300, seed: int = 0):
    """Magnification A and image-centroid shift (dx, dy) in theta_E units for each epoch."""
    rng = np.random.default_rng(seed)
    r = rho * np.sqrt(rng.uniform(0.0, 1.0, n))
    a = rng.uniform(0.0, 2.0 * np.pi, n)
    x = ux[:, None] + r * np.cos(a)
    y = uy[:, None] + r * np.sin(a)
    b = np.hypot(x, y)
    th, mu = images_neg(b)
    amp = mu.sum(0)
    w = (mu * th).sum(0) / np.where(b > 0, b, 1.0)  # centroid along (x, y)/b, times b
    A = amp.mean(1)
    cx = (w * x).sum(1) / np.maximum(amp.sum(1), 1e-30)
    cy = (w * y).sum(1) / np.maximum(amp.sum(1), 1e-30)
    return A, cx - ux, cy - uy  # shift relative to the source-disc centre


# ----------------------------------------------------------------------------- Gaia-like data
def transit_times(data: str) -> np.ndarray:
    from jwst_anomaly import gaia_mulens as gm

    text = (x3.raw_dir() / f"{REF_SOURCE}.csv").read_text()
    t = np.sort(np.asarray(gm.parse_epoch_csv(text)["time"], float))
    if data == "DR4":
        t2 = t + DR3_SPAN_D
        t = np.concatenate([t, t2[t2 < DR4_END]])
    return t


def parallax_al(t: np.ndarray, psi: np.ndarray) -> np.ndarray:
    import astropy.units as u
    from astropy.coordinates import get_body_barycentric
    from astropy.time import Time

    ra, dec = np.radians(x3.M31[0]), np.radians(x3.M31[1])
    e = get_body_barycentric("earth", Time(t, format="jd", scale="tcb"))
    X, Y, Z = (e.x.to(u.au).value, e.y.to(u.au).value, e.z.to(u.au).value)
    pa = X * np.sin(ra) - Y * np.cos(ra)
    pd = X * np.cos(ra) * np.sin(dec) + Y * np.sin(ra) * np.sin(dec) - Z * np.cos(dec)
    return pa * np.sin(psi) + pd * np.cos(psi)


def sigma_al(gmag: float) -> float:
    """Per-FoV-transit AL sigma (mas) from dev catalogue errors, median in a 1-mag bin."""
    from astropy.table import Table

    at = Table.read(x3.out_dir() / "g1_sel_attrs.csv", format="ascii.csv")
    ids = x3.ids_table(19.0)
    gm_ = {str(s): float(g) for s, g in zip(ids["source_id"], ids["phot_g_mean_mag"], strict=True)}
    vals, gs = [], []
    for r in at:
        s = str(r["source_id"])
        if s not in gm_ or np.ma.is_masked(r["ra_error"]):
            continue
        n_tr = float(r["astrometric_n_good_obs_al"]) / 9.0
        vals.append(max(float(r["ra_error"]), float(r["dec_error"])) * np.sqrt(n_tr / 2.0))
        gs.append(gm_[s])
    vals, gs = np.array(vals), np.array(gs)
    sel = np.abs(gs - gmag) <= 0.5
    return float(np.median(vals[sel]))


# ----------------------------------------------------------------------------- one case
def case(p: dict, sig: float, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    t = transit_times(p["data"])
    t0 = 2458300.0 if p["data"] == "DR4" else 2457400.0  # 2018-07 (after DR3) or 2016-01
    out = {
        "dchi2_O1": [],
        "dchi2_pos5": [],
        "term2_raw": [],
        "term1_resid": [],
        "n_lost": [],
        "dchi2_vs_pos_shape": [],
    }
    for _ in range(N_DRAWS):
        psi = rng.uniform(0.0, np.pi, t.size)
        phi_l = rng.uniform(0.0, 2.0 * np.pi)  # lens track direction on the sky
        phi_b = rng.uniform(0.0, 2.0 * np.pi)  # blend direction
        tau = (t - t0) / p["te_d"]
        ux0, uy0 = tau, np.full_like(tau, p["u0"])
        A, dx, dy = lensed(ux0, uy0, p["rho"], seed=int(rng.integers(1 << 30)))
        A = np.where(np.isfinite(A), A, 0.0)
        c, s = np.cos(phi_l), np.sin(phi_l)
        sx, sy = c * dx - s * dy, s * dx + c * dy  # image shift on the sky (theta_E units)
        fb = 1.0 - p["fs"]
        F = np.maximum(p["fs"] * A + fb, 1e-12)  # F -> 0 only below the detection limit (masked)
        w2 = p["fs"] * A / F
        bx, by = p["d_mas"] * np.cos(phi_b), p["d_mas"] * np.sin(phi_b)
        t1x, t1y = bx * fb / F, by * fb / F
        t2x, t2y = p["theta_e_mas"] * sx * w2, p["theta_e_mas"] * sy * w2
        # transits fainter than the detection limit carry no astrometry; sigma grows as F^-1/2
        keep = p["G"] - 2.5 * np.log10(np.maximum(F, 1e-12)) <= DET_LIMIT_G
        sig_i = sig / np.sqrt(np.maximum(F, 1e-12))
        sp, cp = np.sin(psi), np.cos(psi)
        al1 = t1x * sp + t1y * cp
        al2 = t2x * sp + t2y * cp
        ty = (t - t.mean()) / 365.25
        b5 = np.vstack([sp, cp, sp * ty, cp * ty, parallax_al(t, psi)]).T
        bo1 = np.hstack([b5, np.vstack([sp / F, cp / F]).T])

        def resid(y, B, keep=keep, w=1.0 / sig_i[keep]):
            coef = np.linalg.lstsq(B[keep] * w[:, None], y[keep] * w, rcond=None)[0]
            return float(np.sum(((y[keep] - B[keep] @ coef) * w) ** 2))

        out["dchi2_O1"].append(resid(al1 + al2, bo1))
        # positive-mass image-centroid shape on the same |u(t)|, free theta_E (linear), track
        # direction and side (grid): what astrometry adds beyond "a lens is present"
        best = np.inf
        for side in (1.0, -1.0):
            px, py = ux0, side * uy0
            ppos = px * px + py * py + 2.0
            qx, qy = px / ppos, py / ppos  # PSPL centroid shift (theta_E units)
            for ang in np.linspace(0.0, 2.0 * np.pi, 72, endpoint=False):
                ca, sa = np.cos(ang), np.sin(ang)
                col = ((ca * qx - sa * qy) * sp + (sa * qx + ca * qy) * cp) * w2
                best = min(best, resid(al1 + al2, np.hstack([bo1, col[:, None]])))
        out["dchi2_vs_pos_shape"].append(best)
        out["dchi2_pos5"].append(resid(al1 + al2, b5))
        out["term2_raw"].append(float(np.sum((al2[keep] / sig_i[keep]) ** 2)))
        out["n_lost"].append(int(np.sum(~keep)))
        out["term1_resid"].append(resid(al1, bo1))
    med = {k: float(np.median(v)) for k, v in out.items()}
    v = med["dchi2_O1"]
    med["verdict"] = "measurable" if v >= 25 else ("marginal" if v >= 9 else "not measurable")
    med["n_transits"] = int(t.size)
    med["n_transits_with_source_light_lensed"] = int(
        np.sum(np.abs((t - (2458300.0 if p["data"] == "DR4" else 2457400.0)) / p["te_d"]) < 5)
    )
    med["sigma_al_mas"] = sig
    return med


def main() -> int:
    cases = [("base", dict(BASE))]
    for k, vals in VARIATIONS.items():
        for v in vals:
            cases.append((f"{k}={v}", {**BASE, k: v}))
    sig = {g: sigma_al(g) for g in (15.0, 17.0, 19.0)}
    rows = []
    for name, p in cases:
        r = case(p, sig[p["G"]], seed=11)
        rows.append({"case": name, **{k: p[k] for k in BASE}, **r})
        print(
            f"{name:16s} sigma={r['sigma_al_mas']:.2f} mas  O1 dchi2={r['dchi2_O1']:9.1f}  "
            f"pos5-only={r['dchi2_pos5']:10.1f}  term2 raw={r['term2_raw']:9.1f}  "
            f"term1 left={r['term1_resid']:.1e}  vs +mass shape={r['dchi2_vs_pos_shape']:7.1f}  "
            f"lost={r['n_lost']:.0f}/{r['n_transits']}  {r['verdict']}"
        )
    rec = {
        "label": "simulated / derived identifiability check; no observation tested",
        "criterion": "median dchi2_O1 >= 25 measurable, 9-25 marginal, < 9 not measurable",
        "assumptions": __doc__.split("ASSUMPTIONS")[1].split("python scripts")[0].strip(),
        "sigma_al_mas_by_G": sig,
        "n_draws": N_DRAWS,
        "cases": rows,
    }
    out = x3.results_dir() / "g2_ident.json"
    out.write_text(json.dumps(rec, indent=1) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
