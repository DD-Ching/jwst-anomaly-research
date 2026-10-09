"""S2 (D-069 round-1 survivor): SN Ia Hubble residuals vs a flat-kernel foreground column.

Hypothesis (A3 P2b, docs/hypotheses/round-1/): SN brightness depends on the foreground matter along
the line of sight with a *uniform* weight in redshift, not only through the lensing kernel.
Ordinary expectation: only the lensing column matters (overdense foregrounds magnify, negative
residual), and the flat coefficient is zero once the lensing column is in the fit.

Chain: Pantheon+ (observed m_b_corr) -> residual vs flat LCDM (model_prediction, Om = 0.334) ->
Legacy Surveys DR9 galaxies with photo-z inside a disc around each SN (observed) -> per-SN weighted
counts relative to the expected counts from the mean shell densities (derived) -> weighted least
squares on [1, z, X_lens, X_flat] -> the flat coefficient gamma_F with an analytic error and a
z-matched scramble null -> injection-recovery of gamma_F.

    python scripts/s2_flat_kernel.py fetch      # Data Lab TAP, batched boxes, cached per chunk
    python scripts/s2_flat_kernel.py fit        # columns, regression, scramble null, injections
"""

from __future__ import annotations

import argparse
import io
import json
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import requests
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table, vstack

from jwst_anomaly import paths

TAP = "https://datalab.noirlab.edu/tap/sync"
PPLUS_URL = (
    "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
    "Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat"
)
COLS = "t.ra, t.dec, t.mag_z, t.release, p.z_phot_median, p.z_spec"
OUT = paths.repo_root() / "results" / "s2_flat_kernel"


@dataclass
class Params:
    z_min: float = 0.1  # ASSUMPTION: below this the foreground column is too short to matter
    z_max: float = 1.3
    radius_arcsec: float = 120.0  # ASSUMPTION: disc radius of the foreground column
    mag_z_max: float = 21.0  # ASSUMPTION: photo-z reliable (Zhou et al. 2021 DR9 photo-z)
    z_gal_min: float = 0.02
    host_gap: float = 0.05  # ASSUMPTION: drop galaxies with z_g > z_s - gap (host and its group)
    shell_dz: float = 0.05
    om: float = 0.334  # Pantheon+ flat LCDM (Brout et al. 2022)
    batch: int = 40  # boxes per TAP query
    workers: int = 2  # Data Lab: batch, little concurrency
    n_scramble: int = 1000
    n_inject: int = 200
    seed: int = 69
    trim_pct: float = 99.0  # ASSUMPTION: leverage cut for the robustness fit
    min_shell_expect: float = 1.0  # ASSUMPTION: alpha = 2 uses shells expecting >= 1 galaxy
    alpha: int = 1  # flat column on (1 + delta)^alpha; A3's fiducial is alpha >= 2 (1 or 2 here)


def cache() -> Path:
    d = paths.data_root() / "s2_flat_kernel"
    d.mkdir(parents=True, exist_ok=True)
    return d


def load_pantheon(p: Params) -> Table:
    """Pantheon+ SNe Ia, Hubble-flow sample, one row per CID (inverse-variance mean)."""
    f = cache() / "PantheonPlusSH0ES.dat"
    if not f.exists():
        r = requests.get(PPLUS_URL, timeout=120)
        r.raise_for_status()
        f.write_bytes(r.content)
    t = Table.read(f, format="ascii.basic")
    t = t[(t["IS_CALIBRATOR"] == 0) & (t["zHD"] > p.z_min) & (t["zHD"] < p.z_max)]
    rows = []
    for cid in np.unique(t["CID"]):
        s = t[t["CID"] == cid]
        w = 1 / s["m_b_corr_err_DIAG"] ** 2
        rows.append(
            (
                str(cid),
                float(s["RA"][0]),
                float(s["DEC"][0]),
                float(np.mean(s["zHD"])),
                float(np.sum(w * s["m_b_corr"]) / np.sum(w)),
                float(np.sum(w * s["c"]) / np.sum(w)),
                float(
                    np.sqrt(1 / np.sum(w)) * np.sqrt(len(s))
                ),  # duplicates share the SN: no √n gain
            )
        )
    return Table(rows=rows, names=["cid", "ra", "dec", "z", "m", "c", "err"])


def box_terms(ra: float, dec: float, half_arcsec: float) -> list[str]:
    h = half_arcsec / 3600
    dpart = f"t.dec BETWEEN {dec - h:.7f} AND {dec + h:.7f}"
    hr = h / np.cos(np.radians(abs(dec) + h))
    lo, hi = ra - hr, ra + hr
    if lo < 0:
        ranges = [(0.0, hi), (lo + 360, 360.0)]
    elif hi >= 360:
        ranges = [(lo, 360.0), (0.0, hi - 360)]
    else:
        ranges = [(lo, hi)]
    return [f"(t.ra BETWEEN {a:.7f} AND {b:.7f} AND {dpart})" for a, b in ranges]


def query(sne: Table, p: Params) -> Table:
    terms = [
        t
        for r, d in zip(sne["ra"], sne["dec"], strict=True)
        for t in box_terms(r, d, p.radius_arcsec)
    ]
    q = (
        f"SELECT {COLS} FROM ls_dr9.tractor t JOIN ls_dr9.photo_z p ON t.ls_id = p.ls_id "
        f"WHERE t.mag_z < {p.mag_z_max} AND t.type <> 'PSF' AND (" + " OR ".join(terms) + ")"
    )
    err = ""
    for attempt in range(4):
        try:
            r = requests.post(
                TAP,
                data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q},
                timeout=900,
            )
            if r.status_code == 200 and r.text.startswith("ra,"):
                return Table.read(io.BytesIO(r.content), format="ascii.csv", guess=False)
            err = r.text[:300]
        except requests.RequestException as e:
            err = str(e)
        time.sleep(10 * (attempt + 1))
    raise RuntimeError(f"Data Lab TAP failed: {err}")


def fetch(p: Params) -> None:
    sne = load_pantheon(p)
    sne = sne[np.argsort(sne["ra"])]
    chunks = [sne[i : i + p.batch] for i in range(0, len(sne), p.batch)]

    def one(k: int) -> str:
        f = cache() / f"gal_{k:03d}.ecsv"
        if f.exists():
            return f"{k} cached"
        t0 = time.time()
        g = query(chunks[k], p)
        tmp = f.with_suffix(".tmp")
        g.write(tmp, format="ascii.ecsv", overwrite=True)
        tmp.replace(f)
        return f"{k} {len(g)} rows {time.time() - t0:.0f}s"

    with ThreadPoolExecutor(p.workers) as ex:
        for msg in ex.map(one, range(len(chunks))):
            print(msg, flush=True)


def ang_sep_arcsec(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.radians, (ra1, dec1, ra2, dec2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return np.degrees(np.arccos(np.clip(c, -1, 1))) * 3600


def lens_weight(zg: np.ndarray, zs: float, cosmo: FlatLambdaCDM) -> np.ndarray:
    """Convergence kernel (1 + z_l) chi_l (chi_s - chi_l) / chi_s (flat, Gpc); 0 behind the SN."""
    cl = cosmo.comoving_distance(zg).value / 1e3
    cs = cosmo.comoving_distance(zs).value / 1e3
    return np.clip((1 + zg) * cl * (cs - cl) / cs, 0, None)


def columns(sne: Table, gal: Table, p: Params) -> dict:
    """Per-SN weighted foreground counts in a disc, relative to the mean-shell-density expectation.

    Mean shell densities are per photometric region (release 9010 north, 9011 south: different
    depths), estimated from all SN discs of that region. X = sum(w) / E[sum(w)] - 1 (derived).
    """
    edges = np.arange(p.z_gal_min, p.z_max + p.shell_dz, p.shell_dz)
    mid = 0.5 * (edges[1:] + edges[:-1])
    zg_all = np.where(gal["z_spec"] > 0, gal["z_spec"], gal["z_phot_median"])
    counts = np.zeros((len(sne), len(mid)))
    region = np.full(len(sne), -1)
    ngal = np.zeros(len(sne), int)
    # galaxies near each SN (both arrays modest: ~1e3 SNe, ~1e5 galaxies)
    order = np.argsort(gal["dec"])
    gdec = np.asarray(gal["dec"])[order]
    h = p.radius_arcsec / 3600
    for i, s in enumerate(sne):
        lo, hi = np.searchsorted(gdec, [s["dec"] - h, s["dec"] + h])
        idx = order[lo:hi]
        sep = ang_sep_arcsec(s["ra"], s["dec"], np.asarray(gal["ra"])[idx], gdec[lo:hi])
        idx = idx[sep < p.radius_arcsec]
        if len(idx) == 0:
            continue
        rel = np.asarray(gal["release"])[idx]
        region[i] = int(np.bincount(rel - 9010, minlength=2).argmax())
        zg = zg_all[idx]
        keep = (zg >= p.z_gal_min) & (zg < s["z"] - p.host_gap)
        counts[i] = np.histogram(zg[keep], edges)[0]
        ngal[i] = len(idx)
    xl, xf, nbar = xcols(counts, region, sne, p)
    return {
        "xl": xl,
        "xf": xf,
        "counts": counts,
        "region": region,
        "has": region >= 0,
        "ngal": ngal,
        "nbar": nbar,
        "mid": mid.tolist(),
    }


def xcols(counts, region, sne, p: Params, nbar: dict | None = None):
    """Columns from per-shell counts; mean shell densities from the data unless ``nbar``."""
    cosmo = FlatLambdaCDM(H0=70, Om0=p.om)
    edges = np.arange(p.z_gal_min, p.z_max + p.shell_dz, p.shell_dz)
    mid = 0.5 * (edges[1:] + edges[:-1])
    xl = np.full(len(sne), np.nan)
    xf = np.full(len(sne), np.nan)
    nbar = {} if nbar is None else nbar
    for reg in (0, 1):
        m = region == reg
        if m.sum() == 0:
            continue
        if reg in nbar:
            nb = np.asarray(nbar[reg])
        else:  # shell densities: only discs whose SN lies behind the shell contribute
            nb = np.zeros(len(mid))
            for k in range(len(mid)):
                behind = m & (sne["z"] - p.host_gap > edges[k + 1])
                nb[k] = counts[behind, k].mean() if behind.any() else 0.0
            nbar[reg] = nb.tolist()
        for i in np.flatnonzero(m):
            ok = edges[1:] <= sne["z"][i] - p.host_gap
            if not ok.any():
                continue
            wl = lens_weight(mid[ok], sne["z"][i], cosmo)
            el = np.sum(wl * nb[ok])
            ef = np.sum(nb[ok])
            if el > 0 and ef > 0:
                xl[i] = np.sum(wl * counts[i, ok]) / el - 1
                if p.alpha == 1:
                    xf[i] = np.sum(counts[i, ok]) / ef - 1
                else:  # A3: (1 + delta)^2 per shell, Poisson-unbiased as N(N - 1) / E^2
                    # shells expecting < min_shell_expect galaxies are shot noise only: left out
                    n, e = counts[i, ok], nb[ok]
                    use = e >= p.min_shell_expect
                    if use.any():
                        xf[i] = np.sum(n[use] * (n[use] - 1) / e[use]) / np.sum(e[use]) - 1
    return xl, xf, nbar


def attenuation(
    c: dict, sne: Table, p: Params, rng, n_sim: int = 20, xf_max: float = np.inf
) -> tuple[float, float]:
    """Errors-in-variables factor lambda = 1 - var(Poisson-only column) / var(observed column).

    Shot noise in the counts dilutes a true column-brightness relation by lambda (regression
    dilution); a coefficient fitted on the noisy column must be divided by lambda. Per-column
    approximation (ASSUMPTION): ignores the noise covariance between the two columns.
    """
    region = c["region"]
    expect = np.zeros_like(c["counts"])
    for reg, nb in c["nbar"].items():
        expect[region == reg] = np.asarray(nb)
    ok = np.isfinite(c["xl"]) & np.isfinite(c["xf"])
    ok[ok] &= c["xf"][ok] <= xf_max
    vl, vf = [], []
    for _ in range(n_sim):
        sl, sf, _ = xcols(rng.poisson(expect).astype(float), region, sne, p, c["nbar"])
        vl.append(np.nanvar(sl[ok]))
        vf.append(np.nanvar(sf[ok]))
    return (
        float(1 - np.mean(vl) / np.var(c["xl"][ok])),
        float(1 - np.mean(vf) / np.var(c["xf"][ok])),
    )


def residuals(sne: Table, p: Params) -> np.ndarray:
    """m_b_corr minus the flat-LCDM distance modulus (model_prediction); offset fitted later."""
    cosmo = FlatLambdaCDM(H0=70, Om0=p.om)
    return np.asarray(sne["m"]) - cosmo.distmod(np.asarray(sne["z"])).value


def wls(y, err, z, xl, xf):
    """WLS on [1, z, X_lens, X_flat]; coefficients and errors scaled by chi2/dof."""
    a = np.column_stack([np.ones_like(z), z, xl, xf])
    w = 1 / err**2
    cov = np.linalg.inv(a.T @ (a * w[:, None]))
    beta = cov @ (a.T @ (w * y))
    chi2 = np.sum(w * (y - a @ beta) ** 2)
    dof = len(y) - a.shape[1]
    return beta, np.sqrt(np.diag(cov) * chi2 / dof), chi2 / dof


def scramble(z, xl, xf, rng, dz=0.05):
    """Permute column pairs among SNe in one z_s bin (keeps z dependence, breaks the sky match)."""
    out_l, out_f = xl.copy(), xf.copy()
    bins = np.floor(z / dz).astype(int)
    for b in np.unique(bins):
        i = np.flatnonzero(bins == b)
        j = rng.permutation(i)
        out_l[i], out_f[i] = xl[j], xf[j]
    return out_l, out_f


def fit(p: Params) -> dict:
    sne = load_pantheon(p)
    files = sorted(cache().glob("gal_*.ecsv"))
    if not files:
        raise SystemExit("run `fetch` first")
    gal = vstack([Table.read(f) for f in files])
    gal = gal[
        np.unique(
            np.column_stack([np.round(gal["ra"], 6), np.round(gal["dec"], 6)]),
            axis=0,
            return_index=True,
        )[1]
    ]
    c = columns(sne, gal, p)
    y = residuals(sne, p)
    ok = np.isfinite(c["xl"]) & np.isfinite(c["xf"])
    cols = (sne["z"], sne["err"], c["xl"], c["xf"], y, sne["c"])
    z, err, xl, xf, y, col = (np.asarray(v)[ok] for v in cols)
    # robust: drop > 5 sigma Hubble-residual outliers of the base fit (ASSUMPTION; count reported)
    a0 = np.column_stack([np.ones_like(z), z])
    b0 = np.linalg.lstsq(a0 / err[:, None], y / err, rcond=None)[0]
    clip = np.abs(y - a0 @ b0) < 5 * err
    z, err, xl, xf, y, col = z[clip], err[clip], xl[clip], xf[clip], y[clip], col[clip]
    beta, sig, rchi2 = wls(y, err, z, xl, xf)
    rng = np.random.default_rng(p.seed)
    null = np.array([wls(y, err, z, *scramble(z, xl, xf, rng))[0][2:] for _ in range(p.n_scramble)])
    gf_inj = 3 * sig[3]
    rec = []
    for _ in range(p.n_inject):
        # carrier: a scrambled column pair, so the real sky signal in y is uncorrelated with it
        sl, sf = scramble(z, xl, xf, rng)
        bi = wls(y + gf_inj * sf, err, z, sl, sf)[0]
        rec.append(bi[3])
    rec = np.array(rec)
    corr = float(np.corrcoef(xl, xf)[0, 1])
    bc, ec, _ = wls(col, np.full_like(col, 0.05), z, xl, xf)  # errors rescaled by chi2/dof
    res = {
        "params": asdict(p),
        "provenance": {
            "pantheon_plus": "observed (m_b_corr, Brout et al. 2022 / Scolnic et al. 2022)",
            "galaxies": "observed (LS DR9 Tractor + DR9 photo-z, Data Lab TAP)",
            "columns": "derived",
            "residual_model": "model_prediction (flat LCDM Om = 0.334)",
            "thresholds": "assumption (Params)",
        },
        "n_sne_input": len(sne),
        "n_sne_with_columns": int(ok.sum()),
        "n_clipped_5sigma": int((~clip).sum()),
        "n_sne_fit": len(y),
        "n_galaxies": len(gal),
        "corr_xl_xf": corr,
        "std_xl": float(np.std(xl)),
        "std_xf": float(np.std(xf)),
        "reduced_chi2": float(rchi2),
        "gamma_lens": float(beta[2]),
        "gamma_lens_err": float(sig[2]),
        "gamma_flat": float(beta[3]),
        "gamma_flat_err": float(sig[3]),
        # grey-dust mimic (B-on-A3 P2b item 2): SALT colour c regressed on the same columns
        "colour_lens": float(bc[2]),
        "colour_lens_err": float(ec[2]),
        "colour_flat": float(bc[3]),
        "colour_flat_err": float(ec[3]),
        "scramble_std_lens": float(null[:, 0].std()),
        "scramble_std_flat": float(null[:, 1].std()),
        "p_scramble_flat": float(np.mean(np.abs(null[:, 1]) >= abs(beta[3]))),
        "p_scramble_lens": float(np.mean(np.abs(null[:, 0]) >= abs(beta[2]))),
        "inject_gamma_flat": float(gf_inj),
        "inject_recovered_mean": float(rec.mean()),
        "inject_recovered_std": float(rec.std()),
        "inject_fraction_detected_2sigma": float(np.mean(rec / sig[3] > 2)),
    }
    sf_cal = max(sig[3], res["scramble_std_flat"])
    res["gamma_flat_95"] = float(abs(beta[3]) + 1.96 * sf_cal)
    res["per_unit"] = "mag per unit fractional excess of the weighted foreground count (X)"
    # robustness: drop the sightlines above the p.trim_pct percentile of X_flat (leverage)
    keep = xf <= np.percentile(xf, p.trim_pct)
    bt, st, _ = wls(y[keep], err[keep], z[keep], xl[keep], xf[keep])
    nt = np.array(
        [
            wls(y[keep], err[keep], z[keep], *scramble(z[keep], xl[keep], xf[keep], rng))[0][3]
            for _ in range(p.n_scramble)
        ]
    )
    res["trim_n"] = int(keep.sum())
    res["gamma_flat_trim"] = float(bt[3])
    res["gamma_flat_trim_err"] = float(max(st[3], nt.std()))
    # shot noise dilutes the coefficient (regression dilution): divide by lambda
    lam_l, lam_f = attenuation(c, sne, p, rng)
    lam_ft = attenuation(c, sne, p, rng, xf_max=float(np.percentile(xf, p.trim_pct)))[1]
    res["lambda_lens"], res["lambda_flat"], res["lambda_flat_trim"] = lam_l, lam_f, lam_ft
    # A3's sign is +gamma (fainter); its unit is the column normalised to its mean, T / <T>
    norm = float(np.mean(1 + xf))
    up_full = max(beta[3] + 1.645 * sf_cal, 0) / max(lam_f, 1e-3)
    up_trim = max(bt[3] + 1.645 * res["gamma_flat_trim_err"], 0) / max(lam_ft, 1e-3)
    res["norm_mean_1_plus_x"] = norm
    res["gamma_norm_upper95_one_sided"] = float(max(up_full, up_trim) * norm)
    res["gamma_norm_upper95_unit"] = "mag per unit T/<T>, corrected for shot-noise dilution"
    return res


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("cmd", choices=["fetch", "fit"])
    ap.add_argument("--z-min", type=float, default=None, help="fit only SNe above this z")
    ap.add_argument("--alpha", type=int, default=1, choices=[1, 2], help="flat column power")
    a = ap.parse_args()
    p = Params()
    fit_p = Params(alpha=a.alpha, **({} if a.z_min is None else {"z_min": a.z_min}))
    if a.cmd == "fetch":
        fetch(p)
    else:
        res = fit(fit_p)
        OUT.mkdir(parents=True, exist_ok=True)
        tag = ("" if a.z_min is None else f"_zmin{a.z_min:g}") + f"_alpha{a.alpha}"
        (OUT / f"fit{tag}.json").write_text(json.dumps(res, indent=1) + "\n")
        print(json.dumps({k: v for k, v in res.items() if k != "params"}, indent=1))


if __name__ == "__main__":
    main()
