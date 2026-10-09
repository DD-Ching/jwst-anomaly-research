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
import hashlib
import io
import json
import time
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np
import requests
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table, vstack
from w12_lenscats import box_terms

from jwst_anomaly import paths

TAP = "https://datalab.noirlab.edu/tap/sync"
PPLUS_URL = (
    "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
    "Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat"
)
COLS = "t.ra, t.dec, t.dered_mag_z, t.release, p.z_phot_median, p.z_spec"
QUERY_VERSION = 2  # 2: brick_primary, maskbits = 0, no DUP, dereddened z (countmap selection)
FETCH_Z = (0.1, 1.3)  # SN redshift range of the galaxy fetch
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
    batch: int = 60  # boxes per TAP query
    n_scramble: int = 1000
    n_chain: int = 20  # whole-chain injections (thinned counts)
    chain_gamma: float = 0.01  # mag per unit X injected in the chain test (A3's fiducial scale)
    min_chain_ratio: float = 0.3  # ASSUMPTION: below this recovery the test sets no limit
    min_coverage: float = 0.5  # ASSUMPTION: disc galaxy count >= half the region median (edges)
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


def galaxy_dir(p: Params) -> Path:
    """Galaxy chunks keyed by everything that selects them (never reuse stale chunks).

    The SN list of the fetch is always the default z range; fits on a narrower range use a subset.
    """
    key = json.dumps([QUERY_VERSION, FETCH_Z, p.radius_arcsec, p.mag_z_max, p.batch])
    d = cache() / f"gal_{hashlib.sha1(key.encode()).hexdigest()[:10]}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def query(sne: Table, p: Params) -> Table:
    terms = [
        t
        for r, d in zip(sne["ra"], sne["dec"], strict=True)
        for t in box_terms(r, d, p.radius_arcsec)  # unqualified ra/dec: only the tractor has them
    ]
    q = (
        f"SELECT {COLS} FROM ls_dr9.tractor t JOIN ls_dr9.photo_z p ON t.ls_id = p.ls_id "
        f"WHERE t.brick_primary = 1 AND t.maskbits = 0 AND t.type <> 'PSF' AND t.type <> 'DUP' "
        f"AND t.dered_mag_z > 0 AND t.dered_mag_z < {p.mag_z_max} AND (" + " OR ".join(terms) + ")"
    )
    err = ""
    for attempt in range(4):
        if attempt:
            time.sleep(10 * attempt)
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
    raise RuntimeError(f"Data Lab TAP failed: {err}")


def fetch(p: Params) -> None:
    """Sequential batched queries (Data Lab etiquette: batch, no concurrency)."""
    sne = load_pantheon(replace(p, z_min=FETCH_Z[0], z_max=FETCH_Z[1]))
    sne = sne[np.argsort(sne["ra"])]
    out = galaxy_dir(p)
    for k in range(0, len(sne), p.batch):
        f = out / f"gal_{k // p.batch:03d}.ecsv"
        if f.exists():
            continue
        t0 = time.time()
        g = query(sne[k : k + p.batch], p)
        tmp = f.with_suffix(".tmp")
        g.write(tmp, format="ascii.ecsv", overwrite=True)
        tmp.replace(f)
        print(f"{k // p.batch} {len(g)} rows {time.time() - t0:.0f}s", flush=True)


def ang_sep_arcsec(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.radians, (ra1, dec1, ra2, dec2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return np.degrees(np.arccos(np.clip(c, -1, 1))) * 3600


def lens_weight(zg: np.ndarray, zs: float, cosmo: FlatLambdaCDM) -> np.ndarray:
    """Convergence kernel (1 + z_l) chi_l (chi_s - chi_l) / chi_s (flat, Gpc); 0 behind the SN."""
    cl = cosmo.comoving_distance(zg).value / 1e3
    cs = cosmo.comoving_distance(zs).value / 1e3
    return np.clip((1 + zg) * cl * (cs - cl) / cs, 0, None)


REGION = {9010: 0, 9011: 1}  # DR9 south (DECam) and north (BASS/MzLS): different depths


def shell_edges(p: Params) -> tuple[np.ndarray, np.ndarray]:
    edges = np.arange(p.z_gal_min, p.z_max + p.shell_dz, p.shell_dz)
    return edges, 0.5 * (edges[1:] + edges[:-1])


def geometry(sne: Table, p: Params) -> list:
    """Per SN: foreground-shell mask and lensing weights (computed once, reused by every xcols)."""
    cosmo = FlatLambdaCDM(H0=70, Om0=p.om)
    edges, mid = shell_edges(p)
    out = []
    for zs in np.asarray(sne["z"]):
        ok = edges[1:] <= zs - p.host_gap
        out.append((ok, lens_weight(mid[ok], zs, cosmo) if ok.any() else None))
    return out


def columns(sne: Table, gal: Table, p: Params) -> dict:
    """Per-SN weighted foreground counts in a disc, relative to the mean-shell-density expectation.

    Mean shell densities are per photometric region (``REGION``), estimated from all SN discs of
    that region. X = sum(w) / E[sum(w)] - 1 (derived). ``ngal`` (all redshifts) is the coverage
    proxy for discs cut by the footprint edge or masks.
    """
    edges, mid = shell_edges(p)
    zg_all = np.where(gal["z_spec"] > 0, gal["z_spec"], gal["z_phot_median"])
    counts = np.zeros((len(sne), len(mid)))
    region = np.full(len(sne), -1)
    ngal = np.zeros(len(sne), int)
    order = np.argsort(gal["dec"])
    gdec = np.asarray(gal["dec"])[order]
    gra = np.asarray(gal["ra"])[order]
    grel = np.asarray(gal["release"])[order]
    gz = np.asarray(zg_all)[order]
    h = p.radius_arcsec / 3600
    for i, s in enumerate(sne):
        lo, hi = np.searchsorted(gdec, [s["dec"] - h, s["dec"] + h])
        sel = lo + np.flatnonzero(
            ang_sep_arcsec(s["ra"], s["dec"], gra[lo:hi], gdec[lo:hi]) < p.radius_arcsec
        )
        regs = [REGION[r] for r in grel[sel] if r in REGION]
        if not regs:
            continue
        region[i] = int(np.bincount(regs, minlength=2).argmax())
        zg = gz[sel]
        keep = (zg >= p.z_gal_min) & (zg < s["z"] - p.host_gap)
        counts[i] = np.histogram(zg[keep], edges)[0]
        ngal[i] = len(sel)
    geo = geometry(sne, p)
    xl, xf, nbar = xcols(counts, region, sne, p, geo)
    return {
        "xl": xl,
        "xf": xf,
        "counts": counts,
        "region": region,
        "ngal": ngal,
        "nbar": nbar,
        "geo": geo,
        "mid": mid.tolist(),
    }


def xcols(counts, region, sne, p: Params, geo: list, nbar: dict | None = None):
    """Columns from per-shell counts; mean shell densities from the data unless ``nbar``."""
    edges, mid = shell_edges(p)
    xl = np.full(len(sne), np.nan)
    xf = np.full(len(sne), np.nan)
    nbar = {} if nbar is None else nbar
    zs = np.asarray(sne["z"])
    for reg in (0, 1):
        m = region == reg
        if m.sum() == 0:
            continue
        if reg in nbar:
            nb = np.asarray(nbar[reg])
        else:  # shell densities: only discs whose SN lies behind the shell contribute
            nb = np.zeros(len(mid))
            for k in range(len(mid)):
                behind = m & (zs - p.host_gap > edges[k + 1])
                nb[k] = counts[behind, k].mean() if behind.any() else 0.0
            nbar[reg] = nb.tolist()
        for i in np.flatnonzero(m):
            ok, wl = geo[i]
            if wl is None:
                continue
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


def _trimmed_var(x: np.ndarray, trim_pct: float | None) -> float:
    x = x[np.isfinite(x)]
    if trim_pct is not None:
        x = x[x <= np.percentile(x, trim_pct)]
    return float(np.var(x))


def attenuation(
    c: dict, sne: Table, p: Params, rng, sel: np.ndarray, trim_pct: float | None = None, n_sim=20
) -> tuple[float, float]:
    """Errors-in-variables factor lambda = 1 - var(Poisson-only column) / var(observed column).

    Shot noise in the counts dilutes a true column-brightness relation by lambda (regression
    dilution); a coefficient fitted on the noisy column must be divided by lambda. Observed and
    simulated columns get the same selection ``sel`` and the same trimming (each at its own
    percentile). Per-column approximation (ASSUMPTION): ignores the noise covariance between the
    two columns and clustering beyond Poisson in the noise model.
    """
    region = c["region"]
    expect = np.zeros_like(c["counts"])
    for reg, nb in c["nbar"].items():
        expect[region == reg] = np.asarray(nb)
    vl, vf = [], []
    for _ in range(n_sim):
        sim = rng.poisson(expect).astype(float)
        sl, sf, _ = xcols(sim, region, sne, p, c["geo"], c["nbar"])
        vl.append(_trimmed_var(sl[sel], trim_pct))
        vf.append(_trimmed_var(sf[sel], trim_pct))
    return (
        float(1 - np.mean(vl) / _trimmed_var(c["xl"][sel], trim_pct)),
        float(1 - np.mean(vf) / _trimmed_var(c["xf"][sel], trim_pct)),
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


def z_perm(z: np.ndarray, rng, dz=0.05) -> np.ndarray:
    """Index permutation within z_s bins."""
    out = np.arange(len(z))
    bins = np.floor(z / dz).astype(int)
    for b in np.unique(bins):
        i = np.flatnonzero(bins == b)
        out[i] = rng.permutation(i)
    return out


def selection(c: dict, y: np.ndarray, err: np.ndarray, z: np.ndarray, p: Params) -> np.ndarray:
    """Fit sample: finite columns, disc coverage, and no > 5 sigma Hubble outliers (ASSUMPTIONs)."""
    ok = np.isfinite(c["xl"]) & np.isfinite(c["xf"])
    for reg in (0, 1):
        m = ok & (c["region"] == reg)
        if m.any():
            ok[m] &= c["ngal"][m] >= p.min_coverage * np.median(c["ngal"][m])
    a0 = np.column_stack([np.ones_like(z), z])
    b0 = np.linalg.lstsq(a0[ok] / err[ok, None], y[ok] / err[ok], rcond=None)[0]
    ok &= np.abs(y - a0 @ b0) < 5 * err
    return ok


def limit(c, sne, y, err, z, sel, p: Params, rng) -> dict:
    """Coefficients, scramble errors, trimmed refit, dilution and the one-sided 95 % limit."""
    xl, xf = c["xl"][sel], c["xf"][sel]
    ys, es, zs = y[sel], err[sel], z[sel]
    beta, sig, rchi2 = wls(ys, es, zs, xl, xf)
    null = np.array(
        [wls(ys, es, zs, *scramble(zs, xl, xf, rng))[0][2:] for _ in range(p.n_scramble)]
    )
    keep = xf <= np.percentile(xf, p.trim_pct)
    bt, st, _ = wls(ys[keep], es[keep], zs[keep], xl[keep], xf[keep])
    nt = np.array(
        [
            wls(ys[keep], es[keep], zs[keep], *scramble(zs[keep], xl[keep], xf[keep], rng))[0][3]
            for _ in range(p.n_scramble // 4)
        ]
    )
    lam_l, lam_f = attenuation(c, sne, p, rng, sel)
    lam_ft = attenuation(c, sne, p, rng, sel, p.trim_pct)[1]
    sf_cal = max(sig[3], null[:, 1].std())
    st_cal = max(st[3], nt.std())
    norm = float(np.mean(1 + xf))
    up_full = max(beta[3] + 1.645 * sf_cal, 0) / max(lam_f, 1e-3)
    up_trim = max(bt[3] + 1.645 * st_cal, 0) / max(lam_ft, 1e-3)
    return {
        "beta": beta,
        "sig": sig,
        "rchi2": rchi2,
        "null": null,
        "trim_n": int(keep.sum()),
        "bt": bt,
        "st_cal": st_cal,
        "sf_cal": sf_cal,
        "lam": (lam_l, lam_f, lam_ft),
        "norm": norm,
        "upper": float(max(up_full, up_trim) * norm),
    }


def chain_injection(c, sne, y, err, z, p: Params, rng, gamma: float) -> np.ndarray:
    """Whole-chain check of the dilution correction (scripts/CLAUDE.md: inject through the chain).

    Truth: the full-count column of a z-matched permutation of sightlines (so the real sky signal
    in y is uncorrelated with it), y += gamma * X_true. Observed: the same counts binomially thinned
    to half, then the whole chain (shell densities, columns, selection, fit, trimming, dilution).
    Returns the recovered gamma estimates (fit / lambda, full and trimmed).
    """
    pt = replace(p, min_shell_expect=p.min_shell_expect / 2)  # same shells after thinning
    out = []
    for _ in range(p.n_chain):
        # SN i gets the whole sightline (counts, shells, weights) of a donor SN at similar z
        perm = z_perm(z, rng)
        counts, region, sne_p = c["counts"][perm], c["region"][perm], sne[perm]
        geo = [c["geo"][j] for j in perm]
        x_true = xcols(counts, region, sne_p, p, geo)[1]
        thin = rng.binomial(counts.astype(int), 0.5).astype(float)
        xl_t, xf_t, nb_t = xcols(thin, region, sne_p, pt, geo)
        ct = {**c, "counts": thin, "region": region, "xl": xl_t, "xf": xf_t, "nbar": nb_t}
        ct["ngal"], ct["geo"] = c["ngal"][perm], geo
        yi = y + gamma * np.nan_to_num(x_true)
        sel = selection(ct, yi, err, z, pt) & np.isfinite(x_true)
        xs = xf_t[sel]
        keep = xs <= np.percentile(xs, p.trim_pct)
        b = wls(yi[sel], err[sel], z[sel], xl_t[sel], xs)[0][3]
        bt = wls(yi[sel][keep], err[sel][keep], z[sel][keep], xl_t[sel][keep], xs[keep])[0][3]
        lam_f = attenuation(ct, sne_p, pt, rng, sel, n_sim=5)[1]
        lam_ft = attenuation(ct, sne_p, pt, rng, sel, p.trim_pct, n_sim=5)[1]
        out.append((b / lam_f, bt / lam_ft))
    return np.array(out)


def fit(p: Params) -> dict:
    sne = load_pantheon(p)
    files = sorted(galaxy_dir(p).glob("gal_*.ecsv"))
    if not files:
        raise SystemExit("run `fetch` first")
    gal = vstack([Table.read(f) for f in files])
    # a galaxy inside two SN boxes of different chunks comes back twice with identical values
    gal = gal[np.unique(np.column_stack([gal["ra"], gal["dec"]]), axis=0, return_index=True)[1]]
    c = columns(sne, gal, p)
    y = residuals(sne, p)
    err, z, colour = (np.asarray(sne[k]) for k in ("err", "z", "c"))
    sel = selection(c, y, err, z, p)
    rng = np.random.default_rng(p.seed)
    lim = limit(c, sne, y, err, z, sel, p, rng)
    beta, sig = lim["beta"], lim["sig"]
    inj = chain_injection(c, sne, y, err, z, p, rng, p.chain_gamma)
    ratio = float(min(inj[:, 0].mean(), inj[:, 1].mean()) / p.chain_gamma)
    xl, xf = c["xl"][sel], c["xf"][sel]
    bc, ec, _ = wls(colour[sel], np.full(sel.sum(), 0.05), z[sel], xl, xf)  # rescaled errors
    finite = np.isfinite(c["xl"]) & np.isfinite(c["xf"])
    return {
        "params": asdict(p),
        "provenance": {
            "pantheon_plus": "observed (m_b_corr, Brout et al. 2022 / Scolnic et al. 2022)",
            "pantheon_plus_sha256": hashlib.sha256(
                (cache() / "PantheonPlusSH0ES.dat").read_bytes()
            ).hexdigest(),
            "galaxies": "observed (LS DR9 Tractor + DR9 photo-z, Data Lab TAP)",
            "columns": "derived",
            "residual_model": "model_prediction (flat LCDM Om = 0.334)",
            "thresholds": "assumption (Params)",
        },
        "n_sne_input": len(sne),
        "n_sne_with_columns": int(finite.sum()),
        "n_cut_coverage_or_5sigma": int((finite & ~sel).sum()),
        "n_sne_fit": int(sel.sum()),
        "n_galaxies": len(gal),
        "corr_xl_xf": float(np.corrcoef(xl, xf)[0, 1]),
        "std_xl": float(np.std(xl)),
        "std_xf": float(np.std(xf)),
        "reduced_chi2": float(lim["rchi2"]),
        "gamma_lens": float(beta[2]),
        "gamma_lens_err": float(max(sig[2], lim["null"][:, 0].std())),
        "gamma_flat": float(beta[3]),
        "gamma_flat_err": float(lim["sf_cal"]),
        "p_scramble_flat": float(np.mean(np.abs(lim["null"][:, 1]) >= abs(beta[3]))),
        "p_scramble_lens": float(np.mean(np.abs(lim["null"][:, 0]) >= abs(beta[2]))),
        # grey-dust mimic (B-on-A3 P2b item 2): SALT colour c regressed on the same columns
        "colour_lens": float(bc[2]),
        "colour_lens_err": float(ec[2]),
        "colour_flat": float(bc[3]),
        "colour_flat_err": float(ec[3]),
        "trim_n": lim["trim_n"],
        "gamma_flat_trim": float(lim["bt"][3]),
        "gamma_flat_trim_err": float(lim["st_cal"]),
        "lambda_lens": lim["lam"][0],
        "lambda_flat": lim["lam"][1],
        "lambda_flat_trim": lim["lam"][2],
        "chain_inject_gamma": p.chain_gamma,
        "chain_recovered_full": [float(inj[:, 0].mean()), float(inj[:, 0].std())],
        "chain_recovered_trim": [float(inj[:, 1].mean()), float(inj[:, 1].std())],
        "per_unit": "mag per unit fractional excess of the weighted foreground count (X)",
        "norm_mean_1_plus_x": lim["norm"],
        # A3's sign is +gamma (fainter); its unit is the column normalised to its mean, T / <T>
        # divided by the whole-chain recovery ratio when the dilution correction under-recovers
        "chain_recovery_ratio": ratio,
        # no limit when the chain does not recover the signal (ASSUMPTION: ratio < min_chain_ratio)
        "gamma_norm_upper95_one_sided": (
            lim["upper"] / min(ratio, 1.0) if ratio >= p.min_chain_ratio else None
        ),
        "gamma_norm_upper95_unit": "mag per unit T/<T>, corrected for shot-noise dilution",
    }


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
