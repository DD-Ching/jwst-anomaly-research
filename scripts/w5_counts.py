"""W5 count-deficit screen (D-063) on one Legacy Surveys DR10 region.

    python scripts/w5_counts.py fetch  --ra 150.2 153.8 --dec -0.05 3.55
    python scripts/w5_counts.py screen --n-inject 200

`fetch` stores the Tractor extract under $JWST_ANOMALY_DATA/w5/ (gitignored) and writes a manifest.
`screen` bins it, runs the disk-deficit statistic for each Einstein radius, vets flags with the
star counts in the same disk (stars are foreground of any extragalactic lens, masks remove both),
injects predicted holes for efficiency and writes the result table and a contact sheet.
Thresholds are ASSUMPTIONs and are recorded in the output meta.
"""

from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np
from astropy.table import Table

from jwst_anomaly import countmap, paths, schema

CELL = 0.5  # arcmin
MAG_GAL = 21.0  # galaxy flux limit (r); DR10 is complete well below it (ASSUMPTION)
THETA_E = (4.0, 8.0, 11.0, 16.0, 22.0, 32.0)  # arcmin
DISK_FRAC = 0.3  # disk radius / theta_E: the predicted ratio is < 0.05 inside x = 0.3, 0.6 at 0.5
MIN_COV = 0.8  # minimum mean unmasked fraction in the disk
FAP = 0.05  # field-wide false-alarm probability under the clustered (negative-binomial) null
# O/E: deeper than the field's deepest clustering underdensities (O/E 0.24 at theta_E = 8' in the
# pilot with a 0.5 cut, all with normal star counts) and shallower than the predicted hole's
# disk-averaged ratio (``core_ratio`` in the output, ~0.03-0.1)
FLAG_RATIO = 0.15


def _dir():
    d = paths.data_root() / "w5"
    d.mkdir(parents=True, exist_ok=True)
    return d


def fetch(args):
    tab = countmap.fetch_tractor(tuple(args.ra), tuple(args.dec), args.mag_max)
    out = _dir() / "tractor.fits"
    tab.write(out, overwrite=True)
    manifest = {
        "uri": tab.meta["source"],
        "file": out.name,
        "sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
        "size": out.stat().st_size,
        "rows": len(tab),
    }
    man = paths.manifests_dir() / "w5_tractor_pilot.json"
    man.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


def _flags(gal, clean, every, r_cells, k=None):
    """Flags: O/E < FLAG_RATIO and n_indep * P_NB(N <= O) < FAP. Returns (O, E, log10 P_NB, hit, k).

    k (clustering) is fitted on the field without injections and reused for injected maps.
    """
    obs, exp, _ = countmap.disk_deficit(gal, clean, every, r_cells, MIN_COV)
    k_fit, n_ind = countmap.clustering_k(obs, exp, r_cells)
    k = k_fit if k is None else k
    lt = countmap.nb_log10_tail(obs, exp, k)
    with np.errstate(invalid="ignore", divide="ignore"):
        hit = (lt < np.log10(FAP / n_ind)) & (obs < FLAG_RATIO * exp)
    return obs, exp, lt, hit, k


def _peaks(hit, logp, sep):
    """Local minima of logp among flagged cells, at least ``sep`` cells apart (greedy)."""
    idx = np.argwhere(hit)
    idx = idx[np.argsort(logp[hit])]
    kept = []
    for p in idx:
        if all(np.hypot(*(p - q)) >= sep for q in kept):
            kept.append(p)
    return kept


def screen(args):
    tab = Table.read(_dir() / "tractor.fits")
    ra0 = 0.5 * (tab["ra"].min() + tab["ra"].max())
    dec0 = 0.5 * (tab["dec"].min() + tab["dec"].max())
    ra_span = (tab["ra"].max() - tab["ra"].min()) * np.cos(np.radians(dec0))
    half = 0.5 * 60 * min(tab["dec"].max() - tab["dec"].min(), ra_span) - 1  # square inside the box
    gal, star, clean, every = countmap.count_grids(tab, ra0, dec0, half, CELL, MAG_GAL)
    f_mask = clean.sum() / every.sum()
    area = gal.size * f_mask * (CELL / 60) ** 2  # unmasked area of the grid
    n_bar = gal.sum() / (gal.size * f_mask)
    mags = np.arange(16.0, MAG_GAL + 0.01, 0.5)
    psf = np.char.strip(np.asarray(tab["type"]).astype(str)) == "PSF"
    clean_gal = (np.asarray(tab["maskbits"]) & countmap.MASK_BITS) == 0
    clean_gal &= ~psf
    cum = [max((clean_gal & (tab["mag_r"] < m)).sum(), 1) / area for m in mags]
    n_brighter = countmap.power_law_counts(mags, cum)
    # vetting maps: every unmasked row brighter than MAG_GAL whatever its type (a star/galaxy type
    # swap conserves it; a W5 hole removes only the galaxy share) and every row (depth proxy)
    bright_any = ((np.asarray(tab["maskbits"]) & countmap.MASK_BITS) == 0) & (
        tab["mag_r"] < MAG_GAL
    )
    any_grid = countmap.grid(tab, np.asarray(bright_any), ra0, dec0, half, CELL)
    rng = np.random.default_rng(args.seed)
    rows, flags = [], []
    yy, xx = np.indices(gal.shape)
    for te in THETA_E:
        te_c, r_c = te / CELL, DISK_FRAC * te / CELL
        obs, exp, logp, hit, k = _flags(gal, clean, every, r_c)
        sobs, sexp, _ = countmap.disk_deficit(star, clean, every, r_c, MIN_COV)
        aobs, aexp, _ = countmap.disk_deficit(any_grid, clean, every, r_c, MIN_COV)
        ratio_fn = lambda x: countmap.predicted_ratio(x, n_brighter, MAG_GAL)  # noqa: E731
        core = _disk_mean(ratio_fn, DISK_FRAC)
        for p in _peaks(hit, logp, 2 * r_c):
            i, j = p
            ra, dec = _cell_radec(i, j, gal.shape[0], half, ra0, dec0)
            flags.append(
                dict(
                    theta_e=te,
                    ra=ra,
                    dec=dec,
                    obs=obs[i, j],
                    exp=exp[i, j],
                    log10p_nb=logp[i, j],
                    star_obs=sobs[i, j],
                    star_exp=sexp[i, j],
                    any_obs=aobs[i, j],
                    any_exp=aexp[i, j],
                    verdict=_verdict(
                        obs[i, j], exp[i, j], sobs[i, j], sexp[i, j], aobs[i, j], aexp[i, j], core
                    ),
                )
            )
        valid = np.argwhere(np.isfinite(logp))
        rec = 0
        for idx in rng.choice(len(valid), args.n_inject, replace=False):
            c = valid[idx]
            g2 = countmap.inject_hole(gal, n_bar, c, te_c, ratio_fn, rng)
            h2 = _flags(g2, clean, every, r_c, k)[3]
            near = np.hypot(yy - c[0], xx - c[1]) <= r_c
            rec += bool((h2 & near).any())
        eff = rec / args.n_inject
        n_valid_area = np.isfinite(logp).sum() * (CELL / 60) ** 2
        rows.append(
            dict(
                theta_e=te,
                n_flags=sum(f["theta_e"] == te for f in flags),
                recovered=rec,
                n_inject=args.n_inject,
                efficiency=eff,
                area_deg2=n_valid_area,
                core_ratio=core,
                min_oe=float(np.nanmin(obs / exp)),
                k_nb=k,
                min_log10p_nb=float(np.nanmin(logp)),
                n95_per_deg2=3.0 / (eff * n_valid_area) if eff > 0 else np.inf,
            )
        )
        print(rows[-1])
    res = Table(rows=rows)
    res.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source="scripts/w5_counts.py screen on " + str(_dir() / "tractor.fits"),
        cell_arcmin=CELL,
        mag_gal=MAG_GAL,
        disk_frac=DISK_FRAC,
        min_cov=MIN_COV,
        fap=FAP,
        flag_ratio=FLAG_RATIO,
        seed=args.seed,
        ra0=ra0,
        dec0=dec0,
        n_gal=int(gal.sum()),
        area_deg2=float(area),
    )
    out = paths.repo_root() / "results" / "w5" / "w5_counts_pilot.ecsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    res.write(out, overwrite=True)
    ft = Table(rows=flags) if flags else Table(names=("theta_e", "ra", "dec"))
    ft.meta.update(provenance=schema.Provenance.DERIVED.value, source=res.meta["source"])
    ft.write(out.with_name("w5_counts_pilot_flags.ecsv"), overwrite=True)
    print(ft)
    _sheet(gal, clean, every, star, flags, half, ra0, dec0, out.with_name("w5_counts_pilot.png"))


def _verdict(obs, exp, sobs, sexp, any_obs, any_exp, core):
    """Ordinary explanations a W5 hole cannot produce (thresholds are ASSUMPTIONs).

    type_swap: stars in excess (Poisson P(>= obs) < 0.05) where galaxies are missing; a lens behind
    the Milky Way's stars leaves their counts unchanged. not_hole: counts of every type brighter
    than MAG_GAL exceed what a hole leaves (non-galaxies untouched plus ``core`` x galaxies),
    P(>= any_obs) < 0.05. Otherwise open (needs /vet-candidate).
    """
    from scipy import stats

    reasons = []
    if stats.poisson.sf(sobs - 1, sexp) < 0.05:
        reasons.append("type_swap")
    if stats.poisson.sf(any_obs - 1, any_exp - exp + core * exp) < 0.05:
        reasons.append("not_hole")
    return ",".join(reasons) or "open"


def _disk_mean(ratio_fn, x_max, n=400):
    """Area-weighted mean predicted ratio inside x < x_max (model_prediction)."""
    x = (np.arange(n) + 0.5) * x_max / n
    with np.errstate(all="ignore"):
        r = np.nan_to_num(ratio_fn(x), nan=1.0)
    return float((r * x).sum() / x.sum())


def _cell_radec(i, j, n, half, ra0, dec0):
    x = -half + (j + 0.5) * CELL
    y = -half + (i + 0.5) * CELL
    # inverse gnomonic (small field)
    xr, yr = np.radians(x / 60), np.radians(y / 60)
    d0 = np.radians(dec0)
    rho = np.hypot(xr, yr)
    c = np.arctan(rho)
    dec = np.arcsin(np.cos(c) * np.sin(d0) + (yr * np.sin(c) * np.cos(d0) / rho if rho else 0))
    ra = ra0 + np.degrees(
        np.arctan2(xr * np.sin(c), rho * np.cos(d0) * np.cos(c) - yr * np.sin(d0) * np.sin(c))
    )
    return float(ra), float(np.degrees(dec))


def _sheet(gal, clean, every, star, flags, half, ra0, dec0, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    ext = (-half, half, -half, half)
    k = countmap.disk_kernel(4 / CELL)
    from scipy import signal

    obs, exp, _ = countmap.disk_deficit(gal, clean, every, 4 / CELL, 0.0)
    frac = signal.fftconvolve(clean, k, mode="same") / np.maximum(
        signal.fftconvolve(every, k, mode="same"), 1
    )
    for ax, img, title in (
        (axes[0], obs / exp, f"galaxies r<{MAG_GAL}: O/E in a 4' disk"),
        (axes[1], frac, "unmasked fraction of all rows, 4' disk"),
        (axes[2], signal.fftconvolve(star, k, mode="same"), "stars (PSF), 4' disk"),
    ):
        im = ax.imshow(img, origin="lower", extent=ext, cmap="viridis")
        fig.colorbar(im, ax=ax, shrink=0.8)
        ax.set_title(title)
        ax.set_xlabel("x [arcmin, +RA]")
    for f in flags:
        from jwst_anomaly.countmap import tangent_plane

        x, y = tangent_plane(np.array([f["ra"]]), np.array([f["dec"]]), ra0, dec0)
        for ax in axes:
            ax.add_patch(plt.Circle((x[0], y[0]), DISK_FRAC * f["theta_e"], fill=False, color="r"))
    fig.suptitle(f"W5 pilot, centre RA {ra0:.2f} Dec {dec0:.2f}; red: flags (disk = 0.5 theta_E)")
    fig.tight_layout()
    fig.savefig(path, dpi=80)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fetch")
    f.add_argument("--ra", type=float, nargs=2, default=(150.2, 153.8))
    f.add_argument("--dec", type=float, nargs=2, default=(-0.05, 3.55))
    f.add_argument("--mag-max", type=float, default=22.0)
    s = sub.add_parser("screen")
    s.add_argument("--n-inject", type=int, default=200)
    s.add_argument("--seed", type=int, default=63)
    args = ap.parse_args()
    {"fetch": fetch, "screen": screen}[args.cmd](args)


if __name__ == "__main__":
    main()
