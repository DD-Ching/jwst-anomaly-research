"""Multi-epoch, multi-band forced photometry of transient hits in every covering JWST mosaic.

A two-epoch difference hit (``arc_diff.py``) is one band at two dates. Vetting needs more: the
same visit's other bands (a single-band blob is an artifact), and further epochs (an SN rises
and fades once; an AGN varies repeatedly; a persistence or cosmic-ray block appears in one
visit only). This script finds every public level-3 NIRCam ``_i2d.fits`` that covers each
position (MAST, one query per position group and one product list), reads one small box per
mosaic and group (S3 byte ranges; nothing stored), and measures each position at its fixed sky
coordinate:

- ``flux_njy``: aperture sum (radius ``--radius-arcsec``) minus the median of an annulus
  (``ANNULUS``) times the aperture area;
- ``err_njy``: robust std (1.4826 MAD) of the same measurement on ``N_RANDOM`` random positions
  in the box, so source confusion and correlated noise are included (a conservative error);
- ``snr`` = flux / err.

A host galaxy adds a constant offset per band, so compare epochs **within** a band; PSF
differences between programs shift fluxes of extended light (stated per row by ``program``).
Every number is ``derived``; the radii and counts are ASSUMPTIONs. ``stamps_<name>.png`` shows
each mosaic's cutout (rows = bands, columns = dates) for looking at, which a light curve does
not replace.

    python scripts/hit_lightcurve.py --pos h003 3.577785 -30.394141 \\
        --pos h011 3.592116 -30.366745 --out results/transients/a2744_grid_alt/arcdiff_tiles/vet
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arc_diff import read_box  # noqa: E402
from transient_forced import robust_std  # noqa: E402

ANNULUS = (0.35, 0.6)  # ASSUMPTION: background annulus (arcsec), outside the PSF core
N_RANDOM = 300  # ASSUMPTION: random apertures per box for the empirical error
GROUP_ARCSEC = 20.0  # ASSUMPTION: positions closer than this share one box per mosaic
PAD_ARCSEC = 2.5  # ASSUMPTION: box margin beyond the group's positions
NJY_PER_MJYSR_SR = 1e15  # 1 MJy = 1e15 nJy
READ_RETRIES = 3  # transient S3 read failures are retried with back-off before a mosaic is skipped


def aperture_photometry(
    sci: np.ndarray, x: float, y: float, r_px: float, ann_px: tuple[float, float]
) -> float:
    """Aperture sum minus annulus median × aperture area (pixel units of ``sci``); NaN if off."""
    ny, nx = sci.shape
    r_out = ann_px[1]
    if not (r_out <= x < nx - r_out and r_out <= y < ny - r_out):
        return float("nan")
    # local cutout covering the annulus (not the whole box: this runs N_RANDOM times per mosaic)
    x0, y0 = max(int(np.floor(x - r_out)), 0), max(int(np.floor(y - r_out)), 0)
    x1, y1 = min(int(np.ceil(x + r_out)) + 1, nx), min(int(np.ceil(y + r_out)) + 1, ny)
    sub = sci[y0:y1, x0:x1]
    yy, xx = np.ogrid[y0:y1, x0:x1]
    d2 = (xx - x) ** 2 + (yy - y) ** 2
    ap = d2 <= r_px**2
    ann = (d2 >= ann_px[0] ** 2) & (d2 <= ann_px[1] ** 2)
    vals = sub[ap]
    bg = sub[ann]
    if not np.all(np.isfinite(vals)) or np.isfinite(bg).sum() < 10 or np.all(vals == 0):
        return float("nan")
    return float(vals.sum() - np.nanmedian(bg) * ap.sum())


def empirical_error(
    sci: np.ndarray, r_px: float, ann_px: tuple[float, float], rng: np.random.Generator
) -> float:
    """Robust std of :func:`aperture_photometry` at ``N_RANDOM`` random positions in ``sci``."""
    ny, nx = sci.shape
    m = ann_px[1] + 1
    if nx <= 2 * m or ny <= 2 * m:
        return float("nan")
    xs = rng.uniform(m, nx - m, N_RANDOM)
    ys = rng.uniform(m, ny - m, N_RANDOM)
    f = np.array(
        [aperture_photometry(sci, x, y, r_px, ann_px) for x, y in zip(xs, ys, strict=True)]
    )
    return robust_std(f, min_n=30)


def groups(positions: list[tuple[str, float, float]]) -> list[list[tuple[str, float, float]]]:
    """Greedy grouping of positions within ``GROUP_ARCSEC`` of a group's first member."""
    out: list[list[tuple[str, float, float]]] = []
    for p in positions:
        for g in out:
            dra = (p[1] - g[0][1]) * np.cos(np.radians(g[0][2])) * 3600
            if np.hypot(dra, (p[2] - g[0][2]) * 3600) < GROUP_ARCSEC:
                g.append(p)
                break
        else:
            out.append([p])
    return out


def covering_products(group: list[tuple[str, float, float]], only: set[str] | None = None) -> Table:
    """Public level-3 NIRCam image observations covering any group member, with cloud URIs.

    With ``only``, just those members are queried (controls near a hit share its mosaics).
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.table import unique, vstack

    from jwst_anomaly import query

    client = query.mast_client()
    obs = []
    for name, ra, dec in group:
        if only and name not in only:
            continue
        t = client.query_criteria(
            coordinates=SkyCoord(ra, dec, unit="deg"),
            radius=1 * u.arcsec,
            obs_collection="JWST",
            instrument_name="NIRCAM/IMAGE",
            dataproduct_type="image",
            calib_level=3,
            dataRights="PUBLIC",
        )
        if len(t):
            obs.append(t["obsid", "obs_id", "filters", "t_min", "proposal_id"])
    if not obs:
        return Table()
    obs = unique(vstack(obs), keys="obs_id")
    prods = query.list_products(obs, subgroups=("I2D",), cloud_uris=True)
    mjd = dict(zip(query.str_values(obs["obs_id"]), np.asarray(obs["t_min"], float), strict=True))
    rows = []
    for fname, uri in zip(
        query.str_values(prods["productFilename"]),
        query.str_values(prods["cloud_uri"]),
        strict=True,
    ):
        oid = fname.removesuffix("_i2d.fits")
        if uri and oid in mjd:
            band = oid.rsplit("_", 1)[-1].split("-")[-1].upper()
            rows.append((oid, band, mjd[oid], oid[2:7], uri))
    return Table(rows=rows, names=("obs_id", "band", "mjd", "program", "uri"))


def measure(job: dict) -> list[dict]:
    """One mosaic box: forced photometry of every member of the group, plus stamps."""
    g, prod, r_as = job["group"], job["prod"], job["radius"]
    ra0 = float(np.mean([p[1] for p in g]))
    dec0 = float(np.mean([p[2] for p in g]))
    ext = max(
        np.hypot((p[1] - ra0) * np.cos(np.radians(dec0)) * 3600, (p[2] - dec0) * 3600) for p in g
    )
    for attempt in range(READ_RETRIES):
        try:
            box = read_box(prod["uri"], ra0, dec0, ext + PAD_ARCSEC)
            break
        except Exception as exc:  # outside the mosaic (ValueError), or an unreadable file
            if isinstance(exc, ValueError) or attempt == READ_RETRIES - 1:
                return [
                    {"name": p[0], "obs_id": prod["obs_id"], "error": str(exc)[:120]} for p in g
                ]
            time.sleep(2**attempt)
    scale = box["scale"]
    to_njy = box["pixar_sr"] * NJY_PER_MJYSR_SR
    r_px = r_as / scale
    ann = (ANNULUS[0] / scale, ANNULUS[1] / scale)
    # stable seed (str hash() is salted per process, which made err_njy irreproducible)
    rng = np.random.default_rng(zlib.crc32(str(prod["obs_id"]).encode()))
    err = empirical_error(box["sci"], r_px, ann, rng) * to_njy
    out = []
    # registration correction (``--offsets``): a source at (ra, dec) of the reference appears at
    # (ra, dec) + (dx, dy) in this mosaic, so it is measured there
    dx_as, dy_as = job.get("offset", (0.0, 0.0))
    for name, ra, dec in g:
        ra_m = ra + dx_as / 3600.0 / np.cos(np.radians(dec))
        dec_m = dec + dy_as / 3600.0
        x, y = (float(v) for v in box["wcs"].all_world2pix(ra_m, dec_m, 0))
        f = aperture_photometry(box["sci"], x, y, r_px, ann) * to_njy
        h = int(np.ceil(1.2 / scale))
        xi, yi = int(round(x)), int(round(y))
        stamp = box["sci"][max(yi - h, 0) : yi + h + 1, max(xi - h, 0) : xi + h + 1] * to_njy
        out.append(
            {
                "name": name,
                "obs_id": prod["obs_id"],
                "band": prod["band"],
                "mjd": float(prod["mjd"]),
                "program": prod["program"],
                "flux_njy": f,
                "err_njy": err,
                "snr": f / err if err > 0 else float("nan"),
                "stamp": stamp,
            }
        )
    return out


def stamp_sheet(rows: list[dict], name: str, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = [r for r in rows if r["name"] == name and "stamp" in r and np.isfinite(r["flux_njy"])]
    if not rows:
        return
    bands = sorted(
        {r["band"] for r in rows}, key=lambda b: int("".join(c for c in b if c.isdigit()))
    )

    # one column per visit (obs_id without its filter), so two visits on one date don't overwrite
    def visit(r: dict) -> str:
        return str(r["obs_id"]).rsplit("_", 1)[0]

    first = {}
    for r in rows:
        first[visit(r)] = min(first.get(visit(r), np.inf), r["mjd"])
    visits = sorted(first, key=lambda v: (first[v], v))
    fig, axes = plt.subplots(
        len(bands),
        len(visits),
        figsize=(1.1 * len(visits) + 1, 1.1 * len(bands) + 0.6),
        squeeze=False,
    )
    for ax in axes.flat:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis("off")
    for r in rows:
        i, j = bands.index(r["band"]), visits.index(visit(r))
        ax = axes[i, j]
        s = r["stamp"]
        lo, hi = np.nanpercentile(s, [5, 99.5])
        ax.imshow(s, origin="lower", cmap="gray", vmin=lo, vmax=hi)
        ax.axis("on")
        ax.set_xticks([])
        ax.set_yticks([])
        c = (s.shape[1] - 1) / 2, (s.shape[0] - 1) / 2
        ax.plot(*c, "+", color="red", ms=5, mew=0.6)
        ax.set_title(f"{r['program']} {r['snr']:+.1f}", fontsize=5)
        if j == 0 or not axes[i, j - 1].axison:
            ax.set_ylabel(r["band"], fontsize=6)
    for j, v in enumerate(visits):
        axes[-1, j].set_xlabel(str(round(first[v])), fontsize=5)
    fig.suptitle(f'{name}: 2.4" stamps, title = program, forced S/N', fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=90)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pos", nargs=3, action="append", required=True, metavar=("NAME", "RA", "DEC"))
    ap.add_argument("--radius-arcsec", type=float, default=0.15)
    ap.add_argument(
        "--query-only",
        action="append",
        metavar="NAME",
        help="query coverage only at these positions (others, e.g. controls, share the mosaics)",
    )
    ap.add_argument(
        "--offsets",
        type=Path,
        help="JSON {obs_id: [dRA cos(dec), dDec] arcsec}: measured mosaic-minus-reference "
        "registration shifts (e.g. from reg_offsets.py); listed mosaics are measured there",
    )
    ap.add_argument(
        "--only-program",
        action="append",
        help="measure only mosaics of these programs (e.g. 01176)",
    )
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    offsets = json.loads(args.offsets.read_text()) if args.offsets else {}
    positions = [(n, float(r), float(d)) for n, r, d in args.pos]
    args.out.mkdir(parents=True, exist_ok=True)
    jobs = []
    for g in groups(positions):
        prods = covering_products(g, set(args.query_only or []))
        if args.only_program and len(prods):
            prods = prods[np.isin(prods["program"], args.only_program)]
        print(f"group {[p[0] for p in g]}: {len(prods)} mosaics", flush=True)
        jobs += [
            {
                "group": g,
                "prod": dict(zip(prods.colnames, r, strict=True)),
                "radius": args.radius_arcsec,
                "offset": tuple(offsets.get(str(r["obs_id"]), (0.0, 0.0))),
            }
            for r in prods
        ]
    with ThreadPoolExecutor(args.threads) as pool:
        rows = [r for res in pool.map(measure, jobs) for r in res]
    good = [r for r in rows if "error" not in r]
    for n, _, _ in positions:
        stamp_sheet(good, n, args.out / f"stamps_{n}.png")
    cols = ("name", "obs_id", "program", "band", "mjd", "flux_njy", "err_njy", "snr")
    t = Table(rows=[[r[c] for c in cols] for r in good], names=cols)
    t.sort(["name", "band", "mjd"])
    for c in ("mjd", "flux_njy", "err_njy", "snr"):
        t[c].info.format = ".2f"
    t.meta["provenance"] = "derived"
    t.meta["source"] = "hit_lightcurve.py forced photometry on public NIRCam level-3 i2d (S3)"
    t.meta["params_assumption"] = {
        "radius_arcsec": args.radius_arcsec,
        "annulus_arcsec": ANNULUS,
        "n_random": N_RANDOM,
    }
    t.write(args.out / "lightcurves.ecsv", overwrite=True)
    errors = [r for r in rows if "error" in r]
    (args.out / "skipped.json").write_text(json.dumps(errors, indent=1) + "\n")
    print(f"{len(t)} measurements, {len(errors)} skipped; {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
