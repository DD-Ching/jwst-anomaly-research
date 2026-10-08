"""W5 count-deficit screen (D-063): fetch, predict, screen, inject, vet, limit.

A negative-mass lens (n = 1, ε < 0) would leave an under-dense disc of background galaxies inside
about θ_E (``exotic_sim.count_ratio``; docs/exotic_lensing.md, W5). Subcommands, in order (outputs
under ``$JWST_ANOMALY_DATA/derived/w5_counts/``; the small result tables are copied to
``results/w5_counts/``):

- ``fetch``: per-pixel galaxy counts, mask and depth for the two DES-area regions (``REGIONS``)
  from Legacy Surveys DR10 via Data Lab TAP (``countmap.LegacySurveysCountMap``).
- ``numcounts``: the galaxy number counts N(< r) of the same selection (``observed``), which set
  the predicted profile through magnification bias.
- ``predict``: the deficit profile N_obs/N̄ vs x = θ/θ_E (``model_prediction``) and θ_E vs |M| for
  three lens geometries (``derived``; the geometry is an ASSUMPTION).
- ``screen``: matched filter on every scale in both regions; each region's null (its peak-Z
  distribution) sets the threshold for the other region; flags above threshold.
- ``inject``: the predicted deficit painted into the real maps at random positions (``simulated``),
  through the same screen and the automated vetting; recovery efficiency vs θ_E.
- ``vet``: flags through the ordinary explanations, cheapest first (mask, depth, dust, edge, bright
  stars, galaxy clusters).
- ``limit``: 95 % upper limit on the sky density of such lenses vs θ_E and |M|.

Exotic physics is a hypothesis: a flag is an anomaly to vet, never a discovery.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

from jwst_anomaly import countmap as cm
from jwst_anomaly import exotic_sim as es
from jwst_anomaly import paths, schema

# Two disjoint 10° × 10° DES-wide regions of DR10 south (halved from 20° × 10° so one cloud run can
# fetch both, ~1.7 chunks/min from Data Lab) (uniform ~10-epoch depth, r_5σ,gal ≈ 24.8, b < −45°,
# E(B−V) ≲ 0.04). Each calibrates the null of the other (ASSUMPTION: statistically alike).
REGIONS = {
    "desA": cm.Region("desA", 20.0, 30.0, -30.0, -20.0),
    "desB": cm.Region("desB", 50.0, 60.0, -30.0, -20.0),
}
MAG_LIM = 23.5
CELL = 0.25  # arcmin, base raster cell
SCALES = (2.0, 4.0, 8.0, 16.0, 32.0)  # filter θ_E in arcmin
X_MAX = 2.5  # aperture radius in θ_E (profile within 0.6 % of 1 beyond, model_prediction)
# Pixel acceptance (ASSUMPTIONs): unmasked fraction, 5σ galaxy depth margin over the r < 23.5
# limit (≥ 0.8 mag), dust.
W_MIN, DEPTH_MIN, EBV_MAX = 0.5, 24.3, 0.1
FULL_MIN, COVER_MIN = 0.8, 0.8  # aperture and inner-disc coverage for a usable filter position
# Number-count samples (RA/Dec boxes inside desA, deg) and the magnitude below which each is used.
COUNT_BOXES = (
    ("faint", (30.0, 31.0, -25.0, -24.0), 23.5),
    ("mid", (30.0, 34.0, -26.0, -22.0), 20.5),
    ("bright", (24.0, 36.0, -29.0, -21.0), 18.0),
)
COUNT_GRID = np.arange(16.0, MAG_LIM + 0.01, 0.5)
# Lens geometries for θ_E(|M|) (ASSUMPTIONs): Planck18, background galaxies at z_s = 1.
GEOMETRIES = (("D_L=1 kpc", 1e-3), ("D_L=1 Mpc", 1.0), ("z_L=0.3", None))
Z_SOURCE = 1.0
INJ_THETA = (2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0, 32.0)  # injected θ_E, arcmin
SEED = 20261008


def out_dir() -> Path:
    d = paths.data_root() / "derived" / "w5_counts"
    d.mkdir(parents=True, exist_ok=True)
    return d


def results_dir() -> Path:
    d = paths.repo_root() / "results" / "w5_counts"
    d.mkdir(parents=True, exist_ok=True)
    return d


def stamp(t: Table, provenance: schema.Provenance, source: str) -> Table:
    t.meta["provenance"] = provenance.value
    t.meta["source"] = source
    return t


# ----------------------------------------------------------------------------- fetch / counts


def cmd_fetch(args) -> None:
    """Fetch every missing chunk; failed chunks are retried in later passes (the service 502s)."""

    def one(job):
        survey, box = job
        try:
            return survey.fetch_chunk(box).name
        except Exception as exc:  # report and leave for the next pass
            return f"FAILED {box}: {type(exc).__name__}: {str(exc)[:120]}"

    t0 = time.time()
    for attempt in range(args.passes):
        jobs = []
        for name in args.regions:
            survey = cm.LegacySurveysCountMap(REGIONS[name], MAG_LIM)
            jobs += [(survey, b) for b in survey.missing_chunks()]
        if args.reverse:  # a second process can work from the other end
            jobs = jobs[::-1]
        print(f"pass {attempt + 1}: {len(jobs)} chunks to fetch", flush=True)
        if not jobs:
            return
        try:  # cheap probe: wait out service outages instead of failing every chunk
            cm._run_tap(f"SELECT TOP 1 ra FROM {cm.TABLE} WHERE ra BETWEEN 30 AND 30.01", 1)
        except Exception as exc:
            print(f"  service down ({type(exc).__name__}); waiting {args.pause:.0f}s", flush=True)
            time.sleep(args.pause)
            continue
        with ThreadPoolExecutor(args.workers) as ex:
            for i, msg in enumerate(ex.map(one, jobs)):
                print(f"  {i + 1}/{len(jobs)} {msg} {time.time() - t0:.0f}s", flush=True)
        time.sleep(args.pause)
    raise SystemExit("chunks still missing after all passes")


def cmd_numcounts(args) -> None:
    """Cumulative counts N(< r) per deg² of the screen's galaxy selection (``observed``)."""
    rows = {}
    for tag, (r0, r1, d0, d1), mlim in COUNT_BOXES:
        q = (
            f"SELECT dered_mag_r FROM {cm.TABLE} WHERE ra>={r0} AND ra<{r1} AND dec>={d0} "
            f"AND dec<{d1} AND brick_primary=1 AND {cm.galaxy_selection(mlim)}"
        )
        mags = np.asarray(cm._run_tap(q)["dered_mag_r"], float)
        area = cm.Region(tag, r0, r1, d0, d1).area_deg2()
        rows[tag] = (mags, area, mlim)
        print(f"{tag}: {mags.size} galaxies in {area:.2f} deg²", flush=True)
    n, src = [], []
    for m in COUNT_GRID:
        tag = next(t for t, _, ml in COUNT_BOXES[::-1] if m <= ml)
        mags, area, _ = rows[tag]
        n.append(np.sum(mags < m) / area)
        src.append(tag)
    t = Table({"mag_r": COUNT_GRID, "n_brighter_deg2": np.array(n), "sample": src})
    t.meta["boxes"] = json.dumps({k: [list(b), ml] for k, b, ml in COUNT_BOXES})
    stamp(t, schema.Provenance.OBSERVED, f"{cm.REFERENCE}; {cm.galaxy_selection(MAG_LIM)}")
    t.write(results_dir() / "numcounts.ecsv", overwrite=True)
    print(t)


def load_counts() -> callable:
    t = Table.read(results_dir() / "numcounts.ecsv")
    return cm.tabulated_counts(t["mag_r"], t["n_brighter_deg2"])


def profile_fn():
    counts = load_counts()
    x = np.concatenate([np.linspace(1e-3, 0.999, 2000), np.linspace(1.001, X_MAX + 0.5, 2000)])
    y = cm.deficit_profile(x, counts)
    return lambda r: np.interp(r, x, y, left=float(y[0]), right=1.0)


# ----------------------------------------------------------------------------- predict


def theta_e_arcmin(mass_msun, geometry) -> np.ndarray:
    from astropy.cosmology import Planck18

    name, d_l_mpc = geometry
    mpc = es.PC_M * 1e6
    if d_l_mpc is None:  # cosmological lens
        z_l = 0.3
        d_l = Planck18.angular_diameter_distance(z_l).to_value("Mpc") * mpc
        d_s = Planck18.angular_diameter_distance(Z_SOURCE).to_value("Mpc") * mpc
        d_ls = Planck18.angular_diameter_distance(z_l, Z_SOURCE).to_value("Mpc") * mpc
    else:  # local lens, D_S >> D_L so D_LS / D_S ≈ 1
        d_l, d_s, d_ls = d_l_mpc * mpc, 1.0, 1.0
    m = np.atleast_1d(np.asarray(mass_msun, float))
    th = np.array([es.einstein_radius(es.eps_bar_point_mass(mi), 1.0, d_l, d_s, d_ls) for mi in m])
    return np.degrees(th) * 60.0


def mass_for_theta(theta_arcmin, geometry) -> np.ndarray:
    """|M| (M☉) with θ_E = theta (θ_E ∝ M^½ for n = 1)."""
    return (np.asarray(theta_arcmin, float) / theta_e_arcmin(1.0, geometry)[0]) ** 2


def cmd_predict(args) -> None:
    counts = load_counts()
    x = np.array(
        [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.05, 1.1, 1.2, 1.5, 2.0]
    )
    prof = Table({"x": x, "ratio_measured_counts": cm.deficit_profile(x, counts)})
    for a in (0.6, 1.5):  # power-law counts, alpha = 2.5 dlogN/dm
        prof[f"ratio_alpha{a}"] = cm.deficit_profile(x, cm.power_counts(a))
    stamp(
        prof,
        schema.Provenance.MODEL_PREDICTION,
        "exotic_sim.count_ratio (n = 1, eps < 0) with results/w5_counts/numcounts.ecsv",
    )
    prof.write(results_dir() / "profile.ecsv", overwrite=True)
    print(prof)
    dens = float(Table.read(results_dir() / "numcounts.ecsv")["n_brighter_deg2"][-1]) / 3600.0
    pf = profile_fn()
    rows = []
    for th in (10 / 60, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 60.0):
        miss1, tot1 = cm.expected_missing(th, dens, pf, 1.0)
        miss, tot = cm.expected_missing(th, dens, pf, X_MAX)
        row = {"theta_e_arcmin": th, "n_expected_x_lt_1": tot1, "n_missing_x_lt_1": miss1}
        row["n_missing_x_lt_2.5"] = miss
        for g in GEOMETRIES:
            row[f"mass_msun {g[0]}"] = float(mass_for_theta(th, g))
        rows.append(row)
    t = Table(rows=rows)
    stamp(
        t,
        schema.Provenance.DERIVED,
        f"exotic_sim.einstein_radius, Planck18, z_s = {Z_SOURCE}; geometries are ASSUMPTIONs",
    )
    t.write(results_dir() / "theta_mass.ecsv", overwrite=True)
    t.pprint(max_width=200)


# ----------------------------------------------------------------------------- maps


class MapState:
    """One region's count map on the raster with its pixel geometry (cached in derived/)."""

    def __init__(self, name: str):
        from astropy import units as u
        from astropy_healpix import HEALPix

        self.name = name
        self.region = REGIONS[name]
        t = cm.LegacySurveysCountMap(self.region, MAG_LIM).count_map()
        t.sort("pix")
        self.table = t
        self.raster = cm.make_raster(self.region, CELL)
        cache = out_dir() / f"idx_{name}_{len(t)}_{CELL}.npy"  # keyed by the pixel list
        if cache.exists():
            self.idx = np.load(cache)
        else:
            self.idx = cm.pixel_index_image(self.raster, np.asarray(t["pix"], np.int64))
            np.save(cache, self.idx)
        lon, lat = HEALPix(nside=cm.NSIDE, order="nested").healpix_to_lonlat(t["pix"])
        self.pix_ra, self.pix_dec = lon.to_value(u.deg), lat.to_value(u.deg)
        self.pix_xy = self.raster.xy(self.pix_ra, self.pix_dec)
        good = (
            (np.asarray(t["w"]) >= W_MIN)
            & (np.asarray(t["depth_r"]) >= DEPTH_MIN)
            & (np.asarray(t["ebv"]) <= EBV_MAX)
        )
        self.pix_weight = np.where(good, np.asarray(t["w"], float), 0.0)
        self.weight = cm.rasterise(self.pix_weight, self.idx)
        self.n_gal = np.asarray(t["n_gal"], np.int64)
        self.area_deg2 = float(np.sum(self.pix_weight) * cm.PIX_AREA_DEG2)

    def density(self, n_gal=None) -> np.ndarray:
        n = self.n_gal if n_gal is None else n_gal
        return cm.rasterise(np.where(self.pix_weight > 0, n, 0), self.idx, per_cell=True)


def factor(scale: float) -> int:
    return max(1, int(scale / (8 * CELL)))


class Screen:
    """Matched filters of one region on all scales (weights precomputed)."""

    def __init__(self, ms: MapState, pf):
        self.ms = ms
        self.filters = {}
        for s in SCALES:
            f = factor(s)
            wb = cm.block_sum(ms.weight, f) / f**2
            self.filters[s] = (f, cm.MatchedFilter(wb, s / (CELL * f), pf, X_MAX))

    def run(self, density: np.ndarray) -> dict[float, dict[str, np.ndarray]]:
        out = {}
        for s, (f, mf) in self.filters.items():
            amp, sig = mf(cm.block_sum(density, f))  # counts per block (Poisson σ)
            ok = (mf.full >= FULL_MIN) & (mf.cover >= COVER_MIN) & np.isfinite(amp) & (sig > 0)
            with np.errstate(divide="ignore", invalid="ignore"):
                z = np.where(ok, amp / sig, np.nan)
            out[s] = {"amp": amp, "sigma": sig, "z": z, "ok": ok}
        return out

    def peaks(self, res, s: float, z_min: float) -> Table:
        f, mf = self.filters[s]
        r = res[s]
        iy, ix = cm.find_peaks(r["z"], r["ok"], int(round(s / (CELL * f))), z_min)
        xc = self.ms.raster.x0 + (ix + 0.5) * CELL * f
        yc = self.ms.raster.y0 + (iy + 0.5) * CELL * f
        ra, dec = self.ms.raster.radec(xc, yc)
        return Table(
            {
                "scale_arcmin": np.full(ix.size, s),
                "x": xc,
                "y": yc,
                "ra": ra,
                "dec": dec,
                "z": r["z"][iy, ix],
                "amp": r["amp"][iy, ix],
                "sigma": r["sigma"][iy, ix],
                "cover": mf.cover[iy, ix],
            }
        )


# ----------------------------------------------------------------------------- screen


def cmd_screen(args) -> None:
    pf = profile_fn()
    peaks, summary = {}, {}
    for name in REGIONS:
        t0 = time.time()
        ms = MapState(name)
        sc = Screen(ms, pf)
        res = sc.run(ms.density())
        tabs = []
        for s in SCALES:
            p = sc.peaks(res, s, -np.inf)
            p["region"] = name
            tabs.append(p)
            z = res[s]["z"]
            summary[f"{name} {s}"] = {
                "n_usable_cells": int(np.sum(res[s]["ok"])),
                "amp_median": float(np.nanmedian(res[s]["amp"][res[s]["ok"]])),
                "amp_std": float(np.nanstd(res[s]["amp"][res[s]["ok"]])),
                "z_std": float(np.nanstd(z)),
                "z_max": float(np.nanmax(z)),
            }
        peaks[name] = vstack(tabs)
        summary[name] = {
            "area_deg2": ms.area_deg2,
            "box_area_deg2": ms.region.area_deg2(),
            "n_gal": int(np.sum(ms.n_gal[ms.pix_weight > 0])),
            "pixels": len(ms.table),
            "wall_s": time.time() - t0,
        }
        print(name, summary[name], flush=True)
    # Null calibration: each region's peaks (same scale) calibrate the other region (D-063).
    names = list(REGIONS)
    null = {}
    for name, other in ((names[0], names[1]), (names[1], names[0])):
        for s in SCALES:
            zo = np.asarray(peaks[other]["z"][peaks[other]["scale_arcmin"] == s], float)
            n_search = int(np.sum(peaks[name]["scale_arcmin"] == s))
            null[f"{name} {s}"] = null_model(zo, n_search)
    allp = vstack([peaks[n] for n in names])
    key = [f"{r} {s}" for r, s in zip(allp["region"], allp["scale_arcmin"], strict=True)]
    allp["z_flag"] = [null[k]["z_flag"] for k in key]
    allp["n_false_expected"] = [n_false(null[k], z) for k, z in zip(key, allp["z"], strict=True)]
    allp["flag"] = allp["z"] > allp["z_flag"]
    stamp(
        allp,
        schema.Provenance.DERIVED,
        "w5_counts screen: local maxima of the W5 matched-filter Z per scale (LS DR10)",
    )
    allp.write(out_dir() / "peaks.ecsv", overwrite=True)
    flags = allp[allp["flag"]]
    flags.write(results_dir() / "flags.ecsv", overwrite=True)
    summary["null"] = null
    summary["n_flags"] = len(flags)
    (results_dir() / "screen_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(null, indent=1))
    flags.pprint(max_width=200, max_lines=60)


FAP = 0.01  # expected false peaks per region and scale for a detection (ASSUMPTION)
TAIL_Q = 0.01  # fraction of null peaks used for the exponential tail fit (at least 30)


def null_model(z_null: np.ndarray, n_search: int) -> dict:
    """Empirical null of the peak Z from an independent region (no Gaussian assumption).

    ``z_flag`` is the null region's highest peak (anything above it is flagged for vetting).
    The upper tail of the null peak distribution is fitted with an exponential (MLE over the top
    ``TAIL_Q`` of peaks; an ASSUMPTION about the shape beyond the highest null peak) to give the
    expected number of null peaks above z among ``n_search`` peaks, and ``z_det`` where that
    number is ``FAP``.
    """
    z = np.sort(z_null[np.isfinite(z_null)])[::-1]
    q = min(max(30, int(TAIL_Q * z.size)), z.size - 1)
    z0 = float(z[q])
    tau = float(np.mean(z[:q] - z0))
    rate0 = q / z.size
    z_det = z0 + tau * math.log(max(n_search * rate0 / FAP, 1.0))
    return {
        "z_flag": float(z[0]),
        "z0": z0,
        "tau": tau,
        "rate0": rate0,
        "n_null_peaks": int(z.size),
        "n_search_peaks": int(n_search),
        "z_det": float(max(z_det, z[0])),
    }


def n_false(model: dict, z: float) -> float:
    """Expected number of null peaks with Z ≥ z among the search region's peaks."""
    return float(
        model["n_search_peaks"] * model["rate0"] * math.exp(-(z - model["z0"]) / model["tau"])
    )


# ----------------------------------------------------------------------------- injection


def place_centres(ms: MapState, theta: float, rng, n_max: int) -> np.ndarray:
    """Random footprint positions (raster cells with weight > 0), ≥ 2 (X_MAX + 1) θ_E apart."""
    cells = np.flatnonzero(ms.weight.ravel() > 0)
    rng.shuffle(cells)
    min_sep = 2 * (X_MAX + 1) * max(theta, max(SCALES) / 4)
    chosen: list[tuple[float, float]] = []
    for c in cells[: 50 * n_max]:
        iy, ix = divmod(int(c), ms.raster.nx)
        x = ms.raster.x0 + (ix + 0.5) * CELL
        y = ms.raster.y0 + (iy + 0.5) * CELL
        if all(math.hypot(x - a, y - b) >= min_sep for a, b in chosen):
            chosen.append((x, y))
            if len(chosen) >= n_max:
                break
    return np.array(chosen)


def _inject_job(job) -> Table:
    name, theta, rep, null = job
    pf = profile_fn()
    ms = MapState(name)
    sc = Screen(ms, pf)
    rng = np.random.default_rng([SEED, rep, int(theta * 10), list(REGIONS).index(name)])
    centres = place_centres(ms, theta, rng, 200)
    n = ms.n_gal.copy()
    for xy in centres:
        n = cm.inject_deficit(n, ms.pix_xy, tuple(xy), theta, pf, rng, X_MAX)
    res = sc.run(ms.density(n))
    found = []
    for s in SCALES:
        p = sc.peaks(res, s, null[f"{name} {s}"]["z_flag"])
        p["n_false_expected"] = [n_false(null[f"{name} {s}"], z) for z in p["z"]]
        found.append(p)
    found = vstack(found)
    vet = Vetter(ms)
    rows = []
    for x, y in centres:
        d = np.hypot(found["x"] - x, found["y"] - y) if len(found) else np.array([])
        hit = d <= max(0.5 * theta, 1.0)
        # the hit that is most significant against its own scale's null
        best = found[hit][int(np.argmin(found["n_false_expected"][hit]))] if hit.any() else None
        v = None
        if best is not None:
            v = vet(best["x"], best["y"], best["scale_arcmin"], best["n_false_expected"])
        rows.append(
            {
                "region": name,
                "theta_e_arcmin": theta,
                "rep": rep,
                "x": x,
                "y": y,
                "recovered": bool(hit.any()),
                "best_scale": float(best["scale_arcmin"]) if best is not None else np.nan,
                "best_z": float(best["z"]) if best is not None else np.nan,
                "best_n_false": float(best["n_false_expected"]) if best is not None else np.nan,
                "vetoed_by": v["reasons"] if v is not None else "",
                "detected": bool(v is not None and v["survives"]),
            }
        )
    return Table(rows=rows)


def cmd_inject(args) -> None:
    summary = json.loads((results_dir() / "screen_summary.json").read_text())
    null = summary["null"]
    jobs = [(n, th, rep, null) for n in REGIONS for th in INJ_THETA for rep in range(args.reps)]
    t0 = time.time()
    tabs = []
    with Pool(args.workers) as pool:
        for i, t in enumerate(pool.imap_unordered(_inject_job, jobs)):
            tabs.append(t)
            print(f"{i + 1}/{len(jobs)} {time.time() - t0:.0f}s", flush=True)
    inj = vstack(tabs)
    stamp(
        inj,
        schema.Provenance.SIMULATED,
        "W5 deficits (exotic_sim.count_ratio) injected into LS DR10 count maps",
    )
    inj.write(out_dir() / "injections.ecsv", overwrite=True)
    print(f"wall {time.time() - t0:.0f}s, {len(inj)} injections")


# ----------------------------------------------------------------------------- vet

# Vetting thresholds (ASSUMPTIONs, cheapest first; the same rules veto injections).
VET_MASK_FRAC = 0.10  # masked or missing area inside θ_E
VET_DEPTH_DROP = 0.15  # mag shallower than the region median inside θ_E
VET_EBV_EXCESS = 0.02  # mag above the region median inside θ_E
VET_STAR_G = 9.0  # Gaia DR3 G of a star whose halo/mask could empty the disc
VET_STAR_MIN_R = 2.0  # arcmin: stars this close always veto, whatever θ_E
VET_CLUSTER_M500 = 2.0  # 10¹⁴ M☉: clusters massive enough for magnification-bias depletion
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap"
WH24 = "J/ApJS/272/39/table2"  # Wen & Han 2024 DESI Legacy Surveys clusters (VizieR)


def ancillary(name: str) -> tuple[Table, Table]:
    """Bright Gaia DR3 stars and massive WH24 clusters in a region (cached, ``observed``)."""
    reg = REGIONS[name]
    p_star = out_dir() / f"stars_{name}.ecsv"
    p_cl = out_dir() / f"clusters_{name}.ecsv"
    if not p_star.exists():
        import pyvo

        q = (
            "SELECT source_id, ra, dec, phot_g_mean_mag FROM gaiadr3.gaia_source WHERE "
            f"phot_g_mean_mag < {VET_STAR_G} AND ra BETWEEN {reg.ra_min - 1} AND "
            f"{reg.ra_max + 1} AND dec BETWEEN {reg.dec_min - 1} AND {reg.dec_max + 1}"
        )
        t = pyvo.dal.TAPService(GAIA_TAP).run_sync(q, maxrec=100000).to_table()
        stamp(t, schema.Provenance.OBSERVED, f"Gaia DR3 via {GAIA_TAP}: {q}")
        t.write(p_star, overwrite=True)
    if not p_cl.exists():
        from astropy import units as u
        from astropy.coordinates import SkyCoord
        from astroquery.vizier import Vizier

        v = Vizier(
            columns=["Name", "RAJ2000", "DEJ2000", "zCl", "M500", "lam500"],
            column_filters={"M500": f">={VET_CLUSTER_M500}"},
            row_limit=-1,
        )
        c = SkyCoord(0.5 * (reg.ra_min + reg.ra_max), 0.5 * (reg.dec_min + reg.dec_max), unit="deg")
        res = v.query_region(
            c,
            width=(reg.ra_max - reg.ra_min + 2) * math.cos(math.radians(reg.dec_max)) * u.deg,
            height=(reg.dec_max - reg.dec_min + 2) * u.deg,
            catalog=WH24,
        )
        t = res[0] if len(res) else Table(names=("RAJ2000", "DEJ2000", "zCl", "M500"))
        t = Table(t, masked=False)
        stamp(t, schema.Provenance.OBSERVED, f"VizieR {WH24} (Wen & Han 2024), M500 ≥ 2e14")
        t.write(p_cl, overwrite=True)
    return Table.read(p_star), Table.read(p_cl)


class Vetter:
    """Ordinary explanations for a deficit at (x, y) with scale θ_E (arcmin), cheapest first."""

    def __init__(self, ms: MapState):
        from scipy.spatial import cKDTree

        self.ms = ms
        t = ms.table
        self.tree = cKDTree(np.column_stack(ms.pix_xy))
        self.w = np.asarray(t["w"], float)
        self.depth = np.asarray(t["depth_r"], float)
        self.ebv = np.asarray(t["ebv"], float)
        self.depth_med = float(np.median(self.depth))
        self.ebv_med = float(np.median(self.ebv))
        stars, clusters = ancillary(ms.name)
        self.star_xy = np.column_stack(ms.raster.xy(stars["ra"], stars["dec"]))
        self.star_g = np.asarray(stars["phot_g_mean_mag"], float)
        self.cl_xy = np.column_stack(ms.raster.xy(clusters["RAJ2000"], clusters["DEJ2000"]))
        self.pix_side = math.sqrt(cm.PIX_AREA_DEG2) * 60.0

    def __call__(self, x: float, y: float, theta: float, n_false_expected: float) -> dict:
        r = max(theta, self.pix_side)
        i = self.tree.query_ball_point([x, y], r)
        n_exp = math.pi * r**2 / self.pix_side**2
        mask_frac = 1.0 - (np.sum(self.w[i]) / max(n_exp, len(i))) if i else 1.0
        depth = float(np.mean(self.depth[i])) if i else np.nan
        ebv = float(np.mean(self.ebv[i])) if i else np.nan
        d_star = np.hypot(self.star_xy[:, 0] - x, self.star_xy[:, 1] - y)
        near = d_star <= max(theta, VET_STAR_MIN_R)
        d_cl = np.hypot(self.cl_xy[:, 0] - x, self.cl_xy[:, 1] - y) if len(self.cl_xy) else []
        out = {
            "mask_frac": mask_frac,
            "depth_r": depth,
            "ebv": ebv,
            "n_bright_star": int(near.sum()),
            "brightest_g": float(self.star_g[near].min()) if near.any() else np.nan,
            "n_cluster": int(np.sum(np.asarray(d_cl) <= theta)),
        }
        reasons = []
        if mask_frac > VET_MASK_FRAC:
            reasons.append("mask")
        if not depth >= self.depth_med - VET_DEPTH_DROP:
            reasons.append("depth")
        if not ebv <= self.ebv_med + VET_EBV_EXCESS:
            reasons.append("dust")
        if out["n_bright_star"]:
            reasons.append("bright_star")
        if out["n_cluster"]:
            reasons.append("cluster")
        if not n_false_expected < FAP:  # no deeper than the other region's null peaks reach
            reasons.append("cosmic_variance")
        out["reasons"] = ",".join(reasons)
        out["survives"] = not reasons
        return out


def cmd_vet(args) -> None:
    flags = Table.read(results_dir() / "flags.ecsv")
    rows = []
    vetters = {}
    for f in flags:
        v = vetters.setdefault(f["region"], Vetter(MapState(f["region"])))
        rows.append(
            {
                "region": f["region"],
                "ra": f["ra"],
                "dec": f["dec"],
                "scale_arcmin": f["scale_arcmin"],
                "z": f["z"],
                "amp": f["amp"],
                "n_false_expected": f["n_false_expected"],
                **v(f["x"], f["y"], f["scale_arcmin"], f["n_false_expected"]),
            }
        )
    t = Table(rows=rows) if rows else Table(names=("region", "survives"), dtype=(str, bool))
    t.meta["thresholds"] = json.dumps(
        {
            "mask_frac": VET_MASK_FRAC,
            "depth_drop": VET_DEPTH_DROP,
            "ebv_excess": VET_EBV_EXCESS,
            "star_g": VET_STAR_G,
            "star_min_r_arcmin": VET_STAR_MIN_R,
            "cluster_m500_1e14": VET_CLUSTER_M500,
        }
    )
    stamp(t, schema.Provenance.DERIVED, "w5_counts vet: ordinary explanations for each flag")
    t.write(results_dir() / "vetting.ecsv", overwrite=True)
    t.pprint(max_width=250, max_lines=100)
    print(f"{len(t)} flags, {int(np.sum(t['survives'])) if len(t) else 0} survive")


def cmd_sheet(args) -> None:
    """Contact sheet of the top flags: galaxy density, unmasked fraction, LS DR10 colour image."""
    import io
    import urllib.request

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    vet = Table.read(results_dir() / "vetting.ecsv")
    vet.sort("n_false_expected")
    vet = vet[: args.n]
    states = {}
    fig, axes = plt.subplots(len(vet), 3, figsize=(10, 3.2 * len(vet)), squeeze=False)
    for row, ax in zip(vet, axes, strict=True):
        ms = states.setdefault(row["region"], MapState(row["region"]))
        th = float(row["scale_arcmin"])
        x, y = ms.raster.xy(row["ra"], row["dec"])
        half = 3 * th
        f = max(1, int(th / 4 / CELL))
        ix0 = int((x - half - ms.raster.x0) / CELL)
        iy0 = int((y - half - ms.raster.y0) / CELL)
        n = int(2 * half / CELL)
        sl = (slice(max(iy0, 0), iy0 + n), slice(max(ix0, 0), ix0 + n))
        dens = cm.block_sum(ms.density()[sl], f)
        wt = cm.block_sum(ms.weight[sl], f) / f**2
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = dens / wt / np.nanmean(dens[wt > 0.5] / wt[wt > 0.5])
        ext = [-half, half, -half, half]
        im = ax[0].imshow(
            np.where(wt > 0.3, rel, np.nan), origin="lower", extent=ext, cmap="RdBu", vmin=0, vmax=2
        )
        fig.colorbar(im, ax=ax[0], fraction=0.046)
        ax[1].imshow(wt, origin="lower", extent=ext, cmap="gray", vmin=0, vmax=1)
        for a in ax[:2]:
            a.add_patch(plt.Circle((0, 0), th, fill=False, color="k"))
        ax[0].set_title(
            f"{row['region']} θ={th:g}′ Z={row['z']:.1f} N_false={row['n_false_expected']:.2g}",
            fontsize=8,
        )
        ax[1].set_title(f"w; {row['reasons'] or 'survives'}", fontsize=8)
        pixscale = max(0.262, 2 * th * 60 / 512)
        url = (
            "https://www.legacysurvey.org/viewer/cutout.jpg?ra="
            f"{row['ra']:.5f}&dec={row['dec']:.5f}&layer=ls-dr10&pixscale={pixscale:.3f}&size=512"
        )
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                ax[2].imshow(plt.imread(io.BytesIO(r.read()), format="jpg"))
        except Exception as exc:  # cutout service down: leave the panel empty
            ax[2].text(0.1, 0.5, f"cutout failed: {type(exc).__name__}", fontsize=7)
        ax[2].set_title(f"{row['ra']:.4f} {row['dec']:.4f} (2θ wide)", fontsize=8)
        ax[2].axis("off")
    fig.tight_layout()
    out = out_dir() / "contact_sheet.png"
    fig.savefig(out, dpi=80)
    print(out)


# ----------------------------------------------------------------------------- limit


def cmd_limit(args) -> None:
    inj = Table.read(out_dir() / "injections.ecsv")
    summary = json.loads((results_dir() / "screen_summary.json").read_text())
    vet = (
        Table.read(results_dir() / "vetting.ecsv")
        if (results_dir() / "vetting.ecsv").exists()
        else None
    )
    n_surv = int(np.sum(vet["survives"])) if vet is not None and len(vet) else 0
    area = sum(summary[n]["area_deg2"] for n in REGIONS)
    # Poisson 95 % upper limit on the expected number for n observed (0 → 2.996)
    from scipy.stats import chi2

    n95 = 0.5 * chi2.ppf(0.95, 2 * (n_surv + 1))
    rows = []
    for th in INJ_THETA:
        sel = inj["theta_e_arcmin"] == th
        n = int(np.sum(sel))
        k_rec = int(np.sum(inj["recovered"][sel]))
        k = int(np.sum(inj["detected"][sel]))
        eff = k / n
        row = {
            "theta_e_arcmin": th,
            "n_injected": n,
            "n_screen": k_rec,
            "n_after_vetting": k,
            "efficiency": eff,
            "density_95_per_deg2": n95 / (area * eff) if eff > 0 else np.inf,
        }
        for g in GEOMETRIES:
            row[f"mass_msun {g[0]}"] = float(mass_for_theta(th, g))
        rows.append(row)
    t = Table(rows=rows)
    t.meta["area_deg2"] = area
    t.meta["n_survivors"] = n_surv
    t.meta["n95"] = n95
    stamp(t, schema.Provenance.DERIVED, "W5 injection-recovery through the w5_counts screen")
    t.write(results_dir() / "limits.ecsv", overwrite=True)
    t.pprint(max_width=200)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fetch")
    f.add_argument("--regions", nargs="+", default=list(REGIONS))
    f.add_argument("--workers", type=int, default=2)
    f.add_argument("--passes", type=int, default=200)
    f.add_argument("--pause", type=float, default=120.0)
    f.add_argument("--reverse", action="store_true")
    f.set_defaults(func=cmd_fetch)
    sub.add_parser("numcounts").set_defaults(func=cmd_numcounts)
    sub.add_parser("predict").set_defaults(func=cmd_predict)
    sub.add_parser("screen").set_defaults(func=cmd_screen)
    i = sub.add_parser("inject")
    i.add_argument("--reps", type=int, default=2)
    i.add_argument("--workers", type=int, default=3)
    i.set_defaults(func=cmd_inject)
    sub.add_parser("vet").set_defaults(func=cmd_vet)
    sh = sub.add_parser("sheet")
    sh.add_argument("--n", type=int, default=12)
    sh.set_defaults(func=cmd_sheet)
    sub.add_parser("limit").set_defaults(func=cmd_limit)
    args = ap.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
