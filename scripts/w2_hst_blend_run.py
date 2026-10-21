"""Line C measurement stages (frozen plan: results/w2_hst_blend/PLAN.md).

    python -I scripts/w2_hst_blend_run.py download --cut <dir> --log <file>   # HAP cutouts
    python -I scripts/w2_hst_blend_run.py measure --role control --cut <dir>  # controls + gate
    python -I scripts/w2_hst_blend_run.py psf --cut <dir> --log <file>        # PSF-mismatch sims
    python -I scripts/w2_hst_blend_run.py measure --role target --cut <dir>   # only if gate passed
    python -I scripts/w2_hst_blend_run.py sheet

HST pixels are ``observed``; fits and residuals ``derived``; injections and the PSF-mismatch
stamps ``simulated``. Nothing here is evidence of a dark or exotic lens.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w2_hst_blend_chain as ch  # noqa: E402
import w12_niq_hst as hst  # noqa: E402
from w2_hst_blend import RES, log_download  # noqa: E402

from jwst_anomaly import schema  # noqa: E402

DOWNLOAD_STATUS = "download_status.json"
NO_HAP = "no HAP product (service answered: no HAP files)"
FAILED = "download failed"


def hapcut_with_retry(c: SkyCoord, size: int, path: Path, tries: int = 3):
    """Hapcut cutouts with retries and back-off. Returns (table or None, status, message): a
    service error after the retries is "download failed", never "no data"; an explicit empty
    answer (astroquery ``NoResultsWarning``, no files) is ``NO_HAP``."""
    import warnings

    from astroquery.mast import Hapcut

    msg = ""
    for attempt in range(tries):
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                t = Hapcut.download_cutouts(c, size=[size, size], path=str(path))
            warn = "; ".join(f"{type(x.message).__name__}: {x.message}" for x in w)
            if t is None or len(t) == 0:
                if "NoResults" in warn:
                    return None, NO_HAP, warn
                msg = f"empty result without NoResultsWarning ({warn})"
            else:
                return t, "ok", warn
        except Exception as e:  # noqa: BLE001 - retried, then recorded as a failure
            msg = f"{type(e).__name__}: {e}"
        print("hapcut retry", attempt + 1, msg[:200], flush=True)
        time.sleep(10 * (attempt + 1))
    return None, FAILED, msg[:500]


NOT_ATTEMPTED = "not attempted"


def set_status(sf: Path, status: dict, key: str, value: dict) -> None:
    """Record one key in ``status`` and in the status file. The file is re-read and only ``key``
    is replaced, then written atomically, so a stage never erases another stage's entries."""
    status[key] = value
    disk = json.loads(sf.read_text()) if sf.exists() else {}
    disk[key] = value
    tmp = sf.with_suffix(".tmp")
    tmp.write_text(json.dumps(disk, indent=2) + "\n")
    tmp.replace(sf)


def download_set(sf: Path, status: dict, key: str, name: str, c, size: int, d: Path, log: Path):
    """One cutout set with its status record (shared by ``download`` and ``psf``): a partial or
    failed earlier attempt is removed and redone; "in progress" is recorded first, so an
    interruption never leaves "ok"; the final record overwrites any earlier one."""
    import shutil

    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True, exist_ok=True)
    set_status(sf, status, key, {"name": name, "status": "in progress"})
    t, st, msg = hapcut_with_retry(c, size, d)
    for f in t["Local Path"] if t is not None else []:
        log_download(log, f"hapcut ra={c.ra.deg:.6f} dec={c.dec.deg:.6f} size={size}", f)
    set_status(sf, status, key, {
        "name": name, "status": st, "n_files": len(t) if t is not None else 0, "message": msg,
        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    })  # fmt: skip
    return t, st, msg


def cutout_complete(status: dict, sid: str, d: Path) -> bool:
    """A cutout set counts as complete only when the last attempt for it recorded "ok" (a run
    that raised half-way leaves a partial set without that record, so it is redone)."""
    return status.get(sid, {}).get("status") == "ok" and d.exists() and bool(list(d.glob("*.fits")))


def needs_download(status: dict, key: str, d: Path, force: bool = False) -> bool:
    """Download unless the set is complete ("ok" with files) or the service's explicit no-data
    answer (``NO_HAP``) is on record; ``force`` re-queries a recorded no-data system."""
    if cutout_complete(status, key, d):
        return False
    return not (status.get(key, {}).get("status") == NO_HAP and not force)


def cmd_download(args) -> None:

    p = ch.Params()
    sel = Table.read(RES / "selection.ecsv")
    log = Path(args.log)
    sf = RES / DOWNLOAD_STATUS
    status = json.loads(sf.read_text()) if sf.exists() else {}
    for r in sel:
        sid = str(r["system_id"])
        d = Path(args.cut) / sid
        if not needs_download(status, sid, d, args.force):
            continue  # complete, or a recorded no-data answer (final unless --force)
        c = SkyCoord(float(r["ra"]), float(r["dec"]), unit="deg")
        t, st, msg = download_set(sf, status, sid, str(r["name"]), c, p.size_px, d, log)
        print(r["name"], st, len(t) if t is not None else 0, "files", flush=True)
    bad = [v["name"] for v in status.values() if v["status"] == FAILED]
    if bad:  # fail loudly: a service error is not "no data"
        raise SystemExit(f"hapcut download failed for {bad}; re-run download before measuring")


def missing_cutout_status(status: dict, sid: str) -> str:
    """Status of a system without cutout files: the service's explicit no-data answer, a failed
    download, or never attempted (no record, or an interrupted attempt)."""
    st = status.get(sid, {}).get("status")
    if st == NO_HAP:
        return "no HAP cutout"
    if st == FAILED:
        return FAILED
    return NOT_ATTEMPTED


def load_cutout(d: Path, c: SkyCoord, p: ch.Params):
    files = sorted(d.glob("*.fits"))
    f = hst.pick_image(files, c, p.search_arcsec)
    if f is None:
        return None
    with fits.open(f) as h:
        hdu = next(x for x in h if x.data is not None and x.data.ndim == 2)
        data = hdu.data.astype(float)
        w = WCS(hdu.header)
        zp = hst.ab_zeropoint(hdu.header, h[0].header)
    valid = np.isfinite(data) & (data != 0)  # 0 = outside the drizzled footprint
    img = np.where(valid, data, np.median(data[valid]))
    pix = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix))) * 3600)
    x0, y0 = w.world_to_pixel(c)
    return {
        "file": f.name,
        "data": data,
        "valid": valid,
        "img": img,
        "pix": pix,
        "zp": zp,
        "centre": (float(x0), float(y0)),
    }


def noise_function(data, valid):
    """Per-pixel aperture noise for ``npix`` pixels: the larger of the cutout's pixel MAD x the
    drizzle-correlation factor and the same-area empty-sky aperture scatter / sqrt(npix)
    (w12_niq_hst rule; the MAD is taken on the cutout's pixels, see RESULT.md)."""
    med = np.median(data[valid])
    sig = 1.4826 * np.median(np.abs(data[valid] - med))
    corr = hst.noise_correlation(data, valid)
    cache: dict[int, float] = {}

    def fn(npix: int) -> float:
        if npix not in cache:
            emp = hst.empirical_aperture_noise(data, valid, npix)
            cache[npix] = float(np.nanmax([sig * corr, emp / np.sqrt(npix)]))
        return cache[npix]

    return sig, corr, fn


def moffat_flux(amp, gamma, alpha) -> float:
    return float(amp * np.pi * gamma**2 / (alpha - 1))


def run_system(r, cut: Path, p: ch.Params, do_inject: bool = True):
    c = SkyCoord(float(r["ra"]), float(r["dec"]), unit="deg")
    row = {"system_id": str(r["system_id"]), "name": str(r["name"]), "role": str(r["role"])}
    row.update(status_ls=str(r["status_ls"]), t814_s=float(r["t814_s"]))
    inj = []
    d = cut / str(r["system_id"])
    sf = RES / DOWNLOAD_STATUS
    status = json.loads(sf.read_text()) if sf.exists() else {}
    if not cutout_complete(status, str(r["system_id"]), d):  # files without "ok" are partial
        row["status"] = missing_cutout_status(status, str(r["system_id"]))
        return row, None, inj
    cu = load_cutout(d, c, p)
    if cu is None:
        row["status"] = "no F814W cutout"
        return row, None, inj
    sig, corr, nfn = noise_function(cu["data"], cu["valid"])
    m = ch.measure(cu["img"], cu["valid"], cu["centre"], cu["pix"], p, sig, nfn)
    row.update(
        hst_image=cu["file"],
        pixel_arcsec=cu["pix"],
        zp_ab=cu["zp"],
        sigma_pix=sig,
        noise_corr=corr,
        status=m["status"],
        n_img=m["n_img"],
        max_sep_arcsec=m["max_sep_arcsec"],
        resid_flux=m["resid_flux"],
        resid_err=m["resid_err"],
        resid_snr=m["resid_snr"],
        resid_snr_raw=m["resid_snr_raw"],
        n_pix=m["n_pix"],
        pair_check="unchecked (no catalogued separation)",
        img_pos_px=";".join(f"{x:.2f},{y:.2f}" for x, y in m["pos"]),
        img_amps=";".join(f"{a:.6g}" for a in m["amps"]),
    )
    row["mag_resid"] = cu["zp"] - 2.5 * np.log10(m["resid_flux"]) if m["resid_flux"] > 0 else np.nan
    row["mag_5sigma"] = (
        cu["zp"] - 2.5 * np.log10(p.detect_sigma * m["resid_err"])
        if np.isfinite(m["resid_err"])
        else np.nan
    )
    if m["amps"]:
        row["image_mags"] = ";".join(
            f"{cu['zp'] - 2.5 * np.log10(moffat_flux(a, m['psf_gamma'], m['psf_alpha'])):.2f}"
            for a in m["amps"]
        )
        row.update(psf_gamma=m["psf_gamma"], psf_alpha=m["psf_alpha"])
        pred = ch.predicted_f814w(
            m["max_sep_arcsec"] / 2, float(r["z_lens"]), float(r["z_source"]), p
        )
        row.update(theta_e_sis=m["max_sep_arcsec"] / 2, m_pred_typ=pred["m_typ"])
        row["m_pred_faint"] = pred["m_faint"]
    for reff in p.inj_reff_arcsec:
        row[f"m_lim_reff{reff}"] = np.nan
    if do_inject and m["status"] in ch.DECIDED:
        lx, ly = ch.lens_position(m["pos"], m["amps"])
        psf = ch.moffat_kernel(m["psf_gamma"], m["psf_alpha"])
        for reff in p.inj_reff_arcsec:
            g = ch.sersic_unit(
                cu["img"].shape, lx, ly, reff / cu["pix"], p.inj_n, p.inj_ellip, psf=psf
            )
            res = []
            for mag in p.inj_mags:
                flux = 10 ** (-0.4 * (mag - cu["zp"]))
                mi = ch.measure(
                    cu["img"] + flux * g, cu["valid"], cu["centre"], cu["pix"], p, sig, nfn
                )
                ok = ch.recovered(m, mi, p)
                res.append((mag, ok))
                ff = (mi["resid_flux"] - m["resid_flux"]) / flux
                inj.append(
                    {"system_id": row["system_id"], "name": row["name"], "role": row["role"],
                     "reff_arcsec": reff, "inject_mag": mag, "status": mi["status"],
                     "n_img": mi["n_img"], "resid_snr": mi["resid_snr"],
                     "snr_above_baseline": mi["resid_snr"] - m["resid_snr"],
                     "flux_frac": ff, "recovered": ok, "max_sep_arcsec": m["max_sep_arcsec"]}
                )  # fmt: skip
            row[f"m_lim_reff{reff}"] = ch.contiguous_limit(res)
        lim = row["m_lim_reff0.3"]
        if "m_pred_typ" in row:
            row["delta_typ"] = lim - row["m_pred_typ"]
            row["delta_faint"] = lim - row["m_pred_faint"]
    panel = None
    if "_sub" in m:
        panel = {"img": m["_sub"], "res": m["_res"], "clean": m["_clean"]}
    return row, panel, inj


def _write(rows, name, prov, source):
    t = Table(rows=rows)
    t.meta.update(provenance=prov, source=source)
    t.write(RES / name, overwrite=True)
    return t


def cmd_measure(args) -> None:
    p = ch.Params()
    sel = Table.read(RES / "selection.ecsv")
    sel_all = sel
    if args.role == "target":
        gate = json.loads((RES / "gate.json").read_text())
        if not gate["passed"]:
            raise SystemExit("control gate failed: targets are not measured (PLAN.md 4.1)")
    sel = sel[np.asarray(sel["role"]) == args.role]
    if args.only:
        sel = sel[np.isin(np.asarray(sel["system_id"]), args.only.split(","))]
    rows, inj = [], []
    panels_dir = Path(args.cut).parent / "panels"  # untracked stamps for the contact sheet
    panels_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.process_time()
    for r in sel:
        row, panel, ij = run_system(r, Path(args.cut), p, do_inject=not args.no_inject)
        rows.append(row)
        inj += ij
        if panel is not None:
            np.savez_compressed(panels_dir / f"{row['system_id']}.npz", **panel)
        keys = ("name", "status", "n_img", "max_sep_arcsec", "resid_snr", "m_lim_reff0.3")
        print({k: row.get(k) for k in keys}, flush=True)
    cpu = time.process_time() - t0
    src = "scripts/w2_hst_blend_run.py measure on MAST HAP-MVM F814W cutouts (hapcut)"
    role_ids = set(np.asarray(sel_all[np.asarray(sel_all["role"]) == args.role]["system_id"]))
    t = merge_rows(RES / "measurements.ecsv", rows, args.role, role_ids)
    t.meta.update(provenance=schema.Provenance.DERIVED.value, source=src, params=asdict(p))
    t.write(RES / "measurements.ecsv", overwrite=True)
    if inj or (RES / "injections.ecsv").exists():
        # rows of re-measured systems are replaced, or dropped when run with --no-inject
        ti = merge_rows(
            RES / "injections.ecsv", inj, args.role, role_ids, {r["system_id"] for r in rows}
        )
        ti.meta.update(
            provenance=schema.Provenance.SIMULATED.value,
            source="Sersic n=4, ellip 0.2, 5x5 subsampled, convolved with the fitted Moffat PSF, "
            "injected into the real cutout and run through the whole chain (w2_hst_blend_run.py)",
        )
        ti.write(RES / "injections.ecsv", overwrite=True)
    (RES / "params.json").write_text(json.dumps(asdict(p), indent=2) + "\n")
    if args.role == "control":
        full = t[np.asarray(t["role"]) == "control"]
        full_chain = chain_status(full)
        g = ch.control_gate(list(full_chain), p)
        g["provenance"] = "derived; pre-registered control gate (PLAN.md 4.1)"
        g["cpu_s_measure_controls"] = round(cpu, 1)
        g["statuses"] = dict(zip(np.asarray(full["name"]).astype(str), full_chain, strict=True))
        (RES / "gate.json").write_text(json.dumps(g, indent=2) + "\n")
        print(json.dumps(g, indent=1))
    else:
        (RES / "cpu_targets.json").write_text(json.dumps({"cpu_s": round(cpu, 1)}) + "\n")


def merge_rows(path: Path, rows, role: str, role_ids, measured=None) -> Table:
    """New rows for ``role`` replace only the systems just measured; rows of the same role for
    other systems of the current selection are kept (``measure --only``); rows of systems no
    longer selected for that role are dropped; other roles are kept."""
    measured = {r["system_id"] for r in rows} if measured is None else set(measured)
    new = Table(rows=rows) if rows else None
    if not path.exists():
        return new if new is not None else Table(names=("system_id", "role"), dtype=(str, str))
    prev = Table.read(path)
    sid = np.asarray(prev["system_id"]).astype(str)
    rl = np.asarray(prev["role"]).astype(str)
    keep = (rl != role) | (np.isin(sid, list(role_ids)) & ~np.isin(sid, list(measured)))
    prev = prev[keep]
    if new is None:
        return prev
    return vstack([prev, new], metadata_conflicts="silent") if len(prev) else new


PSF_MAX_SYSTEMS = 3  # PLAN.md: up to 3 controls with an empirical PSF
PSF_MAX_TRIES = 5  # wide cutouts tried at most (download budget)
WIDE_PX = 2200  # w12_niq_epsf.Params.size_px


def simulate(shape, pos, fluxes, psf_img, lens_unit, lens_flux, sigma, seed):
    """Noise + N point images made with ``psf_img`` (unit-sum stamp) + an optional lens."""
    import w12_niq_epsf as ep

    rng = np.random.default_rng(seed)
    im = rng.normal(0.0, sigma, shape)
    for (x, y), f in zip(pos, fluxes, strict=True):
        im += f * ep.place(psf_img, shape, x, y)
    if lens_unit is not None:
        im += lens_flux * lens_unit
    return im


def gaia_with_retry(c: SkyCoord, radius_deg: float, tries: int = 3):
    import w12_niq_epsf as ep

    msg = ""
    for attempt in range(tries):
        try:
            return ep.gaia_stars(c, radius_deg)
        except Exception as e:  # noqa: BLE001 - retried, then raised (never "no stars")
            msg = f"{type(e).__name__}: {e}"
            time.sleep(10 * (attempt + 1))
    raise RuntimeError(f"Gaia query failed: {msg[:300]}")


def psf_one(sid, r, m, d: Path, log: Path, p: ch.Params, epp) -> tuple[list, bool]:
    """One control: wide cutout, Gaia-star ePSF, simulated stamps. Returns (rows, has_epsf)."""

    import w12_niq_epsf as ep

    c = SkyCoord(float(r["ra"]), float(r["dec"]), unit="deg")
    sf, key = RES / DOWNLOAD_STATUS, f"wide_{sid}"
    status = json.loads(sf.read_text()) if sf.exists() else {}
    if needs_download(status, key, d):
        _, st, msg = download_set(sf, status, key, str(r["name"]), c, WIDE_PX, d, log)
        if st == FAILED:
            raise RuntimeError(f"{r['name']}: wide hapcut download failed: {msg}")
    cu = load_cutout(d, c, p) if cutout_complete(status, key, d) else None
    base_row = {"system_id": sid, "name": str(r["name"]), "n_psf_stars": 0, "psf": "",
                "reff_arcsec": np.nan, "inject_mag": np.nan, "n_img": 0,
                "resid_snr": np.nan, "recovered": False}  # fmt: skip
    if cu is None:
        return [base_row | {"status": "no wide F814W cutout"}], False
    g = gaia_with_retry(c, WIDE_PX * cu["pix"] / 3600 / 2)
    gc = SkyCoord(g["ra"], g["dec"], unit="deg")
    _, d2, _ = gc.match_to_catalog_sky(gc, nthneighbor=2)
    gm = np.asarray(g["phot_g_mean_mag"], float)
    ok = (gm > epp.star_gmag[0]) & (gm < epp.star_gmag[1])
    ok &= gc.separation(c).arcsec > epp.star_min_dist_arcsec
    ok &= d2.arcsec > epp.star_isolation_arcsec
    w = WCS(fits.getheader(sorted(d.glob(cu["file"]))[0], 1))
    xy = [tuple(map(float, w.world_to_pixel(s))) for s in gc[ok]]
    psf, nstar = ep.build_epsf(cu["img"], cu["valid"], xy, epp.psf_half, epp.flat_core)
    print(r["name"], "PSF stars", nstar, flush=True)
    for f in d.glob("*.fits"):  # cloud disk rule: raw wide cutouts deleted after use
        f.unlink()
    if psf is None or nstar < epp.min_psf_stars:
        return [base_row | {"n_psf_stars": nstar, "status": "too few PSF stars"}], False
    pos = [tuple(map(float, q.split(","))) for q in str(m["img_pos_px"]).split(";")]
    amps = [float(a) for a in str(m["img_amps"]).split(";")]
    gam, alp = float(m["psf_gamma"]), float(m["psf_alpha"])
    fluxes = [moffat_flux(a, gam, alp) for a in amps]
    shape = (p.size_px, p.size_px)
    mk = ch.moffat_kernel(gam, alp, half=epp.psf_half)
    lx, ly = ch.lens_position(pos, amps)
    centre = (float(np.mean([q[0] for q in pos])), float(np.mean([q[1] for q in pos])))
    sig = float(m["sigma_pix"])
    out = []
    for kind, kern in (("epsf", psf), ("moffat", mk)):
        base_im = simulate(shape, pos, fluxes, kern, None, 0.0, sig, 0)
        valid = np.ones(shape, bool)
        _, _, nfn = noise_function(base_im, valid)
        base = ch.measure(base_im, valid, centre, cu["pix"], p, sig, nfn)
        row0 = base_row | {"n_psf_stars": nstar, "psf": kind}
        out.append(row0 | {"status": base["status"], "n_img": base["n_img"],
                           "resid_snr": base["resid_snr"]})  # fmt: skip
        print(r["name"], kind, "no lens:", base["status"], base["resid_snr"], flush=True)
        for reff in p.inj_reff_arcsec:
            lu = ch.sersic_unit(shape, lx, ly, reff / cu["pix"], p.inj_n, p.inj_ellip, psf=kern)
            for mag in p.inj_mags:
                fl = 10 ** (-0.4 * (mag - float(m["zp_ab"])))
                im = simulate(shape, pos, fluxes, kern, lu, fl, sig, 0)
                mi = ch.measure(im, valid, centre, cu["pix"], p, sig, nfn)
                out.append(row0 | {"reff_arcsec": reff, "inject_mag": mag,
                                   "status": mi["status"], "n_img": mi["n_img"],
                                   "resid_snr": mi["resid_snr"],
                                   "recovered": ch.recovered(base, mi, p)})  # fmt: skip
    return out, True


def finite_or_none(x):
    """JSON-safe number: NaN or inf become null (json.dumps would write invalid NaN)."""
    return float(x) if x is not None and np.isfinite(x) else None


def no_lens_max(t: Table) -> tuple[float | None, str | None]:
    """Largest no-lens empirical-PSF S/N and its system, or (None, None) when there is no such
    row or every value is NaN (not measurable; ``np.nanargmax`` would raise)."""
    if not len(t) or "psf" not in t.colnames:
        return None, None
    nolens = (np.asarray(t["psf"]).astype(str) == "epsf") & ~np.isfinite(fcol(t, "inject_mag"))
    v = fcol(t, "resid_snr")[nolens]
    ids = np.asarray(t["system_id"]).astype(str)[nolens]
    fin = np.isfinite(v)  # NaN and +-inf are not measurable values
    if not fin.any():
        return None, None
    i = int(np.argmax(v[fin]))
    return float(v[fin][i]), str(ids[fin][i])


V1_LABELS = {"none between images": "not detected between images"}  # v1 -> v2 wording


def normalise_status(st) -> str:
    return V1_LABELS.get(str(st), str(st))


REUSE_KEYS = ("hst_image", "img_pos_px", "img_amps", "psf_gamma", "psf_alpha", "sigma_pix", "zp_ab")


def reuse_ok(old_row, cur_row) -> bool:
    """Rows of an earlier PSF run may be reused only if the control measurement they were built
    from has the same inputs as now: cutout identity, image positions and amplitudes, PSF shape,
    sigma_pix and zero point (string fields compared exactly, numbers to 1e-9 relative)."""
    for k in REUSE_KEYS:
        try:
            a, b = old_row[k], cur_row[k]
        except (KeyError, ValueError):  # a missing input column is a mismatch: recompute
            return False
        if isinstance(a, str | np.str_) or isinstance(b, str | np.str_):
            if str(a) != str(b):
                return False
        elif not np.isclose(float(a), float(b), rtol=1e-9, atol=0.0, equal_nan=True):
            return False
    return True


def write_psf(out, done, tries, p: ch.Params, reused, refused=()) -> None:
    t = Table(rows=out) if out else Table(names=("system_id", "psf"), dtype=(str, str))
    t.meta.update(
        provenance=schema.Provenance.SIMULATED.value,
        source="w2_hst_blend_run.py psf: white noise at the control's sigma_pix, its fitted "
        "images made with a Gaia-star empirical PSF (w12_niq_epsf.build_epsf) or the fitted "
        "Moffat, optional Sersic lens convolved with the same PSF, measured by the Moffat chain",
    )
    t.write(RES / "psf_mismatch.ecsv", overwrite=True)
    fp, fp_sys = no_lens_max(t)
    psfc = np.asarray(t["psf"]).astype(str) if len(t) else np.zeros(0, str)
    summ = {
        "systems_with_epsf": done,
        "wide_cutouts_tried": tries,
        "rows_reused": reused,
        "reuse_refused_input_mismatch": list(refused),
        "epsf_no_lens_max_snr": fp,  # null = no negative control or not measurable
        "epsf_no_lens_max_snr_system": fp_sys,
        "psf_mismatch_not_excluded": bool(fp is not None and fp >= p.detect_sigma),
        "provenance": "simulated; pre-registered PSF-mismatch test with negative controls "
        "(no-lens empirical-PSF stamps), PLAN.md 4.2",
    }
    for sid in done:
        for kind in ("epsf", "moffat"):
            for reff in p.inj_reff_arcsec:
                sel = (np.asarray(t["system_id"]).astype(str) == sid) & (psfc == kind)
                s = t[sel & (fcol(t, "reff_arcsec") == reff)]
                summ[f"m_lim_{sid}_{kind}_reff{reff}"] = finite_or_none(
                    ch.contiguous_limit([(float(x["inject_mag"]), bool(x["recovered"])) for x in s])
                )
    (RES / "psf_mismatch.json").write_text(json.dumps(summ, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summ, indent=1))


def chain_status(meas: Table) -> np.ndarray:
    """The chain's own status per row (``chain_output_uninterpreted`` where ``verdict`` has
    relabelled targets; ``status`` otherwise)."""
    st = np.asarray(meas["status"]).astype(str)
    if "chain_output_uninterpreted" in meas.colnames:
        c = meas["chain_output_uninterpreted"]
        raw = np.asarray(c.filled("") if hasattr(c, "filled") else c).astype(str)
        st = np.where(np.isin(raw, ("", "--", "None")), st, raw)
    return st


def pick_reuse(sid: str, cur_meas_row, sources) -> tuple[list | None, bool]:
    """First source whose rows for ``sid`` were built from a measurement equal to the current one
    (``reuse_ok``). Returns (rows with v2 labels, or None; whether some source had rows for
    ``sid`` but a mismatching measurement). A refusal counts only when no source is accepted."""
    mismatched = False
    for prev, pmeas in sources:
        if sid not in set(np.asarray(prev["system_id"]).astype(str)):
            continue
        old = pmeas[np.asarray(pmeas["system_id"]).astype(str) == sid]
        if len(old) and reuse_ok(old[0], cur_meas_row):
            rows = [dict(zip(prev.colnames, x, strict=True)) for x in prev
                    if str(x["system_id"]) == sid]  # fmt: skip
            for x in rows:
                x["status"] = normalise_status(x["status"])
            return rows, False
        mismatched = True
    return None, mismatched


def cmd_psf(args) -> None:
    """PSF-mismatch simulation (PLAN.md 4.2): images made with an empirical Gaia-star PSF, fitted
    by the Moffat chain; the same with Moffat-made images as the no-mismatch reference.

    ``--reuse PSF.ecsv=MEAS.ecsv`` (repeatable): rows of a system from an earlier run are reused
    only if the control measurement they were built from (in MEAS.ecsv) has the same inputs as
    the current one (``reuse_ok``); their status labels are normalised to the v2 wording. A
    refused or absent system is recomputed (wide cutout). Partial results are written before an
    error is re-raised; a failure while writing them never masks the original error."""
    import w12_niq_epsf as ep

    p = ch.Params()
    epp = ep.Params()
    meas = Table.read(RES / "measurements.ecsv")
    ctl = meas[(np.asarray(meas["role"]) == "control") & np.isin(chain_status(meas), ch.DECIDED)]
    sel = Table.read(RES / "selection.ecsv")
    sel = sel[np.asarray(sel["role"]) == "control"]
    order = [s for s in sel["system_id"] if s in set(ctl["system_id"])]
    sources = []
    for spec in args.reuse or []:
        pf, mf = spec.split("=")
        sources.append((Table.read(pf), Table.read(mf)))
    out, done, tries, reused, refused = [], [], 0, [], []
    log = Path(args.log)
    try:
        for sid in order:
            if len(done) >= PSF_MAX_SYSTEMS or tries >= PSF_MAX_TRIES:
                break
            tries += 1
            m = ctl[ctl["system_id"] == sid][0]
            rows, mismatched = pick_reuse(sid, m, sources)
            if rows is None and mismatched and sid not in refused:
                refused.append(sid)  # no source accepted: recomputed below
            if rows is not None:
                reused.append(sid)
                out += rows
                if any(str(x["psf"]) == "epsf" for x in rows):
                    done.append(sid)
                continue
            r = sel[sel["system_id"] == sid][0]
            rows, ok = psf_one(sid, r, m, Path(args.cut) / f"wide_{sid}", log, p, epp)
            out += rows
            if ok:
                done.append(sid)
    except BaseException:
        try:
            write_psf(out, done, tries, p, reused, refused)
        except Exception as e2:  # noqa: BLE001 - never mask the original error
            print("writing partial PSF results failed:", type(e2).__name__, e2, flush=True)
        raise
    write_psf(out, done, tries, p, reused, refused)


def cmd_sheet(args) -> None:
    """Contact sheet: F814W stamp, N-Moffat residual, halo-corrected residual (smoothed 1 px)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image
    from scipy import ndimage

    meas = Table.read(RES / "measurements.ecsv")
    tgt_label = verdict_label()  # "" when the verdict passed (targets keep the chain status)
    roles = args.roles.split(",")
    meas = meas[np.isin(np.asarray(meas["role"]), roles)]
    pd = Path(args.cut).parent / "panels"
    chain = dict(zip(np.asarray(meas["system_id"]).astype(str), chain_status(meas), strict=True))
    items = [(r, pd / f"{r['system_id']}.npz") for r in meas]
    items = [(r, f) for r, f in items if f.exists()]
    fig, axes = plt.subplots(len(items), 3, figsize=(7.5, 2.5 * len(items)), squeeze=False)
    for i, (r, f) in enumerate(items):
        z = np.load(f)
        img, res, clean = z["img"], z["res"], z["clean"]
        vmax = np.percentile(img, 99.7)
        rs = np.percentile(np.abs(res), 99)
        sm = ndimage.gaussian_filter(clean, 1.0)
        ss = max(np.percentile(np.abs(sm), 99.5), 1e-12)
        panels = (
            (img, f"{r['name']} [{r['role']}] F814W", {"cmap": "gray_r", "vmin": -0.02 * vmax,
                                                       "vmax": vmax}),
            (res, f"N={r['n_img']} Moffat resid (raw S/N {r['resid_snr_raw']:.1f})",
             {"cmap": "RdBu_r", "vmin": -rs, "vmax": rs}),
            (sm, (f"stored S/N {r['resid_snr']:.1f} ({tgt_label})"
                  if r["role"] == "target" and tgt_label
                  else f"{chain[str(r['system_id'])]} S/N {r['resid_snr']:.1f}"),
             {"cmap": "RdBu_r", "vmin": -ss, "vmax": ss}),
        )  # fmt: skip
        for ax, (im, title, kw) in zip(axes[i], panels, strict=True):
            ax.imshow(im, origin="lower", **kw)
            ax.set_title(title, fontsize=7)
            ax.set_xticks([])
            ax.set_yticks([])
    fig.tight_layout()
    out = RES / args.out
    png = out.with_suffix(".png")
    fig.savefig(png, dpi=60)
    Image.open(png).convert("RGB").save(out, quality=80)
    png.unlink()
    print(out, out.stat().st_size)


CENTRE_AP_ARCSEC = 0.2  # post hoc diagnostic aperture at the images' centroid


def fp_reference(pmt: Table) -> str | None:
    """system_id of the ePSF simulation with the largest no-lens S/N (None if not measurable)."""
    return no_lens_max(pmt)[1]


def verdict_label() -> str:
    """Short target label from ``verdict.json``: "" if the method was validated, else the
    "not interpreted" wording without the reason in brackets (as on the sheets)."""
    vf = RES / "verdict.json"
    if not vf.exists():
        return "verdict not run"
    v = json.loads(vf.read_text())
    return "" if v.get("validated") else str(v["target_interpretation"]).split(" (")[0]


def cmd_posthoc(args) -> None:
    """Post hoc diagnostics (not frozen rules; RESULT.md labels them): (1) the PSF-mismatch
    false-positive S/N of the one ePSF simulation scaled to each system by image flux and
    aperture error (ASSUMPTION: the same fractional mismatch); (2) for each "none between
    images" system, the halo-corrected residual in a 0.2" circle at the images' centroid."""
    p = ch.Params()
    meas = Table.read(RES / "measurements.ecsv")
    pm = json.loads((RES / "psf_mismatch.json").read_text())
    # one reference for all three quantities: the ePSF system with the largest no-lens S/N
    ref_id = pm.get("epsf_no_lens_max_snr_system") or fp_reference(
        Table.read(RES / "psf_mismatch.ecsv")
    )
    if ref_id is None or pm.get("epsf_no_lens_max_snr") is None:
        raise SystemExit("no measurable negative control: the post hoc scaling is undefined")
    ref = meas[meas["system_id"] == ref_id][0]
    # controls only, on purpose (owner decision 2026-10-21): when the method is not validated,
    # no target is interpreted, so the post hoc diagnostics are computed for controls alone.
    # The filter uses the chain's own status (chain_status), never the verdict relabel.
    vf = RES / "verdict.json"
    v = json.loads(vf.read_text()) if vf.exists() else {}
    chain = chain_status(meas)
    use = np.asarray(meas["role"]).astype(str) == "control"
    if v.get("validated"):
        use |= np.asarray(meas["role"]).astype(str) == "target"

    def img_flux(r):
        return sum(
            10 ** (-0.4 * (float(x) - float(r["zp_ab"]))) for x in r["image_mags"].split(";")
        )

    out = {"reference": str(ref["name"]), "fp_snr_ref": finite_or_none(pm["epsf_no_lens_max_snr"]),
           "systems": {}}  # fmt: skip
    for r, cs in zip(
        meas[use & np.isin(chain, ch.DECIDED)], chain[use & np.isin(chain, ch.DECIDED)], strict=True
    ):
        fp = pm["epsf_no_lens_max_snr"] * img_flux(r) / img_flux(ref)
        fp *= float(ref["resid_err"]) / float(r["resid_err"])
        out["systems"][str(r["name"])] = {
            "role": str(r["role"]),
            "status": str(cs),
            "resid_snr": finite_or_none(float(r["resid_snr"])),
            "fp_snr_scaled": finite_or_none(float(fp)),
            "snr_over_fp_scaled": finite_or_none(r["resid_snr"] / fp),
            "psf_alpha_at_bound": bool(float(r["psf_alpha"]) <= 1.2 + 1e-6),
        }
    sel = Table.read(RES / "selection.ecsv")
    for r in meas[use & (chain == ch.DECIDED[1])]:
        s = sel[sel["system_id"] == r["system_id"]][0]
        c = SkyCoord(float(s["ra"]), float(s["dec"]), unit="deg")
        cu = load_cutout(Path(args.cut) / str(r["system_id"]), c, p)
        sig, _corr, nfn = noise_function(cu["data"], cu["valid"])
        m = ch.measure(cu["img"], cu["valid"], cu["centre"], cu["pix"], p, sig, nfn)
        y0, _y1, x0, _x1 = m["stamp"]
        xc = np.mean([q[0] for q in m["pos"]]) - x0
        yc = np.mean([q[1] for q in m["pos"]]) - y0
        yy, xx = np.indices(m["_clean"].shape)
        ap = np.hypot(xx - xc, yy - yc) <= CENTRE_AP_ARCSEC / cu["pix"]
        npix = int(ap.sum())
        flux = float(m["_clean"][ap].sum())
        err = nfn(npix) * np.sqrt(npix)
        out["systems"][str(r["name"])].update(
            centre_ap_arcsec=CENTRE_AP_ARCSEC,
            centre_flux=finite_or_none(flux),
            centre_snr=finite_or_none(flux / err),
            centre_mag=finite_or_none(cu["zp"] - 2.5 * np.log10(flux)) if flux > 0 else None,
        )
    scope = (
        "controls and targets (method validated)"
        if v.get("validated")
        else f"controls only, on purpose: targets {v.get('target_interpretation', 'not judged')} "
        f"(verdict.json: {v.get('verdict', 'verdict not run')})"
    )
    out["provenance"] = (
        "derived; POST HOC diagnostics, not pre-registered in PLAN.md; not used to interpret "
        f"targets; scope: {scope}"
    )
    (RES / "posthoc.json").write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    print(json.dumps(out, indent=1))


def fcol(t: Table, name: str) -> np.ndarray:
    """Float column with masked entries as NaN (a masked value is not 0)."""
    c = t[name]
    return np.asarray(c.filled(np.nan) if hasattr(c, "filled") else c, float)


def role_summary(meas: Table, role: str) -> dict:
    """Status counts, injection limits and separation range for one role; an empty role (e.g.
    targets when the gate failed) gives empty counts and None, never a crash."""
    m = meas[np.asarray(meas["role"]).astype(str) == role] if len(meas) else meas
    st = np.asarray(m["status"]).astype(str) if len(m) else np.zeros(0, str)
    lim = fcol(m, "m_lim_reff0.3") if len(m) and "m_lim_reff0.3" in m.colnames else np.zeros(0)
    dsep = fcol(m, "max_sep_arcsec") if len(m) and "max_sep_arcsec" in m.colnames else np.zeros(0)
    ok = np.isfinite(lim)
    fin = dsep[np.isfinite(dsep)]
    return {
        f"{role}_status_counts": {k: int(np.sum(st == k)) for k in sorted(set(st))},
        f"{role}_m_lim_reff0.3": {
            "n": int(ok.sum()),
            "min": float(np.min(lim[ok])) if ok.any() else None,
            "median": float(np.median(lim[ok])) if ok.any() else None,
            "max": float(np.max(lim[ok])) if ok.any() else None,
        },
        f"{role}_max_sep_range": [float(fin.min()), float(fin.max())] if len(fin) else None,
    }


def cmd_summary(args) -> None:
    """Numbers quoted in RESULT.md, from the tracked tables (``derived``)."""
    meas = Table.read(RES / "measurements.ecsv")
    fi, fp = RES / "injections.ecsv", RES / "psf_mismatch.json"
    inj = (
        Table.read(fi) if fi.exists() else Table(names=("role", "reff_arcsec"), dtype=(str, float))
    )
    pm = json.loads(fp.read_text()) if fp.exists() else None
    gate = json.loads((RES / "gate.json").read_text())
    out = {"gate": {k: gate[k] for k in ("n_controls", "n_decided", "n_light", "efficiency",
                                         "passed")}}  # fmt: skip
    from scipy.stats import beta

    k, n = gate["n_light"], gate["n_decided"]
    # one-sided 95 % Clopper-Pearson lower bound on the control efficiency
    out["gate"]["efficiency_lower95_one_sided"] = float(beta.ppf(0.05, k, n - k + 1)) if k else 0.0
    out["psf_mismatch"] = pm
    rec = {}
    for role in ("control", "target"):
        for reff in sorted(set(fcol(inj, "reff_arcsec"))):
            s = inj[
                (np.asarray(inj["role"]).astype(str) == role) & (fcol(inj, "reff_arcsec") == reff)
            ]
            rec[f"{role}_reff{reff}"] = {
                f"{mag:g}": [int(np.sum(q["recovered"])), len(q)]
                for mag in sorted(set(np.asarray(s["inject_mag"], float)))
                for q in [s[np.asarray(s["inject_mag"], float) == mag]]
            }
            ff = np.asarray(s["flux_frac"], float)
            rec[f"{role}_reff{reff}_flux_frac_median"] = float(np.nanmedian(ff)) if len(s) else None
            rec[f"{role}_reff{reff}_flux_frac_p10_p90"] = (
                [float(x) for x in np.nanpercentile(ff, [10, 90])] if len(s) else None
            )
    out["injections_recovered_of_n"] = rec
    for role in ("control", "target"):
        out.update(role_summary(meas, role))
    a = fcol(meas, "psf_alpha")
    out["psf_alpha_at_lower_bound"] = int(np.sum(np.isfinite(a) & (a <= 1.2 + 1e-6)))
    out["psf_fits"] = int(np.sum(np.isfinite(a)))
    if args.downloads_log:
        lines = [x.split("\t") for x in Path(args.downloads_log).read_text().splitlines()]
        lines = [x for x in lines if len(x) >= 4]  # note lines (e.g. deletions) are skipped
        b = [int(x[2]) for x in lines]
        out["downloads"] = {"files": len(b), "bytes": int(sum(b))}
        for kind, key in (("hapcut 250 px", "size=250"), ("hapcut 2200 px", "size=2200")):
            sel = [int(x[2]) for x in lines if key in x[1]]
            out["downloads"][kind] = {"files": len(sel), "bytes": int(sum(sel))}
    vf = RES / "verdict.json"
    out["verdict"] = json.loads(vf.read_text()) if vf.exists() else None
    out["provenance"] = (
        "derived; gate, injections and psf_mismatch are pre-registered (PLAN.md); verdict is the "
        "owner decision of 2026-10-21 applied after the plan was frozen; nothing here is post hoc "
        "(post hoc diagnostics are only in posthoc.json)"
    )
    (RES / "summary.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=1))


RAW = "chain_output_uninterpreted"  # the chain's status string, stored, never a reading


def apply_verdict(meas: Table, v: dict) -> Table:
    """Idempotent: ``status`` is always re-derived from the chain's own status (``RAW``, created
    or filled from ``status`` for rows that lack it, e.g. rows added later by ``measure``); only
    when the verdict fails do targets get the "not interpreted" status. A passing verdict
    restores the chain status."""
    meas = Table(meas, masked=False) if not len(meas) else meas.copy()
    if "tool_status" in meas.colnames:  # name used before d9dab07's final form
        meas.rename_column("tool_status", RAW)
    chain = [str(x) for x in chain_status(meas)]
    tgt = np.asarray(meas["role"]).astype(str) == "target"
    st = [
        v["target_interpretation"] if (t and not v["validated"]) else c
        for c, t in zip(chain, tgt, strict=True)
    ]
    meas[RAW] = np.array(chain, dtype=str)
    meas["status"] = np.array(st, dtype=str)
    meas.meta["verdict"] = v
    return meas


def cmd_verdict(args) -> None:
    """Owner decision 2026-10-21 (after the frozen plan): method validation needs the gate AND
    every negative control (no-lens empirical-PSF stamp) below the detection threshold."""
    p = ch.Params()
    gate = json.loads((RES / "gate.json").read_text())
    fp = RES / "psf_mismatch.json"
    v = ch.method_validation(gate, json.loads(fp.read_text()) if fp.exists() else None, p)
    (RES / "verdict.json").write_text(json.dumps(v, indent=2, allow_nan=False) + "\n")
    meas = apply_verdict(Table.read(RES / "measurements.ecsv"), v)
    meas.write(RES / "measurements.ecsv", overwrite=True)
    print(json.dumps(v, indent=1))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("download")
    a.add_argument("--cut", required=True)
    a.add_argument("--log", required=True)
    a.add_argument("--force", action="store_true", help="re-query recorded no-data systems")
    a = sub.add_parser("measure")
    a.add_argument("--role", choices=("control", "target"), required=True)
    a.add_argument("--cut", required=True)
    a.add_argument("--only", default="")
    a.add_argument("--no-inject", action="store_true")
    a = sub.add_parser("psf")
    a.add_argument("--cut", required=True)
    a.add_argument("--log", required=True)
    a.add_argument(
        "--reuse",
        action="append",
        default=[],
        help="PSF.ecsv=MEAS.ecsv: earlier PSF rows and the measurement table they were built from",
    )
    a = sub.add_parser("sheet")
    a.add_argument("--cut", required=True)
    a.add_argument("--roles", default="control,target")
    a.add_argument("--out", default="contact_sheet.jpg")
    a = sub.add_parser("posthoc")
    a.add_argument("--cut", required=True)
    sub.add_parser("verdict")
    a = sub.add_parser("summary")
    a.add_argument("--downloads-log", default="")
    args = ap.parse_args(argv)
    cmds = {"download": cmd_download, "measure": cmd_measure, "psf": cmd_psf, "sheet": cmd_sheet}
    cmds["posthoc"] = cmd_posthoc
    cmds["summary"] = cmd_summary
    cmds["verdict"] = cmd_verdict
    cmds[args.cmd](args)


if __name__ == "__main__":
    main()
