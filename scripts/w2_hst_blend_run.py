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


def cmd_download(args) -> None:
    from astroquery.mast import Hapcut

    p = ch.Params()
    sel = Table.read(RES / "selection.ecsv")
    log = Path(args.log)
    for r in sel:
        d = Path(args.cut) / str(r["system_id"])
        if d.exists() and list(d.glob("*.fits")):
            continue
        d.mkdir(parents=True, exist_ok=True)
        c = SkyCoord(float(r["ra"]), float(r["dec"]), unit="deg")
        for attempt in range(3):  # a service error is retried, never recorded as "no data"
            try:
                t = Hapcut.download_cutouts(c, size=[p.size_px, p.size_px], path=str(d))
                break
            except Exception as e:  # noqa: BLE001
                print(r["name"], "hapcut error", type(e).__name__, e, flush=True)
                time.sleep(10 * (attempt + 1))
        else:
            continue
        for f in t["Local Path"]:
            log_download(log, f"hapcut ra={c.ra.deg:.6f} dec={c.dec.deg:.6f} size={p.size_px}", f)
        print(r["name"], len(t), "files", flush=True)


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
    cu = load_cutout(cut / str(r["system_id"]), c, p)
    inj = []
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
    old = RES / "measurements.ecsv"
    t = Table(rows=rows)
    if old.exists():
        prev = Table.read(old)
        prev = prev[np.asarray(prev["role"]) != args.role]
        t = vstack([prev, t], metadata_conflicts="silent") if len(prev) else t
    t.meta.update(provenance=schema.Provenance.DERIVED.value, source=src, params=asdict(p))
    t.write(old, overwrite=True)
    if inj:
        ti = Table(rows=inj)
        oi = RES / "injections.ecsv"
        if oi.exists():
            prev = Table.read(oi)
            prev = prev[np.asarray(prev["role"]) != args.role]
            ti = vstack([prev, ti], metadata_conflicts="silent") if len(prev) else ti
        ti.meta.update(
            provenance=schema.Provenance.SIMULATED.value,
            source="Sersic n=4, ellip 0.2, 5x5 subsampled, convolved with the fitted Moffat PSF, "
            "injected into the real cutout and run through the whole chain (w2_hst_blend_run.py)",
        )
        ti.write(oi, overwrite=True)
    (RES / "params.json").write_text(json.dumps(asdict(p), indent=2) + "\n")
    if args.role == "control":
        g = ch.control_gate([row["status"] for row in rows], p)
        g["cpu_s_measure_controls"] = round(cpu, 1)
        g["statuses"] = {row["name"]: row["status"] for row in rows}
        (RES / "gate.json").write_text(json.dumps(g, indent=2) + "\n")
        print(json.dumps(g, indent=1))
    else:
        (RES / "cpu_targets.json").write_text(json.dumps({"cpu_s": round(cpu, 1)}) + "\n")


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


def cmd_psf(args) -> None:
    """PSF-mismatch simulation (PLAN.md 4.2): images made with an empirical Gaia-star PSF, fitted
    by the Moffat chain; the same with Moffat-made images as the no-mismatch reference."""
    import w12_niq_epsf as ep
    from astroquery.mast import Hapcut

    p = ch.Params()
    epp = ep.Params()
    meas = Table.read(RES / "measurements.ecsv")
    ctl = meas[(np.asarray(meas["role"]) == "control") & np.isin(meas["status"], ch.DECIDED)]
    sel = Table.read(RES / "selection.ecsv")
    order = [s for s in sel["system_id"] if s in set(ctl["system_id"])]
    out, done, tries = [], [], 0
    log = Path(args.log)
    for sid in order:
        if len(done) >= PSF_MAX_SYSTEMS or tries >= PSF_MAX_TRIES:
            break
        tries += 1
        r = sel[sel["system_id"] == sid][0]
        m = ctl[ctl["system_id"] == sid][0]
        c = SkyCoord(float(r["ra"]), float(r["dec"]), unit="deg")
        d = Path(args.cut) / f"wide_{sid}"
        if not (d.exists() and list(d.glob("*.fits"))):
            d.mkdir(parents=True, exist_ok=True)
            t = Hapcut.download_cutouts(c, size=[WIDE_PX, WIDE_PX], path=str(d))
            for f in t["Local Path"]:
                log_download(log, f"hapcut ra={c.ra.deg:.6f} dec={c.dec.deg:.6f} size={WIDE_PX}", f)
        cu = load_cutout(d, c, p)
        if cu is None:
            continue
        g = ep.gaia_stars(c, WIDE_PX * cu["pix"] / 3600 / 2)
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
            out.append({"system_id": sid, "name": str(r["name"]), "n_psf_stars": nstar,
                        "psf": "", "reff_arcsec": np.nan, "inject_mag": np.nan,
                        "status": "too few PSF stars", "n_img": 0, "resid_snr": np.nan,
                        "recovered": False})  # fmt: skip
            continue
        done.append(sid)
        pos = [tuple(map(float, q.split(","))) for q in str(m["img_pos_px"]).split(";")]
        amps = [float(a) for a in str(m["img_amps"]).split(";")]
        gam, alp = float(m["psf_gamma"]), float(m["psf_alpha"])
        fluxes = [moffat_flux(a, gam, alp) for a in amps]
        shape = (p.size_px, p.size_px)
        mk = ch.moffat_kernel(gam, alp, half=epp.psf_half)
        lx, ly = ch.lens_position(pos, amps)
        centre = (float(np.mean([q[0] for q in pos])), float(np.mean([q[1] for q in pos])))
        sig = float(m["sigma_pix"])
        for kind, kern in (("epsf", psf), ("moffat", mk)):
            base_im = simulate(shape, pos, fluxes, kern, None, 0.0, sig, 0)
            valid = np.ones(shape, bool)
            _, _, nfn = noise_function(base_im, valid)
            base = ch.measure(base_im, valid, centre, cu["pix"], p, sig, nfn)
            out.append({"system_id": sid, "name": str(r["name"]), "n_psf_stars": nstar,
                        "psf": kind, "reff_arcsec": np.nan, "inject_mag": np.nan,
                        "status": base["status"], "n_img": base["n_img"],
                        "resid_snr": base["resid_snr"], "recovered": False})  # fmt: skip
            print(r["name"], kind, "no lens:", base["status"], base["resid_snr"], flush=True)
            for reff in p.inj_reff_arcsec:
                lu = ch.sersic_unit(shape, lx, ly, reff / cu["pix"], p.inj_n, p.inj_ellip, psf=kern)
                for mag in p.inj_mags:
                    fl = 10 ** (-0.4 * (mag - float(m["zp_ab"])))
                    im = simulate(shape, pos, fluxes, kern, lu, fl, sig, 0)
                    mi = ch.measure(im, valid, centre, cu["pix"], p, sig, nfn)
                    out.append({"system_id": sid, "name": str(r["name"]), "n_psf_stars": nstar,
                                "psf": kind, "reff_arcsec": reff, "inject_mag": mag,
                                "status": mi["status"], "n_img": mi["n_img"],
                                "resid_snr": mi["resid_snr"],
                                "recovered": ch.recovered(base, mi, p)})  # fmt: skip
    t = Table(rows=out)
    t.meta.update(
        provenance=schema.Provenance.SIMULATED.value,
        source="w2_hst_blend_run.py psf: white noise at the control's sigma_pix, its fitted "
        "images made with a Gaia-star empirical PSF (w12_niq_epsf.build_epsf) or the fitted "
        "Moffat, optional Sersic lens convolved with the same PSF, measured by the Moffat chain",
    )
    t.write(RES / "psf_mismatch.ecsv", overwrite=True)
    fp = t[(np.asarray(t["psf"]) == "epsf") & ~np.isfinite(np.asarray(t["inject_mag"], float))]
    summ = {
        "systems_with_epsf": done,
        "wide_cutouts_tried": tries,
        "epsf_no_lens_max_snr": float(np.nanmax(fp["resid_snr"])) if len(fp) else None,
        "psf_mismatch_not_excluded": bool(len(fp) and np.nanmax(fp["resid_snr"]) >= p.detect_sigma),
    }
    for sid in done:
        for kind in ("epsf", "moffat"):
            for reff in p.inj_reff_arcsec:
                s = t[(t["system_id"] == sid) & (t["psf"] == kind) & (t["reff_arcsec"] == reff)]
                summ[f"m_lim_{sid}_{kind}_reff{reff}"] = ch.contiguous_limit(
                    [(float(x["inject_mag"]), bool(x["recovered"])) for x in s]
                )
    (RES / "psf_mismatch.json").write_text(json.dumps(summ, indent=2) + "\n")
    print(json.dumps(summ, indent=1))


def cmd_sheet(args) -> None:
    """Contact sheet: F814W stamp, N-Moffat residual, halo-corrected residual (smoothed 1 px)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image
    from scipy import ndimage

    meas = Table.read(RES / "measurements.ecsv")
    roles = args.roles.split(",")
    meas = meas[np.isin(np.asarray(meas["role"]), roles)]
    pd = Path(args.cut).parent / "panels"
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
            (sm, f"{r['status']} S/N {r['resid_snr']:.1f}",
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


def cmd_posthoc(args) -> None:
    """Post hoc diagnostics (not frozen rules; RESULT.md labels them): (1) the PSF-mismatch
    false-positive S/N of the one ePSF simulation scaled to each system by image flux and
    aperture error (ASSUMPTION: the same fractional mismatch); (2) for each "none between
    images" system, the halo-corrected residual in a 0.2" circle at the images' centroid."""
    p = ch.Params()
    meas = Table.read(RES / "measurements.ecsv")
    pm = json.loads((RES / "psf_mismatch.json").read_text())
    ref_id = pm["systems_with_epsf"][0]
    ref = meas[meas["system_id"] == ref_id][0]

    def img_flux(r):
        return sum(
            10 ** (-0.4 * (float(x) - float(r["zp_ab"]))) for x in r["image_mags"].split(";")
        )

    out = {"reference": str(ref["name"]), "fp_snr_ref": pm["epsf_no_lens_max_snr"], "systems": {}}
    for r in meas[np.isin(meas["status"], ch.DECIDED)]:
        fp = pm["epsf_no_lens_max_snr"] * img_flux(r) / img_flux(ref)
        fp *= float(ref["resid_err"]) / float(r["resid_err"])
        out["systems"][str(r["name"])] = {
            "role": str(r["role"]),
            "status": str(r["status"]),
            "resid_snr": float(r["resid_snr"]),
            "fp_snr_scaled": float(fp),
            "snr_over_fp_scaled": float(r["resid_snr"] / fp),
            "psf_alpha_at_bound": bool(float(r["psf_alpha"]) <= 1.2 + 1e-6),
        }
    sel = Table.read(RES / "selection.ecsv")
    for r in meas[meas["status"] == ch.DECIDED[1]]:
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
            centre_flux=flux,
            centre_snr=flux / err,
            centre_mag=float(cu["zp"] - 2.5 * np.log10(flux)) if flux > 0 else None,
        )
    out["provenance"] = "derived (post hoc diagnostics, not frozen rules)"
    (RES / "posthoc.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=1))


def fcol(t: Table, name: str) -> np.ndarray:
    """Float column with masked entries as NaN (a masked value is not 0)."""
    c = t[name]
    return np.asarray(c.filled(np.nan) if hasattr(c, "filled") else c, float)


def cmd_summary(args) -> None:
    """Numbers quoted in RESULT.md, from the tracked tables (``derived``)."""
    meas = Table.read(RES / "measurements.ecsv")
    inj = Table.read(RES / "injections.ecsv")
    pm = json.loads((RES / "psf_mismatch.json").read_text())
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
        for reff in sorted(set(np.asarray(inj["reff_arcsec"], float))):
            s = inj[(inj["role"] == role) & (np.asarray(inj["reff_arcsec"], float) == reff)]
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
        m = meas[meas["role"] == role]
        st = np.asarray(m["status"], str)
        out[f"{role}_status_counts"] = {k: int(np.sum(st == k)) for k in sorted(set(st))}
        lim = fcol(m, "m_lim_reff0.3")
        ok = np.isfinite(lim)
        out[f"{role}_m_lim_reff0.3"] = {
            "n": int(ok.sum()),
            "min": float(np.min(lim[ok])) if ok.any() else None,
            "median": float(np.median(lim[ok])) if ok.any() else None,
            "max": float(np.max(lim[ok])) if ok.any() else None,
        }
        dsep = fcol(m, "max_sep_arcsec")
        out[f"{role}_max_sep_range"] = [float(np.nanmin(dsep)), float(np.nanmax(dsep))]
    a = fcol(meas, "psf_alpha")
    out["psf_alpha_at_lower_bound"] = int(np.sum(np.isfinite(a) & (a <= 1.2 + 1e-6)))
    out["psf_fits"] = int(np.sum(np.isfinite(a)))
    if args.downloads_log:
        lines = [x.split("\t") for x in Path(args.downloads_log).read_text().splitlines()]
        b = [int(x[2]) for x in lines]
        out["downloads"] = {"files": len(b), "bytes": int(sum(b))}
        for kind, key in (("hapcut 250 px", "size=250"), ("hapcut 2200 px", "size=2200")):
            sel = [int(x[2]) for x in lines if key in x[1]]
            out["downloads"][kind] = {"files": len(sel), "bytes": int(sum(sel))}
    out["provenance"] = "derived"
    (RES / "summary.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=1))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("download")
    a.add_argument("--cut", required=True)
    a.add_argument("--log", required=True)
    a = sub.add_parser("measure")
    a.add_argument("--role", choices=("control", "target"), required=True)
    a.add_argument("--cut", required=True)
    a.add_argument("--only", default="")
    a.add_argument("--no-inject", action="store_true")
    a = sub.add_parser("psf")
    a.add_argument("--cut", required=True)
    a.add_argument("--log", required=True)
    a = sub.add_parser("sheet")
    a.add_argument("--cut", required=True)
    a.add_argument("--roles", default="control,target")
    a.add_argument("--out", default="contact_sheet.jpg")
    a = sub.add_parser("posthoc")
    a.add_argument("--cut", required=True)
    a = sub.add_parser("summary")
    a.add_argument("--downloads-log", default="")
    args = ap.parse_args(argv)
    cmds = {"download": cmd_download, "measure": cmd_measure, "psf": cmd_psf, "sheet": cmd_sheet}
    cmds["posthoc"] = cmd_posthoc
    cmds["summary"] = cmd_summary
    cmds[args.cmd](args)


if __name__ == "__main__":
    main()
