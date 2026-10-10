"""E1: GW channels with directions from the GWTC sky maps (D-074 next step).

Each GW event gets its reduced nside-32 sky map (`scripts/e1_gw_skymaps.py`). A partner X at
(ra, dec) with 1-sigma error s is classified against the GW map in the GW event's own hour-angle
frame (the map rotates in RA with a scrambled GW time, as RA does for the other catalogues under
the D-074 nulls):

- same: X lies within ``n_sigma * s`` (+ ``PIX_SLOP_DEG``) of the GW 90 % credible region;
- wide: X lies more than ``max(n_sigma * s, 0.1 deg)`` (+ slop) outside the 99 % region;
- antipodal: X's antipode lies within ``n_sigma * s`` (+ slop) of the 90 % region; X not same.

GW-GW pairs: same when the two 90 % regions overlap; wide when the 99 % regions are more than the
slop apart; antipodal when A's 90 % region overlaps the antipode of B's and the pair is not same.
Cells: 4 channels x 5 lag bins x 3 classes, `jit` null (D-074), pooled global p and analytic
Bonferroni. Sensitivity: B events moved to t_GW +- lag and to a position drawn from the GW map
(same), from outside its 99 % region (wide) or from its antipodal map (antipodal), scattered by
B's error; GW-GW is not injected.

Writes ``results/e1_events/gw_directional.json``.

  python scripts/e1_gw_directional.py [--n 1000] [--cpu 4]
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import e1_events as E  # noqa: E402

from jwst_anomaly import event_network as en  # noqa: E402

OUT = ROOT / "results" / "e1_events" / "gw_directional.json"
EVENTS_MANIFEST = ROOT / "data" / "manifests" / "e1_gw_events.ecsv"
SKYMAP_MANIFEST = ROOT / "data" / "manifests" / "e1_gw_skymaps.ecsv"
MAPS = (
    Path(os.environ.get("JWST_ANOMALY_DATA", ROOT / "data"))
    / "e1_events"
    / "gw_skymaps_nside32.npz"
)
NSIDE = 32
CHANNELS = (("GW", "GBM"), ("GW", "ICECAT"), ("GW", "CHIME"), ("GW", "GW"))
CLASSES = ("same", "wide", "antipodal")
#: ASSUMPTION: half the diagonal of an nside-32 pixel (1.83 deg on a side), added to region
#: distances.
PIX_SLOP_DEG = 1.3
#: ASSUMPTION: credible levels that define the GW region (same) and its outside (wide).
CL_SAME, CL_WIDE = 0.90, 0.99
#: ASSUMPTION: GW170817 is not in the GWTC-2.1 / -3 / -4.1 / -5.0 PE sky-map tarballs; its map is a
#: Gaussian at SSS17a (E.NGC4993) with sigma 1.05 deg (90 % area about 16 deg^2, GWTC-1).
#: Positive control only.
GW170817_SIGMA_DEG = 1.05
INJ_N = (3, 10, 30)
INJ_TRIALS = 10
P = en.Params()
NLAG = len(P.lag_edges) - 1


def _hp():
    from astropy_healpix import HEALPix

    return HEALPix(nside=NSIDE, order="nested")


def pix_vectors() -> np.ndarray:
    import astropy.units as u

    lon, lat = _hp().healpix_to_lonlat(np.arange(12 * NSIDE**2))
    lo, la = lon.to_value(u.rad), lat.to_value(u.rad)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=1)


def credible_level(prob: np.ndarray) -> np.ndarray:
    """Credible level per pixel (0 at the peak, 1 at the least probable pixel)."""
    o = np.argsort(prob)[::-1]
    cl = np.empty_like(prob, dtype=float)
    cl[o] = np.cumsum(prob[o])
    return cl


def region_distance(region: np.ndarray, vec: np.ndarray) -> np.ndarray:
    """Great-circle distance (deg) from every pixel centre to the nearest pixel of ``region``."""
    if not region.any():
        return np.full(len(vec), 180.0)
    rv = vec[region].astype(np.float32)
    best = np.full(len(vec), -1.0, dtype=np.float32)
    for k in range(0, len(rv), 512):
        best = np.maximum(best, (vec.astype(np.float32) @ rv[k : k + 512].T).max(axis=1))
    return np.degrees(np.arccos(np.clip(best, -1, 1)))


def gaussian_map(ra0, dec0, sigma_deg, vec) -> np.ndarray:
    r, d = np.radians(ra0), np.radians(dec0)
    v0 = np.array([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)])
    sep = np.degrees(np.arccos(np.clip(vec @ v0, -1, 1)))
    p = np.exp(-0.5 * (sep / sigma_deg) ** 2)
    return p / p.sum()


def long_name(name: str, mjd: float, idx) -> str:
    """The tarball event name (GWyymmdd_hhmmss, UTC) within +-1 s of ``mjd``, else ``name``."""
    from astropy.time import Time

    for ds in (0.0, -1.0, 1.0):
        dt = Time(mjd + ds / en.DAY, format="mjd").to_datetime()
        cand = dt.strftime("GW%y%m%d_%H%M%S")
        if cand in idx:
            return cand
    return name


class GWMaps:
    """Per-GW-event distance maps (deg) to the 90 % and 99 % regions, aligned with a GW Sample."""

    def __init__(self, names, npz_path=MAPS, mjd=None, manifest=None):
        z = np.load(npz_path)
        if manifest is not None:  # pinned reduced maps (data/manifests/e1_gw_skymaps.ecsv)
            from astropy.table import Table
            from e1_gw_skymaps import maps_digest

            want = Table.read(manifest, format="ascii.ecsv").meta["maps_sha256"]
            got = maps_digest(z["name"], z["prob"])
            if got != want:
                raise SystemExit(f"{npz_path}: maps sha256 {got} != pinned {want}")
        idx = {str(n): k for k, n in enumerate(z["name"])}
        if mjd is not None:  # CSV short names (GW150914) -> tarball names (GW150914_095045)
            names = [
                n if n in idx else long_name(n, t, idx) for n, t in zip(names, mjd, strict=True)
            ]
        vec = pix_vectors()
        self.d90, self.d99, self.in90, self.in99, self.prob = [], [], [], [], []
        self.has = np.zeros(len(names), bool)
        for i, n in enumerate(names):
            if n in idx:
                p = z["prob"][idx[n]].astype(float)
            elif n == E.PC_GW:
                p = gaussian_map(*E.NGC4993, GW170817_SIGMA_DEG, vec)
            else:
                p = None
            if p is None:
                p = np.full(12 * NSIDE**2, 1.0 / (12 * NSIDE**2))
            else:
                self.has[i] = True
            cl = credible_level(p)
            r90, r99 = cl <= CL_SAME, cl <= CL_WIDE
            self.in90.append(r90)
            self.in99.append(r99)
            self.prob.append(p / p.sum())
            self.d90.append(region_distance(r90, vec).astype(np.float32))
            self.d99.append(region_distance(r99, vec).astype(np.float32))
        self.d90, self.d99 = np.array(self.d90), np.array(self.d99)
        self.in90, self.in99 = np.array(self.in90), np.array(self.in99)
        self.prob = np.array(self.prob)
        self.hp = _hp()
        lon, lat = self.hp.healpix_to_lonlat(np.arange(12 * NSIDE**2))
        self.pix_ra, self.pix_dec = lon.deg, lat.deg

    def pix(self, ra, dec) -> np.ndarray:
        import astropy.units as u

        return self.hp.lonlat_to_healpix(
            np.mod(np.asarray(ra, float), 360) * u.deg,
            np.clip(np.asarray(dec, float), -90, 90) * u.deg,
        )


def rotation(gw: en.Sample, orig_mjd: np.ndarray) -> np.ndarray:
    """RA rotation (deg) of each GW map: hour angle kept under a time shift."""
    return en.gmst_deg(gw.mjd) - en.gmst_deg(orig_mjd)


def classify_x(m: GWMaps, gi, rot, ra, dec, sig, p: en.Params = P) -> np.ndarray:
    """Classes for GW(gi) x partner pairs: bit 0 same, bit 1 wide, bit 2 antipodal."""
    s = np.where(np.isfinite(sig), sig, 0.0)
    ok = np.isfinite(ra) & np.isfinite(dec)
    ra, dec = np.where(ok, ra, 0.0), np.where(ok, dec, 0.0)
    ra_g = ra - rot[gi]
    pix = m.pix(ra_g, dec)
    apix = m.pix(ra_g + 180.0, -dec)
    tol = p.n_sigma * s + PIX_SLOP_DEG
    same = m.d90[gi, pix] <= tol
    wide = m.d99[gi, pix] > np.maximum(p.n_sigma * s, p.min_wide_deg) + PIX_SLOP_DEG
    anti = (m.d90[gi, apix] <= tol) & ~same
    return np.where(ok, same * 1 + wide * 2 + anti * 4, 0)


def classify_gw_gw(m: GWMaps, i, j, rot) -> np.ndarray:
    out = np.zeros(len(i), np.int64)
    for k, (a, b) in enumerate(zip(i, j, strict=True)):
        # B's regions in A's frame: rotate B's pixel centres by rot_b - rot_a.
        dr = rot[b] - rot[a]
        r90b, r99b = m.in90[b], m.in99[b]
        pb90 = m.pix(m.pix_ra[r90b] + dr, m.pix_dec[r90b])
        pb99 = m.pix(m.pix_ra[r99b] + dr, m.pix_dec[r99b])
        pa90 = m.pix(m.pix_ra[r90b] + dr + 180.0, -m.pix_dec[r90b])
        same = bool(m.in90[a][pb90].any())
        wide = bool(m.d99[a][pb99].min() > PIX_SLOP_DEG) if len(pb99) else False
        anti = bool(m.in90[a][pa90].any()) and not same
        out[k] = same * 1 + wide * 2 + anti * 4
    return out


def count_channel(
    m, gw, orig, rot, b: en.Sample | None, p: en.Params = P, skip=(), signed=False
) -> np.ndarray:
    """Counts (n_lag, 3 classes) for GW x b (b None = GW x GW). Only GW events with a map count;
    ``skip`` holds (GW index, b index) pairs left out (known ordinary pairs, reported apart).
    ``signed``: D = N(t_b > t_GW) - N(t_b < t_GW) per cell instead of counts (b not None)."""
    e = np.asarray(p.lag_edges)
    same_cat = b is None
    bb = gw if same_cat else b
    i, j = en.pairs_within(gw.mjd, bb.mjd, e[-1], same_cat)
    keep = m.has[i] & (m.has[j] if same_cat else True)
    for gi, bj in skip:
        keep &= ~((i == gi) & (j == bj))
    i, j = i[keep], j[keep]
    dt = (bb.mjd[j] - gw.mjd[i]) * en.DAY
    lag = np.abs(dt)
    w = np.sign(dt) if signed else None
    k = np.clip(np.searchsorted(e, lag, side="right") - 1, 0, NLAG - 1)
    c = (
        classify_gw_gw(m, i, j, rot)
        if same_cat
        else classify_x(m, i, rot, bb.ra[j], bb.dec[j], bb.sigma[j], p)
    )
    out = np.zeros((NLAG, 3), np.int64)
    for q in range(3):
        sel = (c >> q) & 1 == 1
        out[:, q] = np.bincount(
            k[sel], weights=None if w is None else w[sel], minlength=NLAG
        ).astype(np.int64)
    return out


def count_all(m, s, orig, skip=None) -> np.ndarray:
    """``skip``: {partner catalogue: [(GW index, partner index), ...]} left out of the counts."""
    gw = s["GW"]
    rot = rotation(gw, orig)
    skip = skip or {}
    return np.stack(
        [
            count_channel(m, gw, orig, rot, None if b == "GW" else s[b], skip=skip.get(b, ()))
            for _, b in CHANNELS
        ]
    )


_W: dict = {}


def _init(m, s, orig, skip=None):
    os.environ["OMP_NUM_THREADS"] = "1"
    _W.update(m=m, s=s, orig=orig, skip=skip)


def _null_chunk(seeds):
    m, s, orig, skip = _W["m"], _W["s"], _W["orig"], _W["skip"]
    return np.stack(
        [
            count_all(
                m,
                {k: en.scramble_jit(v, np.random.default_rng(sd), P) for k, v in s.items()},
                orig,
                skip,
            )
            for sd in seeds
        ]
    )


def run_null(m, s, orig, n, cpu, base_seed=8_000_000, skip=None):
    seeds = np.arange(n) + base_seed
    chunks = [seeds[i : i + 25] for i in range(0, n, 25)]
    if cpu <= 1:
        _init(m, s, orig, skip)
        return np.concatenate([_null_chunk(c) for c in chunks])
    with cf.ProcessPoolExecutor(cpu, initializer=_init, initargs=(m, s, orig, skip)) as ex:
        return np.concatenate(list(ex.map(_null_chunk, chunks)))


def inject(m, gw, b, n, lag_lo, lag_hi, cls, rng, p: en.Params = P, sign=None) -> en.Sample:
    """Copy of ``b`` with up to n events moved next to distinct mapped GW anchors (provenance:
    simulated). ``sign`` +1 puts B after the GW event, -1 before, None either at random."""
    mjd, ra, dec = b.mjd.copy(), b.ra.copy(), b.dec.copy()
    used: set[int] = set()
    moved = 0
    npix = m.prob.shape[1]
    for j in rng.permutation(len(b.mjd)):
        if moved >= n:
            break
        near = [
            i
            for i in np.flatnonzero((np.abs(gw.mjd - b.mjd[j]) <= p.inject_local_days) & m.has)
            if i not in used
        ]
        if not near:
            continue
        i = int(rng.choice(near))
        lag = np.exp(rng.uniform(np.log(max(lag_lo, 1e-3)), np.log(lag_hi)))
        sgn = rng.choice((-1.0, 1.0)) if sign is None else float(sign)
        mjd[j] = gw.mjd[i] + sgn * lag / en.DAY
        if cls == "wide":
            out = np.flatnonzero(
                m.d99[i]
                > PIX_SLOP_DEG + 2 * (b.sigma[j] if np.isfinite(b.sigma[j]) else 1) * p.n_sigma
            )
            if not len(out):
                continue
            q = int(rng.choice(out))
        else:
            q = int(rng.choice(npix, p=m.prob[i]))
        r0, d0 = m.pix_ra[q], m.pix_dec[q]
        if cls == "antipodal":
            r0, d0 = r0 + 180.0, -d0
        # injections are counted against the unscrambled GW maps (rotation 0): map frame = sky
        s_ = b.sigma[j] if np.isfinite(b.sigma[j]) else 1.0
        dra, ddec = rng.normal(0, s_, 2)
        dec[j] = np.clip(d0 + ddec, -90, 90)
        ra[j] = np.mod(r0 + dra / max(np.cos(np.radians(dec[j])), 0.05), 360)
        used.add(i)
        moved += 1
    return en.Sample(b.cat, mjd, ra, dec, b.sigma, b.year)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args(argv)
    t0 = time.time()
    # This run pins its own catalogue snapshot: the GBM TAP table grows daily, so the D-074 pin
    # (data/manifests/e1_events.ecsv) no longer matches a fresh fetch.
    E.MANIFEST = EVENTS_MANIFEST
    ev = E.load(E.fetch(False))
    s = {k: en.Sample.from_table(v) for k, v in ev.items()}
    names = [str(x) for x in ev["GW"]["name"]]
    m = GWMaps(names, mjd=s["GW"].mjd, manifest=SKYMAP_MANIFEST)
    orig = s["GW"].mjd.copy()
    # positive control: GW170817 x GRB 170817A, a known ordinary pair, is left out of the family
    # and reported apart (it would otherwise dominate the 0-10 s GW-GBM cell)
    gi = names.index(E.PC_GW)
    bj = [str(x) for x in ev["GBM"]["name"]].index(E.PC_GBM)
    skip = {"GBM": [(gi, bj)]}
    obs = count_all(m, s, orig, skip)
    obs_with_pc = count_all(m, s, orig)
    gb = s["GBM"]
    pc = int(
        classify_x(
            m, np.array([gi]), np.zeros(len(names)), gb.ra[[bj]], gb.dec[[bj]], gb.sigma[[bj]]
        )[0]
    )
    null = run_null(m, s, orig, a.n, a.cpu, skip=skip)
    mask = np.ones(obs.shape, bool)
    pmin, pglob = en.global_p(obs, null, mask)
    pa = en.analytic_p(obs, null)
    ncell = int(mask.sum())
    mu, sd = null.mean(axis=0), null.std(axis=0)
    labels = en.lag_labels(P)
    rng = np.random.default_rng(8_500_000)
    rot0 = np.zeros(len(names))
    cells = []
    for c, (_, cb) in enumerate(CHANNELS):
        for k in range(NLAG):
            for q, cls in enumerate(CLASSES):
                n50 = None
                if cb != "GW":
                    for ninj in INJ_N:
                        det = 0
                        for _ in range(INJ_TRIALS):
                            sb = inject(
                                m,
                                s["GW"],
                                s[cb],
                                ninj,
                                P.lag_edges[k],
                                P.lag_edges[k + 1],
                                cls,
                                rng,
                            )
                            o = count_channel(m, s["GW"], orig, rot0, sb, skip=skip.get(cb, ()))[
                                k, q
                            ]
                            p1 = en.analytic_p(np.array([o]), null[:, c, k, q][:, None])[0]
                            det += p1 * ncell <= E.DETECT_P
                        if det >= INJ_TRIALS / 2:
                            n50 = ninj
                            break
                cells.append(
                    {
                        "channel": f"GW-{cb}",
                        "lag": labels[k],
                        "class": cls,
                        "obs": int(obs[c, k, q]),
                        "null_mean": round(float(mu[c, k, q]), 2),
                        "null_sd": round(float(sd[c, k, q]), 2),
                        "z": round(float(en.z_score(obs[c, k, q], null[:, c, k, q])), 2),
                        "analytic_p": float(pa[c, k, q]),
                        "ul95_extra_pairs": round(
                            float(
                                en.upper_limit_95(obs[c, k, q : q + 1], null[:, c, k, q : q + 1])[0]
                            ),
                            1,
                        ),
                        "n50_injected": n50,
                    }
                )
    out = {
        "test": "E1 GW channels with sky maps (jit null)",
        "n_gw": len(names),
        "n_gw_with_map": int(m.has.sum()),
        "gw_without_map": [n for n, h in zip(names, m.has, strict=True) if not h],
        "cl_same": CL_SAME,
        "cl_wide": CL_WIDE,
        "pix_slop_deg": PIX_SLOP_DEG,
        "n_sigma": P.n_sigma,
        "positive_control_GW170817_GRB170817A": {
            "same": bool(pc & 1),
            "wide": bool(pc & 2),
            "antipodal": bool(pc & 4),
            "left_out_of_family": True,
            "gw_gbm_counts_with_pair": obs_with_pc[0].tolist(),
        },
        "n_scrambles": a.n,
        "n_cells": ncell,
        "min_cell_p_pooled": pmin,
        "global_p_pooled": pglob,
        "min_analytic_p_bonferroni": min(1.0, float(pa.min()) * ncell),
        "cells": cells,
        "runtime_s": round(time.time() - t0),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print(json.dumps({k: v for k, v in out.items() if k != "cells"}, indent=1))
    for r in cells:
        print(r)
    return 0


if __name__ == "__main__":
    sys.exit(main())
