"""Spike checks for El Gordo ``arc_diff`` review entries (vetting; one-box pilot, 2026-10-16).

Subcommands (all outputs ``derived``; thresholds are ASSUMPTIONs fixed in
``results/transients/elgordo_arcdiff/PLAN.md`` and its addenda):

``headers``: PA_V3, CAL_VER and CRDS context of the four i2d mosaics (S3 headers) and the finite
fractions of the ``arc_diff.py --maps`` arrays → ``headers.json``.

``spike``: is each review entry on a diffraction spike of the box's brightest point source?
- The source is the brightest F444W pixel, located in both bands.
- Ray angles are the strongest azimuthal medians 1–3″ from it (``ray_angles``). A ray counts as
  a spike only if a ray within 3° of the opposite direction exists in either band
  (``symmetric_spikes``): spikes come in opposite pairs, an arc or a galaxy does not.
- ``on_spike``: the angle to the entry is within 2° of a spike and the distance is ≤ 6″
  (``classify``).
- ``rotation_dipole``: S sampled ±1.5 px across the matched spike line, at the line's point
  nearest the entry (``across_line``), has opposite signs with each |S| ≥ 2.
- The PA_V3 difference ``dpa`` moves a spike sideways by d·tan(dpa) at distance d.

``onspike-grid`` / ``onspike-classify``: on-spike controls for an entry (PLAN addendum 1).
- Grid: 10 equal steps over 2.0–5.0″ from the source along the entry's direction; the first 8
  points not within 0.3″ of the entry are used.
- The entry is ``explained by the spike shift`` if ≥ 4 of the 8 controls change with the entry's
  sign at ≥ 2σ in the chosen band (``onspike_class``), else ``undetermined``.

    python scripts/elgordo_spike_check.py spike --maps M.npz --review review.json \\
        --dpa 0.842 --out spike.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

PEAK_HALF_PX = 3  # local-peak search half-width for the spike source in each band
SPIKE_TOL_DEG = 2.0  # ASSUMPTION: entry on a spike
PAIR_TOL_DEG = 3.0  # ASSUMPTION: opposite ray of a spike
MAX_DIST_ARCSEC = 6.0  # ASSUMPTION
ACROSS_PX = 1.5  # ASSUMPTION
DIPOLE_MIN_S = 2.0  # ASSUMPTION


def _angdiff(a: float, b: float) -> float:
    d = abs(a - b) % 360
    return min(d, 360 - d)


def ray_angles(img: np.ndarray, x0: float, y0: float, r_px: tuple[float, float], n: int = 8):
    """Angles (deg, pixel frame) of the ``n`` strongest 1°-bin azimuthal medians, ≥ 10° apart."""
    yy, xx = np.indices(img.shape, dtype=float)
    r = np.hypot(xx - x0, yy - y0)
    m = (r >= r_px[0]) & (r <= r_px[1]) & np.isfinite(img)
    ang = (np.degrees(np.arctan2(yy[m] - y0, xx[m] - x0)) % 360).astype(int)
    vals = img[m]
    prof = np.full(360, -np.inf)
    for k in np.unique(ang):
        prof[k] = np.median(vals[ang == k])
    out: list[float] = []
    for i in np.argsort(prof)[::-1]:
        if not np.isfinite(prof[i]):
            break
        a = float(i) + 0.5
        if all(_angdiff(a, b) > 10 for b in out):
            out.append(a)
        if len(out) == n:
            break
    return sorted(out)


def symmetric_spikes(
    rays_by_band: dict[str, list[float]], tol: float = PAIR_TOL_DEG
) -> list[float]:
    """Rays (any band) that have a ray within ``tol`` of the opposite direction in any band."""
    allr = [a for v in rays_by_band.values() for a in v]
    keep = []
    for a in sorted(set(allr)):
        if any(_angdiff(a + 180.0, b) <= tol for b in allr):
            keep.append(a)
    return keep


def classify(angle_to_entry: float, dist_arcsec: float, spikes: list[float]) -> tuple[bool, float]:
    if not spikes:
        return False, float("nan")
    d = min(_angdiff(angle_to_entry, s) for s in spikes)
    return bool(d <= SPIKE_TOL_DEG and dist_arcsec <= MAX_DIST_ARCSEC), float(d)


def nearest_spike(angle: float, spikes: list[float]) -> float:
    return min(spikes, key=lambda s: _angdiff(angle, s))


def across_line(S, x0, y0, spike_deg, x, y, off_px: float = ACROSS_PX):
    """S at ±``off_px`` across the spike line (through the source at ``spike_deg``), at the line's
    point nearest (x, y). Returns None when a sample falls outside the map."""
    ux, uy = np.cos(np.radians(spike_deg)), np.sin(np.radians(spike_deg))
    t = (x - x0) * ux + (y - y0) * uy
    lx, ly = x0 + t * ux, y0 + t * uy
    nx, ny = -uy, ux
    vals = []
    for s in (-off_px, off_px):
        xi, yi = int(round(lx + s * nx)), int(round(ly + s * ny))
        if not (0 <= yi < S.shape[0] and 0 <= xi < S.shape[1]):
            return None
        vals.append(float(S[yi, xi]))
    return vals


def onspike_grid(src_ra, src_dec, ent_ra, ent_dec) -> list[tuple[float, float]]:
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    s0, e = SkyCoord(src_ra, src_dec, unit="deg"), SkyCoord(ent_ra, ent_dec, unit="deg")
    pa, dist = s0.position_angle(e), s0.separation(e).arcsec
    pts = []
    for d in np.linspace(2.0, 5.0, 10):
        if abs(d - dist) < 0.3:
            continue
        c = s0.directional_offset_by(pa, d * u.arcsec)
        pts.append((float(c.ra.deg), float(c.dec.deg)))
        if len(pts) == 8:
            break
    return pts


def onspike_class(entry_z: float, control_z: list[float], k_min: int = 4, z_min: float = 2.0):
    sign = np.sign(entry_z)
    n_same = int(sum(1 for z in control_z if np.sign(z) == sign and abs(z) >= z_min))
    cls = "explained by the spike shift" if n_same >= k_min else "undetermined"
    return n_same, cls


def _source(m):
    w4 = WCS(fits.Header.fromstring(str(m["F444W_wcs"])))
    j4, i4 = np.unravel_index(np.nanargmax(m["F444W_e1img"]), m["F444W_e1img"].shape)
    return w4.pixel_to_world(i4, j4)


def cmd_spike(a) -> dict:
    m = np.load(a.maps, allow_pickle=True)
    entries = json.loads(a.review.read_text())["review"]
    src_sky = _source(m)
    geo, rays = {}, {}
    for b in ("F150W", "F444W"):
        w = WCS(fits.Header.fromstring(str(m[b + "_wcs"])))
        img = m[b + "_e1img"]
        xs, ys = w.world_to_pixel(src_sky)
        xi, yi = int(round(float(xs))), int(round(float(ys)))
        h = PEAK_HALF_PX
        y_lo, x_lo = max(yi - h, 0), max(xi - h, 0)
        win = img[y_lo : yi + h + 1, x_lo : xi + h + 1]
        dj, di = np.unravel_index(np.nanargmax(win), win.shape)
        y0, x0 = y_lo + dj, x_lo + di
        px = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))) * 3600
        geo[b] = (w, float(x0), float(y0), px)
        rays[b] = ray_angles(img, float(x0), float(y0), (1.0 / px, 3.0 / px))
    spikes = symmetric_spikes(rays)
    out = {"dpa_deg": a.dpa, "spikes_pixel_deg": spikes, "bands": {}, "entries": {}}
    for b, (w, x0, y0, _) in geo.items():
        src = w.pixel_to_world(x0, y0)
        out["bands"][b] = {
            "source_radec": [round(src.ra.deg, 6), round(src.dec.deg, 6)],
            "ray_angles_pixel_deg": rays[b],
        }
    for e in entries:
        rec = out["entries"].setdefault(e["id"], {})
        for b, (w, x0, y0, px) in geo.items():
            S = m[b + "_s"]
            x, y = w.world_to_pixel_values(e["ra"], e["dec"])
            ang = float(np.degrees(np.arctan2(y - y0, x - x0)) % 360)
            dist = float(np.hypot(x - x0, y - y0) * px)
            on, dang = classify(ang, dist, spikes)
            rb = {
                "angle_from_source_deg": round(ang, 1),
                "dist_arcsec": round(dist, 2),
                "nearest_spike_offset_deg": round(dang, 1),
                "on_spike": on,
                "expected_shift_arcsec": round(dist * np.tan(np.radians(abs(a.dpa))), 3),
            }
            if on:
                v = across_line(S, x0, y0, nearest_spike(ang, spikes), float(x), float(y))
                rb["S_across_spike"] = None if v is None else [round(t, 2) for t in v]
                rb["rotation_dipole"] = bool(
                    v is not None and v[0] * v[1] < 0 and min(map(abs, v)) >= DIPOLE_MIN_S
                )
            rec[b] = rb
        rec["on_spike_any_band"] = any(rec[b]["on_spike"] for b in geo)
        rec["rotation_dipole_any_band"] = any(rec[b].get("rotation_dipole", False) for b in geo)
    return out


def cmd_headers(a) -> dict:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
    from jwst_anomaly import cutouts

    out = {"mosaics": {}, "maps_finite_fraction": {}}
    for uri in a.uri:
        with cutouts._open_image(cutouts._resolve_uri(uri), None, 2**16) as (h, _):
            p, s = h[0].header, h["SCI"].header
            out["mosaics"][uri.rsplit("/", 1)[1]] = {
                "PA_V3": p.get("PA_V3", s.get("PA_V3")),
                "CAL_VER": p.get("CAL_VER"),
                "CRDS_CTX": p.get("CRDS_CTX"),
                "DATE-BEG": p.get("DATE-BEG"),
            }
    if a.maps:
        m = np.load(a.maps, allow_pickle=True)
        for k in m.keys():
            if m[k].ndim == 2:
                out["maps_finite_fraction"][k] = round(float(np.isfinite(m[k]).mean()), 4)
    return out


def cmd_onspike_classify(a) -> dict:
    from astropy.table import Table

    t = Table.read(a.lc)
    out = {}
    for ent in a.entry:

        def dz(name, band):
            r1 = t[(t["name"] == name) & (t["band"] == band) & (t["program"] == a.programs[0])]
            r2 = t[(t["name"] == name) & (t["band"] == band) & (t["program"] == a.programs[1])]
            d = float(r2["flux_njy"][0] - r1["flux_njy"][0])
            e = float(np.hypot(r2["err_njy"][0], r1["err_njy"][0]))
            return {"delta_njy": round(d, 1), "err": round(e, 1), "z": round(d / e, 2)}

        names = sorted({str(n) for n in t["name"] if str(n).startswith(ent + "_c")})
        rows = {n: {b: dz(n, b) for b in ("F150W", "F444W")} for n in [ent] + names}
        n_same, cls = onspike_class(rows[ent][a.band]["z"], [rows[n][a.band]["z"] for n in names])
        out[ent] = rows | {
            "_rule": {"band": a.band, "n_controls": len(names), "n_same_sign_ge2sigma": n_same}
            | {"class": cls}
        }
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("spike")
    s.add_argument("--maps", type=Path, required=True)
    s.add_argument("--review", type=Path, required=True)
    s.add_argument("--dpa", type=float, required=True, help="PA_V3 difference (deg), headers")
    s.add_argument("--out", type=Path, required=True)
    h = sub.add_parser("headers")
    h.add_argument("--uri", nargs="+", required=True)
    h.add_argument("--maps", type=Path)
    h.add_argument("--out", type=Path, required=True)
    g = sub.add_parser("onspike-grid")
    g.add_argument("--source", type=float, nargs=2, required=True)
    g.add_argument(
        "--entry", nargs=3, action="append", required=True, metavar=("NAME", "RA", "DEC")
    )
    c = sub.add_parser("onspike-classify")
    c.add_argument("--lc", type=Path, required=True)
    c.add_argument("--entry", action="append", required=True)
    c.add_argument("--band", default="F444W")
    c.add_argument("--programs", nargs=2, default=["01176", "06882"])
    c.add_argument("--out", type=Path, required=True)
    a = p.parse_args(argv)
    if a.cmd == "onspike-grid":
        args = []
        for name, ra, dec in a.entry:
            pts = onspike_grid(a.source[0], a.source[1], float(ra), float(dec))
            args += [f"--pos {name}_c{k + 1} {r:.7f} {d:.7f}" for k, (r, d) in enumerate(pts)]
            args.append(f"--pos {name} {ra} {dec}")
        print(" ".join(args))
        return 0
    out = {"spike": cmd_spike, "headers": cmd_headers, "onspike-classify": cmd_onspike_classify}[
        a.cmd
    ](a)
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
