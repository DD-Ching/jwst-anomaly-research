"""Residual morphology of ``arc_diff.py`` hits: edges, streaks and PSF rosettes (D-027 addendum).

Outer MACS0416 tiles (2026-10-15) showed artifact classes that the bright-core flag misses (and
residuals at the mosaic edge):

- **streak**: a diffraction spike of a bright star, rotated between the epochs, leaves a long thin
  residual of one sign; hits are strung along it;
- **rosette**: a compact source (star, compact core) fainter than the bright-core level still leaves
  alternating-sign lobes around itself after PSF matching.

A point transient changes one PSF-sized spot. For each hit and band this measures, on the
significance map S of an ``arc_diff.py --maps`` file:

- ``extent_arcsec``: the length of the connected region of sign(S_hit) × S' ≥ ``--streak-level``
  that contains the hit, where S' is S smoothed by a Gaussian of ``--smooth-px`` (it joins the
  broken segments of a spike); the length is √(12 λ_max) of the region's pixel second moments
  (the full length of a uniform line);
- ``n_lobes``, ``both_signs``: connected regions of |S| ≥ ``--level`` (at least ``--min-lobe-px``
  pixels) within ``--lobe-arcsec`` of the hit, and whether both signs occur among them;
- ``dipole_ratio``: the largest |S| within ``--lobe-arcsec`` of the hit of an opposite-sign region
  (|S| ≥ ``--level``, at least ``--min-lobe-px`` pixels) that touches the disk of
  ``--dipole-arcsec`` around the hit (its local peak, not a far peak of a long streak), over
  |S_hit| (0 if none, or if |S_hit| < ``--level``). A registration or PSF residual on a compact
  core leaves a +/− pair of similar strength; a brightening leaves one sign;
- ``edge_arcsec``: the distance to the nearest pixel without data (mosaic edge or a gap), or to the
  map border.

Only **detected** bands (|S_hit| ≥ ``--level``; the top band when none is) enter the verdict, so
noise lobes or a nearby edge in a band without the hit cannot veto it. A hit is at the **edge**
when ``edge_arcsec`` < ``--edge-arcsec`` in any detected band (a target within the search radius
of the map border is judged, not dropped); **extended** when ``extent_arcsec`` >
``--max-extent-arcsec`` in any detected band (spike streaks, also when broken up in the top band,
and the smoothed lobes of bright rosettes); a **rosette** when ``n_lobes`` ≥ ``--min-lobes`` with
both signs in at least ``--min-rosette-bands`` detected bands (with ``--rosette-lobe-ratio`` > 0,
only lobes whose peak reaches that fraction of |S_hit| count; 0 keeps the original rule); a
**dipole** when ``dipole_ratio`` ≥ ``--dipole-ratio`` in at least ``--min-dipole-bands`` bands;
otherwise it passes. ``--rejudge`` re-applies the verdict to the per-band values of an existing
output JSON (no maps needed): measurement thresholds must match the file's; the stored values keep
the definitions of the code that wrote them (``dipole_ratio`` absent in early files counts as 0;
bands dropped at a map border stay dropped). Every number is
``derived``; the thresholds are ASSUMPTIONs. Validation: every Yan+2023 arc-transient position in
the 20″ Warhol box passes in o211/o212; one (D21-S2, outside its discovery pair) is a rosette in
o211/o213.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
from pathlib import Path

import numpy as np
from arc_diff import _wcs as _wcs_header
from astropy.io import fits
from astropy.table import Table
from scipy import ndimage

_CROSS = np.ones((3, 3))


@dataclasses.dataclass(frozen=True)
class Params:
    """Morphology thresholds (every value an ASSUMPTION; CLI flags of the same name)."""

    peak_arcsec: float = 0.06
    level: float = 3.0
    smooth_px: float = 2.0
    streak_level: float = 3.0
    max_extent_arcsec: float = 3.0
    edge_arcsec: float = 0.5
    lobe_arcsec: float = 0.8
    min_lobe_px: int = 4
    min_lobes: int = 4
    min_rosette_bands: int = 2
    dipole_arcsec: float = 0.3
    dipole_ratio: float = 0.5
    min_dipole_bands: int = 2
    rosette_lobe_ratio: float = 0.0


def _wcs(text: str):
    return _wcs_header(fits.Header.fromstring(text))


def robust_sigma(v: np.ndarray) -> float:
    """1.4826 × MAD of the finite values of ``v``; 1.0 when it is not positive."""
    v = v[np.isfinite(v)]
    mad = 1.4826 * float(np.median(np.abs(v - np.median(v)))) if v.size else 0.0
    return mad if mad > 0 else 1.0


def component_extent(mask: np.ndarray, y: int, x: int, lab: np.ndarray | None = None) -> float:
    """√(12 λ_max) in pixels of the connected region of ``mask`` that contains (y, x); 0 if none.

    ``lab`` may hold a precomputed 8-connected labelling of ``mask``.
    """
    if not mask[y, x]:
        return 0.0
    if lab is None:
        lab, _ = ndimage.label(mask, structure=_CROSS)
    yy, xx = np.nonzero(lab == lab[y, x])
    if yy.size < 2:
        return 0.0
    cov = np.cov(np.vstack([xx, yy]).astype(float))
    return float(np.sqrt(12.0 * max(np.linalg.eigvalsh(cov)[-1], 0.0)))


class Regions:
    """8-connected regions of sign × S ≥ level (each sign) with their sizes and peaks."""

    def __init__(self, s: np.ndarray, level: float):
        self.level = level
        self.shape = s.shape
        self.by_sign = {}
        for sign in (1, -1):
            v = np.nan_to_num(sign * s)
            lab, n = ndimage.label(v >= level, structure=_CROSS)
            idx = np.arange(1, n + 1)
            sizes = np.bincount(lab.ravel(), minlength=n + 1)
            peaks = np.zeros(n + 1)
            if n:
                peaks[1:] = ndimage.maximum(v, lab, idx)
            self.by_sign[sign] = (lab, sizes, peaks)

    def _disk(self, y: float, x: float, r_px: float):
        y0, y1 = max(0, int(np.floor(y - r_px))), min(self.shape[0], int(np.ceil(y + r_px)) + 1)
        x0, x1 = max(0, int(np.floor(x - r_px))), min(self.shape[1], int(np.ceil(x + r_px)) + 1)
        yy, xx = np.mgrid[y0:y1, x0:x1]
        return (slice(y0, y1), slice(x0, x1)), np.hypot(yy - y, xx - x) <= r_px

    def near(self, sign: int, y: float, x: float, r_px: float, min_px: int) -> np.ndarray:
        """Labels of the ``sign`` regions with ≥ min_px pixels touching the disk r_px at (y, x)."""
        lab, sizes, _ = self.by_sign[sign]
        win, disk = self._disk(y, x, r_px)
        sub = lab[win][disk]
        k = np.unique(sub[sub > 0])
        return k[sizes[k] >= min_px]

    def local_peak(self, s: np.ndarray, sign: int, labels, y: float, x: float, r_px: float):
        """Largest sign × S of the ``labels`` regions inside the disk r_px at (y, x)."""
        lab = self.by_sign[sign][0]
        win, disk = self._disk(y, x, r_px)
        v = np.nan_to_num(sign * s[win])[disk & np.isin(lab[win], labels)]
        return float(v.max()) if v.size else 0.0

    def check(self, level: float) -> None:
        if level != self.level:
            raise ValueError(f"Regions built at level {self.level}, used at {level}")


def lobes(
    s: np.ndarray, y: float, x: float, r_px: float, level: float, min_px: int, reg=None
) -> list[int]:
    """Signs (+1/-1) of the regions of |S| ≥ level with ≥ min_px pixels touching the disk r_px."""
    return [int(np.sign(v)) for v in lobe_peaks(s, y, x, r_px, level, min_px, reg)]


def lobe_peaks(
    s: np.ndarray, y: float, x: float, r_px: float, level: float, min_px: int, reg=None
) -> list[float]:
    """Signed peak S of each region of |S| ≥ level with ≥ min_px pixels touching the disk r_px."""
    reg = reg or Regions(s, level)
    reg.check(level)
    out = []
    for sign in (1, -1):
        peaks = reg.by_sign[sign][2]
        out += [sign * float(peaks[k]) for k in reg.near(sign, y, x, r_px, min_px)]
    return out


def dipole_ratio(
    s: np.ndarray,
    y: int,
    x: int,
    r_px: float,
    level: float,
    min_px: int,
    reg=None,
    peak_r_px: float | None = None,
) -> float:
    """Largest |S| within peak_r_px (default 2 r_px) of the opposite-sign regions (≥ level,
    ≥ min_px px) that touch the disk r_px, over |S_hit|."""
    s_hit = float(s[y, x])
    # A band where the hit itself is below the level would turn any nearby noise lobe into a large
    # ratio (3 of 11 Yan+2023 transients were flagged that way), so it never counts.
    if not np.isfinite(s_hit) or abs(s_hit) < level:
        return 0.0
    reg = reg or Regions(s, level)
    reg.check(level)
    sign = -1 if s_hit >= 0 else 1
    k = reg.near(sign, y, x, r_px, min_px)
    if not k.size:
        return 0.0
    peak_r = 2 * r_px if peak_r_px is None else max(peak_r_px, r_px)
    return reg.local_peak(s, sign, k, y, x, peak_r) / abs(s_hit)


class BandMaps:
    """Per-band products shared by every target: WCS, smoothed-S labels, edge distance, regions."""

    def __init__(self, s: np.ndarray, wcs_text: str, p: Params):
        self.s = s
        self.w = _wcs(wcs_text)
        self.scale = float(np.sqrt(abs(np.linalg.det(self.w.pixel_scale_matrix)))) * 3600
        good = np.isfinite(s)
        sm = ndimage.gaussian_filter(np.where(good, s, 0.0), p.smooth_px)
        norm = ndimage.gaussian_filter(good.astype(float), p.smooth_px)
        sm = np.where(norm > 0.5, sm / np.maximum(norm, 1e-9), 0.0)
        # S is already a matched-filter output, so its noise is correlated and smoothing lowers it
        # far less than white noise (1 / (2 √π σ)). Normalise by the smoothed map's own robust
        # spread so the level is in σ of the smoothed map.
        sm /= robust_sigma(sm[good])
        self.streak = {}
        for sign in (1, -1):
            m = sign * sm >= p.streak_level
            self.streak[sign] = (m, ndimage.label(m, structure=_CROSS)[0])
        self.edt = ndimage.distance_transform_edt(np.pad(good, 1, constant_values=False))
        self.regions = Regions(s, p.level)


def morph_band(bm: BandMaps, ra: float, dec: float, p: Params) -> dict | None:
    """One band's morphology at (ra, dec); None when off the map or the S window is all NaN."""
    x0, y0 = (float(v) for v in bm.w.all_world2pix(ra, dec, 0))
    s, scale = bm.s, bm.scale
    if not (np.isfinite(x0) and np.isfinite(y0)):
        return None
    yi, xi = int(round(y0)), int(round(x0))
    if not (0 <= yi < s.shape[0] and 0 <= xi < s.shape[1]):
        return None
    r = max(1, int(round(p.peak_arcsec / scale)))
    ya, xa = max(0, yi - r), max(0, xi - r)
    win = s[ya : yi + r + 1, xa : xi + r + 1]
    if not np.isfinite(win).any():
        return None
    k = int(np.nanargmax(np.abs(win)))
    py, px = ya + k // win.shape[1], xa + k % win.shape[1]
    s_hit = float(s[py, px])
    mask, lab = bm.streak[1 if s_hit >= 0 else -1]
    pk = lobe_peaks(s, py, px, p.lobe_arcsec / scale, p.level, p.min_lobe_px, bm.regions)
    lb = [int(np.sign(v)) for v in pk]
    return {
        "s_hit": s_hit,
        "lobe_peaks": [round(v, 2) for v in pk],
        "extent_arcsec": component_extent(mask, py, px, lab) * scale,
        "n_lobes": len(lb),
        "both_signs": bool(1 in lb and -1 in lb),
        "dipole_ratio": dipole_ratio(
            s,
            py,
            px,
            p.dipole_arcsec / scale,
            p.level,
            p.min_lobe_px,
            bm.regions,
            p.lobe_arcsec / scale,
        ),
        "edge_arcsec": float(bm.edt[py + 1, px + 1]) * scale,
    }


def _rosette_band(v: dict, p: Params) -> bool:
    """≥ min_lobes lobes of both signs; with rosette_lobe_ratio > 0 only lobes whose peak reaches
    that fraction of |S_hit| count (needs ``lobe_peaks``; older outputs use n_lobes)."""
    if p.rosette_lobe_ratio <= 0 or "lobe_peaks" not in v:
        return v["n_lobes"] >= p.min_lobes and v["both_signs"]
    strong = [x for x in v["lobe_peaks"] if abs(x) >= p.rosette_lobe_ratio * abs(v["s_hit"])]
    return len(strong) >= p.min_lobes and min(strong) < 0 < max(strong)


def classify(res: dict[str, dict], p: Params) -> tuple[str, str | None]:
    """(verdict, top band); only bands with |S_hit| ≥ level (else the top band) take part."""
    if not res:
        return "not_covered", None
    top = max(res, key=lambda b: abs(res[b]["s_hit"]))
    det = [v for v in res.values() if abs(v["s_hit"]) >= p.level] or [res[top]]
    if min(v["edge_arcsec"] for v in det) < p.edge_arcsec:
        return "edge", top
    if max(v["extent_arcsec"] for v in det) > p.max_extent_arcsec:
        return "extended", top
    n_ros = sum(_rosette_band(v, p) for v in det)
    if n_ros >= p.min_rosette_bands:
        return "rosette", top
    # Outputs written before the dipole test have no dipole_ratio; they count as no dipole.
    n_dip = sum(v.get("dipole_ratio", 0.0) >= p.dipole_ratio for v in res.values())
    if n_dip >= p.min_dipole_bands:
        return "dipole", top
    return "pass", top


def _counts(rows: list[dict]) -> dict[str, int]:
    return {v: sum(r["verdict"] == v for r in rows) for v in sorted({r["verdict"] for r in rows})}


# Thresholds that shaped the stored per-band values; a rejudge may change only the others.
MEASURED = (
    "peak_arcsec",
    "level",
    "smooth_px",
    "streak_level",
    "lobe_arcsec",
    "min_lobe_px",
    "dipole_arcsec",
)


def rejudge(path: Path, p: Params) -> dict:
    """Re-apply ``classify`` to the stored per-band values of an output JSON (in place).

    Refuses when ``p`` differs from the file's thresholds in a measurement parameter."""
    d = json.loads(path.read_text())
    old = d.get("thresholds_assumption", {})
    bad = [k for k in MEASURED if k in old and old[k] != getattr(p, k)]
    if bad:
        raise ValueError(f"{path.name}: measured with other {bad}; re-run with --maps instead")
    bands = [b for r in d["targets"] for b in r["bands"].values()]
    changed = []
    for row in d["targets"]:
        verdict, top = classify(row["bands"], p)
        if verdict != row["verdict"]:
            changed.append({"id": row["id"], "old": row["verdict"], "new": verdict})
        row["verdict"], row["top_band"] = verdict, top
    d["counts"] = _counts(d["targets"])
    d["thresholds_assumption"] = dataclasses.asdict(p)
    d["rejudged"] = True
    d["dipole_measured"] = bool(bands) and all("dipole_ratio" in b for b in bands)
    path.write_text(json.dumps(d, indent=1) + "\n")
    return {"file": path.name, "counts": d["counts"], "changed": changed}


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--hits", type=Path, help="arc_diff.py hits.ecsv")
    src.add_argument("--known", type=Path, help="CSV with name, ra, dec (validation targets)")
    src.add_argument("--rejudge", type=Path, nargs="+", help="output JSONs to re-classify in place")
    ap.add_argument("--maps", type=Path, help="arc_diff.py --maps file (with --hits / --known)")
    ap.add_argument("--frame-offset", type=float, nargs=2, default=[0.0, 0.0])
    ap.add_argument("--all-hits", action="store_true", help="include near_bright_core hits")
    for f in dataclasses.fields(Params):
        kind = int if f.type == "int" else float
        ap.add_argument("--" + f.name.replace("_", "-"), type=kind, default=f.default)
    ap.add_argument("--out", type=Path)
    a = ap.parse_args(argv)
    p = Params(**{f.name: getattr(a, f.name) for f in dataclasses.fields(Params)})
    if a.rejudge:
        for path in a.rejudge:
            print(json.dumps(rejudge(path, p)))
        return
    if a.maps is None or a.out is None:
        ap.error("--maps and --out are required with --hits / --known")
    if a.hits:
        t = Table.read(a.hits)
        if not a.all_hits and "near_bright_core" in t.colnames:
            t = t[~t["near_bright_core"]]
        targets = [(str(r["hit_id"]), float(r["ra"]), float(r["dec"])) for r in t]
    else:
        from arc_diff import load_known

        targets = [
            (k["name"], k["ra"], k["dec"]) for k in load_known(a.known, tuple(a.frame_offset))
        ]
    per = [{} for _ in targets]
    with np.load(a.maps) as m:
        bands = sorted({k.rsplit("_", 1)[0] for k in m.files if k.endswith("_wcs")})
        for b in bands:  # one band's full-size products in memory at a time
            bm = BandMaps(m[f"{b}_s"], str(m[f"{b}_wcs"]), p)
            for r, (_, ra, dec) in zip(per, targets, strict=True):
                if (v := morph_band(bm, ra, dec, p)) is not None:
                    r[b] = v
            del bm
    rows = []
    for (name, ra, dec), res in zip(targets, per, strict=True):
        verdict, top = classify(res, p)
        row = {"id": name, "ra": ra, "dec": dec, "verdict": verdict, "top_band": top}
        rows.append(row | {"bands": res})
    out = {
        "provenance": "derived",
        "maps": a.maps.name,
        "thresholds_assumption": dataclasses.asdict(p),
        "counts": _counts(rows),
        "targets": rows,
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["counts"]))


if __name__ == "__main__":
    main()
