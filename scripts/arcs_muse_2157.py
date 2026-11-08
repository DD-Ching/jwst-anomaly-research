"""Redshift of the `arcs` anti row MACS0308 2157 from the public ESO MUSE cube.

Plan: results/arcs/muse-2157/PLAN.md. The ESO data portal ignores HTTP Range requests, so the
cube is streamed once and only small boxes are kept in memory (target box, control box, a
white-light image); nothing is written to disk (Cloud disk decision).
"""

from __future__ import annotations

import argparse
import gzip
import json
from collections.abc import Iterable, Iterator
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

BLOCK = 2880
CUBE_URL = "https://dataportal.eso.org/dataportal_new/file/ADP.2020-01-04T23:56:39.902"
TARGET = {"name": "MACS0308 2157", "ra": 47.23285635900545, "dec": 26.755093856038467}
# Control: the central galaxy (NED WISEA J030855.93+264537.9);
# cluster z 0.356 (NED "MACS J0308.9+2645").
CONTROL = {"name": "MACS0308 central galaxy", "ra": 47.2330795, "dec": 26.7605368, "z_ref": 0.356}
Z_CLUSTER = 0.356
Z_PHOT = 0.756

# Rest-frame air wavelengths (Angstrom; standard values, ASSUMPTION: no line-ratio priors).
EMISSION = {
    "[OII]3727": 3727.4,
    "Hbeta": 4861.3,
    "[OIII]4959": 4958.9,
    "[OIII]5007": 5006.8,
    "Halpha": 6562.8,
    "[NII]6583": 6583.4,
    "[SII]6716": 6716.4,
}
ABSORPTION = {"CaK": 3933.7, "CaH": 3968.5, "Gband": 4304.4, "Mgb": 5175.4, "NaD": 5892.9}
SNR_FEATURE = 5.0  # ASSUMPTION: a feature counts at |S/N| >= 5 in its window


class _Reader:
    """Exact-size reads from an iterator of byte chunks."""

    def __init__(self, chunks: Iterable[bytes]):
        self._it: Iterator[bytes] = iter(chunks)
        self._buf = bytearray()

    def read(self, n: int) -> bytes:
        while len(self._buf) < n:
            try:
                self._buf += next(self._it)
            except StopIteration:
                raise EOFError(f"stream ended: wanted {n} bytes, have {len(self._buf)}") from None
        out = bytes(self._buf[:n])
        del self._buf[:n]
        return out


def _read_header(reader: _Reader) -> fits.Header:
    raw = b""
    while True:
        block = reader.read(BLOCK)
        raw += block
        cards = [block[i : i + 80] for i in range(0, BLOCK, 80)]
        if any(c.startswith(b"END     ") or c.rstrip() == b"END" for c in cards):
            return fits.Header.fromstring(raw.decode("ascii"))


def stream_cube(chunks: Iterable[bytes], boxes: dict[str, tuple[int, int, int, int]]) -> dict:
    """Stream a MUSE-like FITS (empty primary + 3-D DATA extension) plane by plane.

    ``boxes`` maps a name to 0-based (x0, x1, y0, y1) half-open pixel ranges; returns the box cubes
    (n_lambda, ny, nx), the white-light image (nanmean over planes) and the DATA header.
    """
    reader = _Reader(chunks)
    primary = _read_header(reader)
    if primary.get("NAXIS", 0) != 0:
        raise ValueError("expected an empty primary HDU")
    hdr = _read_header(reader)
    if hdr.get("NAXIS") != 3 or hdr.get("BITPIX") != -32:
        raise ValueError("expected a float32 3-D DATA extension")
    nx, ny, nl = hdr["NAXIS1"], hdr["NAXIS2"], hdr["NAXIS3"]
    out = {k: np.empty((nl, y1 - y0, x1 - x0), np.float32) for k, (x0, x1, y0, y1) in boxes.items()}
    wsum = np.zeros((ny, nx), np.float64)
    wn = np.zeros((ny, nx), np.int32)
    plane_bytes = nx * ny * 4
    for k in range(nl):
        plane = np.frombuffer(reader.read(plane_bytes), ">f4").reshape(ny, nx)
        for name, (x0, x1, y0, y1) in boxes.items():
            out[name][k] = plane[y0:y1, x0:x1]
        good = np.isfinite(plane)
        wsum[good] += plane[good]
        wn += good
    with np.errstate(invalid="ignore", divide="ignore"):
        white = (wsum / wn).astype(np.float32)
    return {"boxes": out, "white": white, "header": hdr}


def wavelengths(hdr: fits.Header) -> np.ndarray:
    k = np.arange(hdr["NAXIS3"], dtype=float)
    return hdr["CRVAL3"] + (k + 1 - hdr["CRPIX3"]) * hdr["CD3_3"]


def aperture_spectrum(cube: np.ndarray, xc: float, yc: float, r_px: float) -> np.ndarray:
    ny, nx = cube.shape[1:]
    yy, xx = np.mgrid[0:ny, 0:nx]
    m = (xx - xc) ** 2 + (yy - yc) ** 2 <= r_px**2
    return np.nansum(cube[:, m], axis=1)


def noise_spectrum(
    cube: np.ndarray,
    r_px: float,
    exclude: list[tuple[float, float, float]],
    n: int = 200,
    seed: int = 0,
) -> np.ndarray:
    """Per-wavelength robust scatter (1.4826 MAD) of aperture sums at random positions in the box.

    Positions within an ``exclude`` (x, y, radius) circle are skipped; sources elsewhere are
    suppressed by the MAD.
    """
    rng = np.random.default_rng(seed)
    ny, nx = cube.shape[1:]
    sums = []
    tries = 0
    while len(sums) < n and tries < 50 * n:
        tries += 1
        x, y = rng.uniform(r_px, nx - 1 - r_px), rng.uniform(r_px, ny - 1 - r_px)
        if any((x - ex) ** 2 + (y - ey) ** 2 < (er + r_px) ** 2 for ex, ey, er in exclude):
            continue
        sums.append(aperture_spectrum(cube, x, y, r_px))
    s = np.asarray(sums)
    med = np.nanmedian(s, axis=0)
    return 1.4826 * np.nanmedian(np.abs(s - med), axis=0)


def continuum_subtract(lam: np.ndarray, flux: np.ndarray, width_a: float = 150.0) -> np.ndarray:
    step = np.median(np.diff(lam))
    half = max(1, int(round(width_a / step / 2)))
    f = np.asarray(flux, float)
    cont = np.array([np.nanmedian(f[max(0, i - half) : i + half + 1]) for i in range(f.size)])
    return f - cont


def line_snr(
    lam: np.ndarray, resid: np.ndarray, noise: np.ndarray, center: float, half_a: float
) -> float:
    m = np.abs(lam - center) <= half_a
    if m.sum() < 2 or not np.all(np.isfinite(noise[m])):
        return float("nan")
    return float(np.nansum(resid[m]) / np.sqrt(np.nansum(noise[m] ** 2)))


def scan_redshift(
    lam: np.ndarray,
    resid: np.ndarray,
    noise: np.ndarray,
    z_grid: np.ndarray,
    sky_mask: np.ndarray | None = None,
) -> dict:
    """Score each z: summed positive emission S/N plus absorption depth S/N (lines in range)."""
    if sky_mask is not None:
        resid = np.where(sky_mask, np.nan, resid)
    lo, hi = lam[0] + 10, lam[-1] - 10
    score = np.zeros(z_grid.size)
    for i, z in enumerate(z_grid):
        s = 0.0
        for rest in EMISSION.values():
            c = rest * (1 + z)
            if lo < c < hi:
                v = line_snr(lam, resid, noise, c, 3.0 * (1 + z))
                s += max(0.0, v) if np.isfinite(v) else 0.0
        for rest in ABSORPTION.values():
            c = rest * (1 + z)
            if lo < c < hi:
                v = line_snr(lam, resid, noise, c, 6.0 * (1 + z))
                s += max(0.0, -v) if np.isfinite(v) else 0.0
        score[i] = s
    k = int(np.argmax(score))
    return {"z_best": float(z_grid[k]), "score_best": float(score[k]), "score": score}


def features_at(
    lam: np.ndarray,
    resid: np.ndarray,
    noise: np.ndarray,
    z: float,
    sky_mask: np.ndarray | None = None,
) -> dict[str, float]:
    if sky_mask is not None:
        resid = np.where(sky_mask, np.nan, resid)
    lo, hi = lam[0] + 10, lam[-1] - 10
    out = {}
    for name, rest in EMISSION.items():
        c = rest * (1 + z)
        if lo < c < hi:
            out[name] = line_snr(lam, resid, noise, c, 3.0 * (1 + z))
    for name, rest in ABSORPTION.items():
        c = rest * (1 + z)
        if lo < c < hi:
            out[name] = -line_snr(lam, resid, noise, c, 6.0 * (1 + z))  # positive = absorption
    return out


def n_features(feats: dict[str, float]) -> int:
    return sum(1 for v in feats.values() if np.isfinite(v) and v >= SNR_FEATURE)


def classify(z_best: float, n_feat: int) -> str:
    """Frozen reading of PLAN.md."""
    if n_feat >= 2 and abs(z_best - Z_CLUSTER) < 0.02:
        return "member"
    if n_feat >= 2 and z_best > 0.40:
        return "background"
    if n_feat >= 2:
        return "foreground"
    if n_feat == 1:
        return "single_line"
    return "undetermined"


def sky_line_mask(lam: np.ndarray, noise: np.ndarray, factor: float = 3.0) -> np.ndarray:
    """Pixels whose noise exceeds ``factor`` x the running median noise (sky-line residuals)."""
    base = np.array([np.nanmedian(noise[max(0, i - 40) : i + 41]) for i in range(noise.size)])
    return ~np.isfinite(noise) | (noise > factor * base)


def _http_chunks(url: str, chunk: int = 4 << 20) -> Iterator[bytes]:
    import requests

    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        yield from r.iter_content(chunk)


def _pix(w: WCS, ra: float, dec: float) -> tuple[float, float]:
    x, y = w.celestial.world_to_pixel_values(ra, dec)
    return float(x), float(y)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=Path("results/arcs/muse-2157"))
    ap.add_argument("--r-arcsec", type=float, default=0.6)
    ap.add_argument("--half-box", type=int, default=40)
    args = ap.parse_args()
    # Header first (small read of a stream we then close), to place the boxes.
    reader = _Reader(_http_chunks(CUBE_URL, 1 << 16))
    _read_header(reader)
    hdr = _read_header(reader)
    w = WCS(hdr)
    boxes, centres = {}, {}
    for key, src in (("target", TARGET), ("control", CONTROL)):
        x, y = _pix(w, src["ra"], src["dec"])
        h = args.half_box
        x0, y0 = max(0, int(round(x)) - h), max(0, int(round(y)) - h)
        x1, y1 = min(hdr["NAXIS1"], int(round(x)) + h), min(hdr["NAXIS2"], int(round(y)) + h)
        boxes[key] = (x0, x1, y0, y1)
        centres[key] = (x - x0, y - y0, x, y)
    started = datetime.now(UTC).isoformat(timespec="seconds")
    res = stream_cube(_http_chunks(CUBE_URL), boxes)
    lam = wavelengths(res["header"])
    pix_arcsec = abs(hdr["CD2_2"]) * 3600
    r_px = args.r_arcsec / pix_arcsec
    z_grid = np.arange(0.0, 1.6, 0.0002)
    report = {
        "plan": "results/arcs/muse-2157/PLAN.md",
        "cube": CUBE_URL,
        "streamed_utc": started,
        "aperture_r_arcsec": args.r_arcsec,
        "snr_feature": SNR_FEATURE,
        "provenance": {
            "observed": "MUSE cube DATA values (ESO phase-3 product)",
            "derived": "aperture spectra, noise, z scan, features, class",
            "assumption": "aperture radius, S/N 5 feature cut, line windows, sky mask factor 3",
        },
        "rows": {},
    }
    spectra = {"lam": lam.tolist()}
    for key, src in (("control", CONTROL), ("target", TARGET)):
        cube = res["boxes"][key]
        xc, yc, xf, yf = centres[key]
        spec = aperture_spectrum(cube, xc, yc, r_px)
        noise = noise_spectrum(cube, r_px, exclude=[(xc, yc, 3 * r_px)])
        mask = sky_line_mask(lam, noise)
        resid = continuum_subtract(lam, spec)
        scan = scan_redshift(lam, resid, noise, z_grid, mask)
        feats = features_at(lam, resid, noise, scan["z_best"], mask)
        nf = n_features(feats)
        row = {
            "name": src["name"],
            "ra": src["ra"],
            "dec": src["dec"],
            "pixel_xy": [xf, yf],
            "z_best": scan["z_best"],
            "score_best": scan["score_best"],
            "features_at_best": feats,
            "n_features": nf,
            "class": classify(scan["z_best"], nf),
            "features_at_z_cluster": features_at(lam, resid, noise, Z_CLUSTER, mask),
            "features_at_z_phot": features_at(lam, resid, noise, Z_PHOT, mask),
            "continuum_median_snr_per_pix": float(np.nanmedian(spec / noise)),
            "n_sky_masked": int(mask.sum()),
        }
        if key == "control":
            row["z_ref"] = CONTROL["z_ref"]
            row["gate_m_passed"] = bool(nf >= 2 and abs(scan["z_best"] - CONTROL["z_ref"]) < 0.01)
        report["rows"][key] = row
        spectra[key] = {
            "flux": spec.tolist(),
            "noise": noise.tolist(),
            "score": scan["score"].tolist(),
        }
    report["z_grid"] = [float(z_grid[0]), float(z_grid[-1]), 0.0002]
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(report, indent=1, default=float) + "\n")
    with gzip.open(args.out / "spectra.json.gz", "wt") as f:
        json.dump(spectra, f)
    _sheet(args.out, res, centres, lam, spectra, r_px)
    print(
        json.dumps(
            {
                k: {kk: v[kk] for kk in ("z_best", "n_features", "class")}
                for k, v in report["rows"].items()
            }
        )
    )


def _sheet(
    out: Path, res: dict, centres: dict, lam: np.ndarray, spectra: dict, r_px: float
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 2, figsize=(12, 10), gridspec_kw={"width_ratios": [1, 3]})
    for row, key in enumerate(("control", "target")):
        img = np.nanmean(res["boxes"][key], axis=0)
        ax = axes[row, 0]
        lo, hi = np.nanpercentile(img, [5, 99.5])
        ax.imshow(img, origin="lower", vmin=lo, vmax=hi, cmap="gray")
        xc, yc = centres[key][:2]
        ax.add_patch(plt.Circle((xc, yc), r_px, fill=False, color="r"))
        ax.set_title(f"{key}: white light")
        ax = axes[row, 1]
        f = np.asarray(spectra[key]["flux"])
        ax.plot(lam, f, lw=0.5, color="k")
        ax.plot(lam, spectra[key]["noise"], lw=0.5, color="0.6")
        for z, col in ((Z_CLUSTER, "b"), (Z_PHOT, "r")):
            for rest in list(EMISSION.values()) + list(ABSORPTION.values()):
                c = rest * (1 + z)
                if lam[0] < c < lam[-1]:
                    ax.axvline(c, color=col, lw=0.4, alpha=0.5)
        ax.set_ylim(np.nanpercentile(f, 1) - 50, np.nanpercentile(f, 99.5) + 50)
        ax.set_title(f"{key}: aperture spectrum (blue z {Z_CLUSTER}, red z {Z_PHOT})")
    ax = axes[2, 1]
    zg = np.arange(0.0, 1.6, 0.0002)
    for key, col in (("control", "b"), ("target", "r")):
        ax.plot(zg, spectra[key]["score"], color=col, lw=0.6, label=key)
    ax.legend()
    ax.set_xlabel("z")
    ax.set_title("feature score")
    axes[2, 0].imshow(
        res["white"],
        origin="lower",
        cmap="gray",
        vmin=np.nanpercentile(res["white"], 5),
        vmax=np.nanpercentile(res["white"], 99.5),
    )
    for key, col in (("control", "c"), ("target", "r")):
        axes[2, 0].plot(*centres[key][2:], "o", mfc="none", color=col)
    axes[2, 0].set_title("cube white light")
    fig.tight_layout()
    fig.savefig(out / "contact_sheet.png", dpi=80)


if __name__ == "__main__":
    main()
