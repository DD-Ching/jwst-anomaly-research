"""Measure how a cross-band colour depends on source size (vetting finding, 2026-10-07).

Compares the aperture (``aper50``) and isophotal colour ``band_a - band_b`` across bins of the
``band_a`` isophotal area, for sources detected in both bands with ``band_a`` aper50 S/N above
``--min-snr``, and prints the bin medians and the Spearman correlation between log area and the
aperture colour. Re-run it after changing feature definitions to check the bias is gone.

    python scripts/feature_size_bias.py --run-dir $JWST_ANOMALY_OUTPUTS/runs/<run_id> \\
        --sample smacs0723_nircam --band-a F200W --band-b F277W
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy.stats import spearmanr

from jwst_anomaly import schema
from jwst_anomaly.features import column_as_float

MAG_ERR_TO_SNR = 2.5 / np.log(10)  # S/N ~ 1.0857 / sigma_mag


def size_bias(
    sources: Table, band_a: str, band_b: str, min_snr: float, bins: list[float]
) -> tuple[list[tuple[float, float, int, float, float]], float, float]:
    """Bin rows ``(lo, hi, n, median aperture colour, median isophotal colour)`` and Spearman."""

    def col(band: str, q: str) -> np.ndarray:
        return column_as_float(sources, schema.band_column(band, q))

    detected = (col(band_a, "detected") > 0) & (col(band_b, "detected") > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = MAG_ERR_TO_SNR / col(band_a, "aper50_abmag_err")
    area = col(band_a, "isophotal_area")
    ap = col(band_a, "aper50_abmag") - col(band_b, "aper50_abmag")
    iso = col(band_a, "isophotal_abmag") - col(band_b, "isophotal_abmag")
    ok = detected & (snr > min_snr) & np.isfinite(ap) & np.isfinite(area) & (area > 0)
    rows = []
    for lo, hi in zip(bins[:-1], bins[1:], strict=True):
        k = ok & (area >= lo) & (area < hi)
        rows.append((lo, hi, int(k.sum()), float(np.nanmedian(ap[k])), float(np.nanmedian(iso[k]))))
    rho, p = spearmanr(np.log10(area[ok]), ap[ok])
    return rows, float(rho), float(p)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--band-a", default="F200W")
    ap.add_argument("--band-b", default="F277W")
    ap.add_argument("--min-snr", type=float, default=10.0)
    ap.add_argument("--bins", nargs="+", type=float, default=[0, 50, 200, 1000, 1e9])
    args = ap.parse_args(argv)
    sources = Table.read(args.run_dir / args.sample / "sources.ecsv")
    rows, rho, p = size_bias(sources, args.band_a, args.band_b, args.min_snr, args.bins)
    a, b = args.band_a.upper(), args.band_b.upper()
    print(f"{args.run_dir.name} {args.sample}: {a}-{b}, {a} aper50 S/N > {args.min_snr:g}")
    print(f"{'isophotal area [pix]':>22} {'n':>5} {'median aper50':>14} {'median isophotal':>17}")
    for lo, hi, n, ma, mi in rows:
        print(f"{f'[{lo:g}, {hi:g})':>22} {n:5d} {ma:+14.2f} {mi:+17.2f}")
    print(f"Spearman(log area, aper50 colour) = {rho:.3f}, p = {p:.1e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
