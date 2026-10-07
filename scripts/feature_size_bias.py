"""Measure how a cross-band colour depends on source size (vetting finding, 2026-10-07).

For sources detected in both bands with ``band_a`` pipeline aper50 S/N above ``--min-snr``, prints
the median colour ``band_a - band_b`` per bin of ``band_a`` isophotal area and the Spearman
correlation between log area and colour, for the pipeline ``aper50`` colour and the isophotal
colour and, with ``--compare <label>``, the matched-aperture colour joined from the run's
``photometry`` table (D-013). With ``--compare`` every statistic uses the same rows (finite in all
colours), so differences come from the photometry, not from the selection.

    python scripts/feature_size_bias.py --run-dir $JWST_ANOMALY_OUTPUTS/runs/<run_id> \\
        --sample smacs0723_nircam --compare dja05
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table, join
from scipy.stats import spearmanr

from jwst_anomaly import schema
from jwst_anomaly.features import column_as_float

MAG_ERR_TO_SNR = 2.5 / np.log(10)  # S/N ~ 1.0857 / sigma_mag


def colour_table(
    sources: Table, band_a: str, band_b: str, min_snr: float, labels: list[str]
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Area and colours (``aper50``, ``isophotal`` and each label) on the common selected rows."""

    def col(band: str, q: str) -> np.ndarray:
        return column_as_float(sources, schema.band_column(band, q))

    detected = (col(band_a, "detected") > 0) & (col(band_b, "detected") > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = MAG_ERR_TO_SNR / col(band_a, "aper50_abmag_err")
    area = col(band_a, "isophotal_area")
    colours = {
        name: col(band_a, f"{name}_abmag") - col(band_b, f"{name}_abmag")
        for name in ["aper50", "isophotal", *labels]
    }
    ok = detected & (snr > min_snr) & np.isfinite(area) & (area > 0)
    for values in colours.values():
        ok &= np.isfinite(values)
    return area[ok], {k: v[ok] for k, v in colours.items()}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--band-a", default="F200W")
    ap.add_argument("--band-b", default="F277W")
    ap.add_argument("--min-snr", type=float, default=10.0)
    ap.add_argument("--compare", action="append", default=[], help="matched-photometry label")
    ap.add_argument("--bins", nargs="+", type=float, default=[0, 50, 200, 1000, 1e9])
    args = ap.parse_args(argv)
    sources = Table.read(args.run_dir / args.sample / "sources.ecsv")
    if args.compare:
        phot = Table.read(args.run_dir / args.sample / "photometry.ecsv")
        sources = join(sources, phot, keys="source_uid", join_type="left")
    area, colours = colour_table(sources, args.band_a, args.band_b, args.min_snr, args.compare)
    a, b = args.band_a.upper(), args.band_b.upper()
    print(
        f"{args.run_dir.name} {args.sample}: {a}-{b}; {a} aper50 S/N > {args.min_snr:g}; "
        f"{len(area)} common rows"
    )
    names = list(colours)
    print(f"{'isophotal area [pix]':>22} {'n':>5} " + " ".join(f"{n:>10}" for n in names))
    for lo, hi in zip(args.bins[:-1], args.bins[1:], strict=True):
        k = (area >= lo) & (area < hi)
        meds = " ".join(f"{np.median(colours[n][k]) if k.any() else np.nan:+10.2f}" for n in names)
        print(f"{f'[{lo:g}, {hi:g})':>22} {int(k.sum()):5d} {meds}")
    for n in names:
        rho, p = spearmanr(np.log10(area), colours[n])
        print(f"Spearman(log area, {n} colour) = {rho:+.3f}, p = {p:.1e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
