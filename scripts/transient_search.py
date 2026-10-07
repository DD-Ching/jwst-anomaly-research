"""All-source two-epoch search for variable, appearing and disappearing sources (time domain, M3).

Compares two level-3 pipeline catalogs of the same filter taken at different epochs (e.g. program
2736 in 2022 against the VENUS program 6882 in 2026 for SMACS 0723). Every result is ``derived``.
- **Frame and zero-point tie:** the global median position offset of bright pairs is removed
  before the final match. Each matched source's magnitude is tied to the median offset of
  bright pairs within ``--tie-radius`` arcsec (local, so zero-point gradients cancel). Pairs
  without enough local references cannot be judged; they are counted in
  ``meta['n_without_local_tie']``, not silently passed.
- **variable:** a matched source whose magnitude change exceeds ``--min-dmag`` and ``--min-sigma``
  times its error (catalog errors in quadrature, plus ``--sys-floor`` mag). Candidate caustic
  crossings, AGN or supernovae.
- **appeared:** a source of the later catalog with S/N >= ``--min-snr`` and no counterpart within
  ``--match-radius`` in the earlier one, inside the earlier footprint and bright enough to have been
  detected there.
- **disappeared:** the reverse, judged against the later catalog's depth.

Footprint and depth are catalog-based proxies, so every candidate needs a look at both images
(``--cutouts``). Point-like sources and lensed regions are where caustic transients are expected.

    python scripts/transient_search.py \\
        --epoch1 $JWST_ANOMALY_DATA/cache/mast/<obs1>/<obs1>_cat.ecsv \\
        --epoch2-obs-id jw06882-o057_t057_nircam_clear-f444w --out outputs/transients/f444w
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
from epoch_compare import _fetch_cats, mutual_matches  # noqa: E402

from jwst_anomaly import catalog, paths  # noqa: E402
from jwst_anomaly.features import snr_from_mag_err  # noqa: E402


def local_offsets(
    c1: SkyCoord,
    c2: SkyCoord,
    dmag: np.ndarray,
    ref: np.ndarray,
    tie_radius_arcsec: float,
    min_refs: int = 10,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per matched pair: the local median (dRA*cosDec, dDec) offset in mas and Δmag of the
    reference pairs within ``tie_radius_arcsec`` (excluding the pair itself); NaN when fewer than
    ``min_refs`` references are near. ``c1``/``c2``/``dmag``/``ref`` are aligned per pair."""
    dra, ddec = c1.spherical_offsets_to(c2)
    off = np.column_stack([dra.to_value(u.mas), ddec.to_value(u.mas)])
    ir = np.flatnonzero(ref)
    med_off = np.full((len(c1), 2), np.nan)
    med_dm = np.full(len(c1), np.nan)
    if ir.size == 0:
        return off, med_off, med_dm
    idx_pair, idx_ref, _, _ = c1[ir].search_around_sky(c1, tie_radius_arcsec * u.arcsec)
    groups: dict[int, list[int]] = {}
    for p, r in zip(idx_pair, idx_ref, strict=True):
        if ir[r] != p:
            groups.setdefault(int(p), []).append(int(ir[r]))
    for p, refs in groups.items():
        if len(refs) >= min_refs:
            med_off[p] = np.median(off[refs], axis=0)
            med_dm[p] = np.median(dmag[refs])
    return off, med_off, med_dm


def search(
    cat1: Table,
    cat2: Table,
    *,
    match_radius_arcsec: float = 0.3,
    tie_radius_arcsec: float = 60.0,
    ref_snr: float = 30.0,
    min_snr: float = 10.0,
    min_dmag: float = 0.3,
    min_sigma: float = 5.0,
    sys_floor: float = 0.03,
    footprint_radius_arcsec: float = 5.0,
    footprint_min: int = 3,
) -> Table:
    """Candidates (``variable``, ``appeared``, ``disappeared``) from two same-filter catalogs."""
    c1 = SkyCoord(cat1["ra"], cat1["dec"], unit="deg")
    c2 = SkyCoord(cat2["ra"], cat2["dec"], unit="deg")
    m1 = np.asarray(cat1["aper50_abmag"], float)
    m2 = np.asarray(cat2["aper50_abmag"], float)
    e1 = np.asarray(cat1["aper50_abmag_err"], float)
    e2 = np.asarray(cat2["aper50_abmag_err"], float)
    s1, s2 = snr_from_mag_err(e1), snr_from_mag_err(e2)
    # Frame tie: remove the median offset of bright pairs before the final match, so a systematic
    # offset between programs cannot turn matches into appeared/disappeared pairs.
    j1, j2 = mutual_matches(c1, c2, match_radius_arcsec)
    bright = (s1[j1] >= ref_snr) & (s2[j2] >= ref_snr)
    shift = np.zeros(2)
    if bright.sum() >= 10:
        dra, ddec = c1[j1[bright]].spherical_offsets_to(c2[j2[bright]])
        shift = np.array([np.median(dra.to_value(u.arcsec)), np.median(ddec.to_value(u.arcsec))])
        c2 = c2.spherical_offsets_by(-shift[0] * u.arcsec, -shift[1] * u.arcsec)
    i1, i2 = mutual_matches(c1, c2, match_radius_arcsec)
    dm = m2[i2] - m1[i1]
    ref = (s1[i1] >= ref_snr) & (s2[i2] >= ref_snr) & np.isfinite(dm)
    _, _, med_dm = local_offsets(c1[i1], c2[i2], dm, ref, tie_radius_arcsec)
    n_no_tie = int((~np.isfinite(med_dm)).sum())  # matched pairs without a local zero point
    rows: list[dict[str, Any]] = []

    # variable: matched, significant zero-point-corrected magnitude change
    dm_corr = dm - med_dm
    sigma = np.hypot(np.hypot(e1[i1], e2[i2]), sys_floor)
    with np.errstate(invalid="ignore"):
        var = (np.abs(dm_corr) >= min_dmag) & (np.abs(dm_corr) / sigma >= min_sigma)
        var &= (s1[i1] >= min_snr) | (s2[i2] >= min_snr)
    for k in np.flatnonzero(var):
        rows.append(
            _row(
                "variable",
                c1[i1[k]],
                m1[i1[k]],
                m2[i2[k]],
                s1[i1[k]],
                s2[i2[k]],
                dm_corr[k],
                dm_corr[k] / sigma[k],
                int(cat1["label"][i1[k]]),
                int(cat2["label"][i2[k]]),
            )
        )

    # appeared / disappeared: unmatched, bright, inside the other catalog's footprint and depth
    matched1 = np.zeros(len(c1), bool)
    matched1[i1] = True
    matched2 = np.zeros(len(c2), bool)
    matched2[i2] = True
    for kind, ca, cb, ma, sa, matched_a, mb, sb in (
        ("appeared", c2, c1, m2, s2, matched2, m1, s1),
        ("disappeared", c1, c2, m1, s1, matched1, m2, s2),
    ):
        # depth of catalog b: magnitude where its detections reach S/N 10 (median of 8-12)
        near10 = (sb >= 8) & (sb <= 12) & np.isfinite(mb)
        depth_b = float(np.median(mb[near10])) if near10.any() else np.inf
        cand = np.flatnonzero(~matched_a & (sa >= min_snr) & np.isfinite(ma) & (ma < depth_b - 0.5))
        if cand.size == 0:
            continue
        # footprint proxy: catalog b has sources around the position
        _, counts = _neighbour_counts(ca[cand], cb, footprint_radius_arcsec)
        for k, n in zip(cand, counts, strict=True):
            if n >= footprint_min:
                rows.append(
                    _row(
                        kind,
                        ca[k],
                        ma[k] if kind == "disappeared" else np.nan,
                        ma[k] if kind == "appeared" else np.nan,
                        sa[k] if kind == "disappeared" else np.nan,
                        sa[k] if kind == "appeared" else np.nan,
                        np.nan,
                        np.nan,
                        int(cat1["label"][k]) if kind == "disappeared" else -1,
                        int(cat2["label"][k]) if kind == "appeared" else -1,
                        depth=depth_b,
                    )
                )
    out = Table(rows=rows or None, names=_COLUMNS, dtype=_DTYPES)
    out.meta.update(
        provenance="derived",
        n_matched=len(i1),
        n_references=int(ref.sum()),
        frame_shift_arcsec=[float(x) for x in shift],
        n_without_local_tie=n_no_tie,
        thresholds={
            "match_radius_arcsec": match_radius_arcsec,
            "tie_radius_arcsec": tie_radius_arcsec,
            "ref_snr": ref_snr,
            "min_snr": min_snr,
            "min_dmag": min_dmag,
            "min_sigma": min_sigma,
            "sys_floor": sys_floor,
            "footprint": [footprint_radius_arcsec, footprint_min],
            "provenance": "assumption",
        },
    )
    return out


_COLUMNS = (
    "kind",
    "ra",
    "dec",
    "mag1",
    "mag2",
    "snr1",
    "snr2",
    "dmag",
    "dmag_sigma",
    "label1",
    "label2",
    "depth_other",
)
_DTYPES = (str, float, float, float, float, float, float, float, float, int, int, float)


def _row(kind, c, mag1, mag2, snr1, snr2, dmag, sig, label1, label2, depth=np.nan) -> dict:
    return dict(
        zip(
            _COLUMNS,
            (
                kind,
                float(c.ra.deg),
                float(c.dec.deg),
                float(mag1),
                float(mag2),
                float(snr1),
                float(snr2),
                float(dmag),
                float(sig),
                label1,
                label2,
                float(depth),
            ),
            strict=True,
        )
    )


def _neighbour_counts(pos: SkyCoord, other: SkyCoord, radius_arcsec: float):
    """Number of ``other`` sources within ``radius_arcsec`` of each ``pos`` (footprint proxy)."""
    # SkyCoord.search_around_sky(arg) returns indices into ``arg`` first, then into ``self``.
    idx_pos, _, _, _ = other.search_around_sky(pos, radius_arcsec * u.arcsec)
    return idx_pos, np.bincount(idx_pos, minlength=len(pos))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--epoch1", type=Path, required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--epoch2", type=Path)
    g.add_argument("--epoch2-obs-id")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, default=None)
    ap.add_argument("--min-dmag", type=float, default=0.3)
    ap.add_argument("--min-sigma", type=float, default=5.0)
    ap.add_argument("--min-snr", type=float, default=10.0)
    args = ap.parse_args(argv)
    manifest = args.manifest or paths.outputs_dir() / "vet" / "epoch_manifest.ecsv"
    cat1 = catalog.load_pipeline_catalog(args.epoch1)
    path2 = args.epoch2 or _fetch_cats(args.epoch2_obs_id, True, manifest)[0]
    cat2 = catalog.load_pipeline_catalog(path2)
    if cat1.meta["band"] != cat2.meta["band"]:
        ap.error(f"filters differ: {cat1.meta['band']} vs {cat2.meta['band']}")
    res = search(cat1, cat2, min_dmag=args.min_dmag, min_sigma=args.min_sigma, min_snr=args.min_snr)
    res.meta.update(
        band=cat1.meta["band"],
        epoch1={"file": Path(args.epoch1).name, "sha256": cat1.meta.get("input_sha256")},
        epoch2={"file": Path(path2).name, "sha256": cat2.meta.get("input_sha256")},
    )
    args.out.mkdir(parents=True, exist_ok=True)
    res.write(args.out / "transients.ecsv", overwrite=True)
    summary = {k: int((res["kind"] == k).sum()) for k in ("variable", "appeared", "disappeared")}
    print(
        json.dumps(
            {
                "band": res.meta["band"],
                "n_matched": res.meta["n_matched"],
                "n_references": res.meta["n_references"],
                "n_without_local_tie": res.meta["n_without_local_tie"],
                "frame_shift_arcsec": res.meta["frame_shift_arcsec"],
                **summary,
            },
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
