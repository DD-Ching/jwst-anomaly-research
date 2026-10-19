"""Frozen rules of results/dimming/edge_coverage/PLAN.md and the weight-map sampler.

- ``wht_grid``: the i2d ``WHT`` extension on a ``step``-pixel grid (S3 row reads).
- ``wht_rel``: smallest WHT of the 3 x 3 grid nodes around each position / median positive grid WHT.
- ``edge_pair``: proxy-covered and ``wht_rel < 0.5`` (NaN counts as 0).
- ``gate_v`` and ``decision``: the plan's validation gate and decision rule.
"""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

WHT_REL_MIN = 0.5  # cut_low_weight threshold of dimming_screen forced (ASSUMPTION)
GATE_V_MIN = 9  # of the 11 NEXUS edge targets (ASSUMPTION)
EXPOSURE_DROP_MAX = 0.10  # decision threshold (ASSUMPTION)


def wht_grid(uri: str, step: int = 32, workers: int = 16) -> tuple[np.ndarray, WCS]:
    """Every ``step``-th row and column of ``WHT``; the WCS is the full image's (SCI header)."""
    kw = {"anon": True, "default_block_size": 2**18}
    with fits.open(uri, use_fsspec=True, fsspec_kwargs=kw) as h:
        wcs = WCS(h["SCI"].header)
        ny = h["WHT"].header["NAXIS2"]
    rows = np.arange(0, ny, step)

    def chunk(rs):
        with fits.open(uri, use_fsspec=True, fsspec_kwargs=kw) as h:
            sec = h["WHT"].section
            return [np.asarray(sec[int(r), :])[::step] for r in rs]

    parts = np.array_split(rows, workers)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        out = [row for part in pool.map(chunk, parts) for row in part]
    return np.array(out, dtype=np.float32), wcs


def wht_rel(grid: np.ndarray, wcs: WCS, ra: np.ndarray, dec: np.ndarray, step: int) -> np.ndarray:
    """Smallest WHT of the 3 x 3 grid nodes around each position over the median positive grid
    WHT; NaN WHT counts as 0 and positions off the grid get 0."""
    g = np.nan_to_num(grid, nan=0.0)
    med = float(np.median(g[g > 0])) if (g > 0).any() else np.nan
    x, y = wcs.all_world2pix(ra, dec, 0)
    gi = np.rint(np.asarray(y) / step).astype(int)
    gj = np.rint(np.asarray(x) / step).astype(int)
    out = np.full(len(gi), np.inf)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            i, j = gi + di, gj + dj
            ok = (i >= 0) & (i < g.shape[0]) & (j >= 0) & (j < g.shape[1])
            v = np.zeros(len(gi))
            v[ok] = g[i[ok], j[ok]]
            out = np.minimum(out, v)
    return out / med


def edge_pair(covered: np.ndarray, rel: np.ndarray) -> np.ndarray:
    """Proxy-covered and ``rel < WHT_REL_MIN`` (NaN counts as 0)."""
    return np.asarray(covered, bool) & (
        np.nan_to_num(np.asarray(rel, float), nan=0.0) < WHT_REL_MIN
    )


def gate_v(n_called: int) -> bool:
    """Gate V: at least GATE_V_MIN of the 11 NEXUS edge targets are called edge pairs."""
    return n_called >= GATE_V_MIN


def decision(exp_masked: np.ndarray, exp_unmasked: np.ndarray) -> str:
    """Per t_E joint veto-off exposures (rho 0.1): 'over_counts' if any drops by more than
    EXPOSURE_DROP_MAX, else 'gain_stands'."""
    m, u = np.asarray(exp_masked, float), np.asarray(exp_unmasked, float)
    drop = 1.0 - m / u
    return "over_counts" if bool(np.any(drop > EXPOSURE_DROP_MAX)) else "gain_stands"


GATE_V2_MIN_POS = 10  # reference-positive pairs needed (ASSUMPTION)
GATE_V2_SENS_MIN = 0.80  # P(edge_pair | stamp_edge) (ASSUMPTION)
GATE_V2_FALSE_MAX = 0.05  # P(edge_pair | not stamp_edge) (ASSUMPTION)


def stamp_edge(quality_flag, on_edge, frac_nan, wht_rel) -> np.ndarray:
    """Per-(source, epoch) ``forced`` stamp cuts (PLAN.md addendum): outside, or on_edge, or
    frac_nan > 0.1, or wht_rel < WHT_REL_MIN (NaN counts as low)."""
    outside = np.array([str(q) == "outside" for q in np.atleast_1d(quality_flag)])
    fn = np.nan_to_num(np.asarray(frac_nan, float), nan=1.0)
    wr = np.nan_to_num(np.asarray(wht_rel, float), nan=0.0)
    return outside | np.asarray(on_edge, bool) | (fn > 0.1) | (wr < WHT_REL_MIN)


def gate_v2(edge: np.ndarray, ref: np.ndarray) -> dict:
    """Gate V2 on paired boolean arrays: 'passed' iff >= GATE_V2_MIN_POS reference positives,
    sensitivity >= GATE_V2_SENS_MIN and false-edge rate <= GATE_V2_FALSE_MAX."""
    e, r = np.asarray(edge, bool), np.asarray(ref, bool)
    npos, nneg = int(r.sum()), int((~r).sum())
    sens = float((e & r).sum() / npos) if npos else float("nan")
    false = float((e & ~r).sum() / nneg) if nneg else float("nan")
    if npos < GATE_V2_MIN_POS:
        verdict = "undetermined"
    elif sens >= GATE_V2_SENS_MIN and false <= GATE_V2_FALSE_MAX:
        verdict = "passed"
    else:
        verdict = "failed"
    return {
        "n_ref_pos": npos,
        "n_ref_neg": nneg,
        "sensitivity": sens,
        "false_edge": false,
        "verdict": verdict,
    }
