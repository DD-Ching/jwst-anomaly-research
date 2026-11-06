# Result: DJA v7.4 0.5″ colours vs CEERS DR1.0 (2026-11-06)

Plan: `PLAN.md` (written first). Script: `scripts/ceers_dr1_colours.py`. Numbers: `summary.json` (derived);
figure: `colour_offsets.png` (colour offset vs DJA `flux_radius`, viewed before writing this).

**Inputs** (observed catalogs, streamed to scratch and deleted after the run; cloud disk rule):
- DJA `ceers-full-grizli-v7.4-fix_phot.fits`, sha256 `741ec72f…1c28` (the D-022 pin), aperture 1 (0.5″ diameter),
  SEP bad flags removed.
- CEERS DR1.0 `ceers_cat_v1.0.fits.gz` (Cox et al. arXiv:2510.08743), sha256
  `ce321809bec14afbdc48ce1aa97bda651f25ef36d26119a8fa8cf7f1904d3504`, 203,963,363 bytes, 87,345 rows;
  `FLUX_APER` index 5 (0.5″ in the readme's aperture list), empirical errors, STAR_FLAG = 0 and BAD_REGION_FLAG = 0.
- One-to-one match within 0.2″; S/N ≥ 10 in both bands in both catalogs (ASSUMPTION thresholds, plan).

**Result** (Δ = colour_DJA − colour_CEERS, median and NMAD; σ = median combined colour error):

| colour | n | median Δ | NMAD Δ | NMAD(Δ/σ) | bright (F200W < 25): median Δ, NMAD(Δ/σ) | plan gate |
|---|---|---|---|---|---|---|
| F115W−F150W | 18,936 | +0.018 | 0.047 | 0.92 | +0.016, 1.34 | pass |
| F150W−F200W | 20,813 | −0.006 | 0.043 | 0.87 | −0.009, 1.37 | pass |
| F200W−F277W | 25,418 | −0.069 | 0.063 | 1.54 | −0.074, 5.50 | **fail** |
| F277W−F356W | 28,315 | −0.003 | 0.049 | 1.00 | −0.003, 4.27 | pass (all); bright scatter fails |
| F356W−F444W | 22,084 | −0.054 | 0.039 | 1.40 | −0.055, 5.11 | **fail** |

The t021 footprint alone (the control's pointing, approximate box) gives the same medians within 0.004 mag.
CEERS 0.5″ and Kron colours agree to ≤ 0.02 mag, so the CEERS side is aperture-stable, as PSF matching predicts.

**Per band** (`summary.json` → `per_band_dja_minus_ceers_mag`; DJA − CEERS 0.5″ magnitude, F200W 23–25,
flux_radius 2.5–4 pix, 592 sources; DJA is brighter because CEERS is convolved to the F444W PSF before
measuring): F115W −0.137,
F150W −0.159, F200W −0.154, F277W −0.083, F356W −0.074, F444W −0.014. The SW bands agree with each other; the LW
bands are 0.07 mag fainter relative to them and F444W another 0.06 mag. The figure shows the offsets shrinking
towards zero for sources with flux_radius ≳ 10 pix.

**Interpretation** (derived; the cause is the plan's first ordinary explanation, supported but not proven):
the offsets have the sign and band pattern of DJA's missing PSF homogenisation (D-010): a fixed 0.5″ aperture
loses more light in the broader LW PSFs, most in F444W, and the effect fades for large sources. It is not a
zero-point problem confined to one band, and it is not magnitude dependent (faint and bright medians agree within
0.006 mag). DJA `{band}_tot_corr` is one factor per source, identical in every band, so it cannot correct a colour.

**Consequence for the search.** H0 fails as written for two of five colours, by 0.004–0.024 mag beyond the
0.05 mag gate, and at high S/N the catalogue errors under-state the scatter by a factor ~4–5. For ranking, a
common, size-dependent offset is shared by SMACS and CEERS because both use the same DJA v7.4 recipe and similar
NIRCam PSFs; it shifts colours, not ranks, as long as the size distributions of the compared samples match.
Absolute colours (SED or photo-z interpretation, comparison with PSF-matched literature catalogs) need either
these offsets or PSF-matched photometry (CEERS DR1.0 in CEERS). No candidate or verdict depends on a 0.07 mag
colour shift (none was checked here; post hoc statement, not tested).

**Deviations from the plan (post hoc):** the bright/faint split at F200W 25, the per-band table and the
`tot_corr` check were added after the first run; the gate verdicts use only the planned quantities.

**Limitations.** Only CEERS; SMACS has no PSF-matched reference in this check. The compact/extended split at
2.5 pix (ASSUMPTION) leaves few compact sources at S/N ≥ 10 (576–1,225). The aperture index relies on the
readme's list (11 values for a 12-element column; the CEERS 0.5″ vs Kron agreement supports index 5).
