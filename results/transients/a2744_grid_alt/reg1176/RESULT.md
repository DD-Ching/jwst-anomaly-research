# PEARLS 1176 registration check on Abell 2744: result (2026-10-16)

The plan and its rules are in `PLAN.md`, committed at 04:55 UTC before any measurement. Every number below is derived.

## 1. Catalogue registration (`scripts/reg_offsets.py`)
1176 o221 t010 against UNCOVER 2561 o001 t003. Controls are compact, isolated, AB 21–25.5 sources inside both
footprints, split by crc32 parity into a calibration set A and a held-out set B.

| band | controls | M0 shift ΔRA cosδ, ΔDec (″) | B median after M0 (″) | B robust RMS (″) |
|---|---|---|---|---|
| F115W | 17 | +0.025, +0.004 | −0.000, −0.006 | 0.003, 0.007 |
| F150W | 35 | +0.022, −0.003 | +0.007, −0.001 | 0.009, 0.008 |
| F200W | 159 | +0.026, −0.016 | −0.000, +0.001 | 0.012, 0.007 |
| F277W | 70 | +0.022, −0.007 | −0.002, −0.002 | 0.012, 0.012 |
| F356W | 228 | +0.025, −0.012 | −0.000, +0.002 | 0.008, 0.009 |
| F410M | 297 | +0.018, −0.010 | +0.003, +0.000 | 0.012, 0.010 |
| F444W | 323 | +0.027, −0.013 | −0.002, −0.001 | 0.006, 0.008 |

**Baseline pair.** ALT 3516 o002 against UNCOVER in F356W (883 controls): shift +0.013, −0.001; B RMS 0.009, 0.011.

**Pre-written rules:**
- **Rule 1 (position dependence).** M1 does not lower B's robust RMS by ≥ 20 % on both axes in any band where it matters: F356W x 17 %, y worse. **M0, a global shift, is used.** The M1 slopes reach 3–5σ in F356W Dec, but their size is ≤ 0.005″ per arcmin.
- **Rule 2 (pass).** Every band passes: |B median| ≤ 0.010″, and B RMS ≤ 1.5 × baseline (0.013″ and 0.016″).
  - F115W and F150W rest on only 10 and 17 held-out controls.
- **Rule 3 (one shift for all bands).** The band shifts agree within 0.02″ (spread 0.009″ in RA, 0.020″ in Dec).
  - Per-band shifts are used in `offsets_1176.json`.
  - F090W has no UNCOVER counterpart, so it takes the median of the other bands.

**Conclusion.** 1176 o221 sits about **0.027″ (≈ 0.7 pixel at 0.04″) east-south-east** of UNCOVER. It is a global shift, with no
sign of a ≈ 0.3″ offset anywhere.
- No matched control is off by more than 0.11″.
- Within 30″ and 60″ of h182, the local medians equal the global shift.

## 2. Image check of the tool's WCS (`scripts/reg_imgcheck.py`; added after §1, not in the plan's main route)
`hit_lightcurve.py` uses the FITS header WCS (`read_box`), not the gwcs behind the catalogues. So centroids were
measured at the control positions in an 80″ box around h182, in F356W.

| mosaic | controls | image − catalogue, median (″) | max (″) |
|---|---|---|---|
| 1176 | 27 | −0.0005, −0.0024 | 0.076 |
| 2561 | 22 | −0.0013, −0.0074 | 0.059 |

The tool's WCS therefore agrees with the catalogues. **The ≈ 0.3″ offset that #263 saw by eye on the cutouts of
control c63 is not a registration offset of 1176 o221.**
- That c63 came from #263's first, unsaved random draw, so it could not be re-examined.
- Its cause remains unexplained. Possibilities are PSF and morphology differences on galaxy outskirts, or the stamp
  cross, which sits at the rounded pixel; that rounding is ≤ 0.02″.

## 3. Effect on the photometry that uses 1176 (`--offsets` in `hit_lightcurve.py`)
h182 and its 78 arc controls of #263 were measured in the 8 mosaics of 1176 o221.

**Reproduction.** The uncorrected re-run reproduces #263's stored 1176 fluxes exactly: 632 / 632 rows, maximum
difference 0.0 nJy.

**Corrected minus uncorrected,** in units of the empirical error:

| band | median | median of the absolute value | max |
|---|---|---|---|
| F200W | +0.9 | 2.5 | 8.6 |
| F277W | +0.7 | 2.4 | 23 |
| F356W | +1.0 | 3.0 | 29 |
| F444W | +2.0 | 2.7 | 12 |
| F090W–F150W and F410M | — | 1.1–1.6 | — |

A 0.7-pixel shift therefore moves on-arc aperture fluxes by several σ. Steep extended light is sensitive to
registration, as #263 suspected.

**Recalibrated h182 pairs (`control_z_corrected.ecsv`).** Pairs not involving 1176 are unchanged.

| pair (UNCOVER → 1176 unless noted) | z_cal | control robust std |
|---|---|---|
| F356W | +0.7 → −0.5 | 3.15 → 3.61 |
| F444W | +1.6 → +1.0 | 3.84 → 3.65 |
| F410M | +1.4 → +1.2 | 2.72 → 1.99 |
| F277W | +0.9 → +0.2 | 3.50 → 3.34 |
| F200W | +0.4 → +0.1 | 3.13 → 2.06 |
| F150W | −0.2 → −0.3 | 2.49 → 2.31 |
| 1176 → 6434, F444W | −1.0 → −0.3 | — |

**What this changes:**
- **h182.** It stays closed under #263's rule; no pair reaches 3. h182 and its controls had been looked at before
  this check, so this result is **exploratory** and not a blind test.
- **Control spread.** The correction narrows the spread in F200W and F410M by about 30 %, and in F150W, F277W and
  F444W by less than 10 %. It does not narrow it in F356W. Registration is therefore **not** the main remaining
  on-arc systematic in F356W.
- **Other work.** #261, #264 and #265 used UNCOVER and ALT only, so they are not affected.

## 4. Limits
- One 1176 visit (o221 t010) and one reference (UNCOVER o001 t003) were checked; other 1176 visits are unmeasured.
- The catalogue footprint is approximated by the convex hull of each catalogue.
- Catalogue centroids carry each mosaic's own source-detection noise.
- The image check covers F356W in one 80″ box only.
- No threshold of #261–#265 was changed.
