# PEARLS 1176 registration check on Abell 2744: result (2026-10-16)

- **Plan and rules:** `PLAN.md`, committed at 04:55 UTC before any measurement. Every number below is derived.
- **Revision after review.** The first version of this record applied rule 1 more loosely than written and let the F115W clip collapse to one source. Both are fixed here; the changes are listed in §5.

## 1. Catalogue registration (`scripts/reg_offsets.py`)
1176 o221 t010 was compared with UNCOVER 2561 o001 t003 in every shared band.
- **Controls:** compact, isolated sources at AB 21–25.5, inside both footprints.
- **Split:** calibration set A and held-out set B, by crc32 parity.
- **Model choice:** rule 1 picks the global shift M0 or the affine model M1, per axis as written.

| band | controls | model (rule 1) | shift at h182, ΔRA cosδ, ΔDec (″) | B median after model (″) | B robust RMS (″) |
|---|---|---|---|---|---|
| F115W | 17 | M0 | +0.025, +0.006 | −0.000, −0.008 | 0.003, 0.007 |
| F150W | 35 | M1 (Dec RMS −40 %, slope 4.6σ) | +0.024, −0.005 | +0.004, +0.003 | 0.008, 0.005 |
| F200W | 159 | M1 (Dec RMS −23 %, 8.5σ) | +0.026, −0.014 | +0.002, +0.002 | 0.011, 0.006 |
| F277W | 70 | M1 (RMS −47 %, −45 %; 4.8σ) | +0.025, −0.006 | −0.001, +0.000 | 0.007, 0.006 |
| F356W | 228 | M0 | +0.025, −0.012 | −0.000, +0.002 | 0.008, 0.009 |
| F410M | 297 | M0 | +0.018, −0.010 | +0.003, +0.000 | 0.012, 0.010 |
| F444W | 323 | M0 | +0.027, −0.013 | −0.002, −0.001 | 0.006, 0.008 |

**Baseline pair.** ALT 3516 o002 against UNCOVER in F356W, 883 controls.
- It is offset by **+0.013″ in RA** and −0.001″ in Dec.
- The plan assumed this pair had no offset. That assumption is wrong. The pair is still a valid noise reference, because its B RMS (0.009″, 0.011″) is measured after M0.

**The rules as written:**
- **Rule 1 (M1 or M0).** M1 is chosen in F150W, F200W and F277W. Its slopes are small, at most 0.007″ per arcmin. Across the 32″ h182 control box that changes the shift by ≤ 0.004″.
- **Rule 2 (pass).** Every band passes: |B median| ≤ 0.010″, and B RMS ≤ 1.5 × baseline (0.013″, 0.016″).
  - F115W rests on 7 calibration sources (5 kept by the clip) and 10 held-out sources, so it is a weak measurement.
- **Rule 3 (one shift for all bands).** The M0 shifts spread by 0.010″ in RA and 0.020″ in Dec. The Dec spread sits at
  the 0.02″ limit, so per-band shifts are used.
  - F090W has no UNCOVER counterpart, so it takes the median of the measured bands.
  - The evaluated shifts are in `offsets_1176.json`, with their sources in `offsets_1176_source.json`.

**Conclusion.** 1176 o221 sits about **0.027″ east-south-east** of UNCOVER.
- That is ≈ 0.4 pixel at the 0.063″ pixel scale of the long-wavelength (LW) mosaics, and ≈ 0.9 pixel at 0.031″ in the short-wavelength (SW) bands.
- Residual variation with position is small (≤ 0.007″ per arcmin) but significant in three bands.
- Local medians within 30″ of h182 differ from the global shift by up to about 0.01″.
- No matched control is off by more than 0.11″. **Nothing comes near 0.3″.**

## 2. Image check of the tool's WCS (`scripts/reg_imgcheck.py`)
This check was added after §1; it is not in the plan's main route. `hit_lightcurve.py` uses the FITS header WCS
(`read_box`), not the gwcs that the catalogues use. So centroids were measured at the control positions in an 80″
box around h182, in F356W.

| mosaic | controls | image − catalogue, median (″) | max (″) |
|---|---|---|---|
| 1176 | 27 | −0.0005, −0.0024 | 0.076 |
| 2561 | 22 | −0.0013, −0.0074 | 0.059 |

The tool's WCS agrees with the catalogues to ≤ 0.008″. **The ≈ 0.3″ offset that #263 saw by eye on control c63's
cutouts is not a registration offset of 1176 o221.**
- That c63 came from #263's first random draw, which was not saved, so it could not be re-examined.
- The cause is unexplained.

## 3. Effect on photometry that uses 1176 (`hit_lightcurve.py --offsets`)
h182 and its 78 arc controls from #263 were re-measured in the 8 mosaics of 1176 o221.

**Reproduction.** The uncorrected re-run reproduces #263's stored 1176 fluxes exactly: 632 of 632 rows, maximum
difference 0.0 nJy.

**Corrected minus uncorrected**, in units of the empirical error:

| bands | median | median of the absolute value | max |
|---|---|---|---|
| F200W, F277W, F356W, F444W | +0.8 to +2.0 | 2.3–3.0 | 9–29 |
| F090W–F150W, F410M | — | 1.0–1.6 | — |

So a shift of under a pixel moves on-arc aperture fluxes by several σ, as #263 suspected.

**h182 pairs, recalibrated** (`control_z_corrected.ecsv`). Pairs without 1176 are unchanged.

| pair (UNCOVER → 1176 unless noted) | z_cal | control spread |
|---|---|---|
| F356W | +0.7 → −0.5 | 3.15 → 3.61 |
| F444W | +1.6 → +1.0 | 3.84 → 3.65 |
| F410M | +1.4 → +1.2 | 2.72 → 1.99 |
| F277W | +0.9 → +0.1 | 3.50 → 3.23 |
| F200W | +0.4 → +0.3 | 3.13 → 1.99 |
| F150W | −0.2 → −0.4 | 2.49 → 2.29 |
| 1176 → 6434, F444W (only the 1176 side corrected) | −1.0 → −0.3 | — |

**What this changes:**
- **h182** stays closed under #263's rule; no pair reaches 3. This is **exploratory**: h182 and these controls had
  already been looked at, so it is not a blind test.
- **Control spread** falls by 27–36 % in F200W and F410M, and by 5–8 % in F150W, F277W and F444W.
  - In F356W it **rises** by 15 %, so registration does not explain the F356W on-arc spread.
- **Earlier work is possibly affected** at a smaller level. #261, #264 and #265 used UNCOVER → ALT pairs, and ALT o002
  carries its own uncorrected +0.013″ shift (about half the 1176 shift).
  - That includes h182's largest z_cal, +2.3 in F356W UNCOVER → ALT o002.
  - Those pairs were not re-measured.

## 4. Limits
- **Coverage.** One 1176 visit (o221 t010) and one reference (UNCOVER o001 t003) were checked; other 1176 visits are
  unmeasured.
- **Method.** Footprints are approximated by the convex hull of each catalogue, and catalogue centroids include each
  mosaic's own detection noise.
- **Image check.** It covers F356W only, in one 80″ box.
- **M1 and F115W.** M1 is evaluated at h182, which is valid within the control box. F115W is weak.
- **Thresholds.** No threshold of #261–#265 was changed.

## 5. Changes after review (post hoc, labelled)
- **Rule 1 now applied per axis, as written.** This changes F150W, F200W and F277W from M0 to M1. Their shifts at h182
  move by ≤ 0.003″, apart from F277W RA, which moves +0.003″, from +0.022″ to +0.025″.
- **Clip floor.** The clip now keeps at least max(5, ½ of A) sources (`MIN_USED`). Before, F115W's M0 had collapsed to
  1 source.
- **Corrections.**
  - The pixel scale is now 0.063″ (LW), not 0.04″.
  - The local-shift wording is now "within about 0.01″", not "equal".
  - The "not affected" claim is replaced by "possibly affected", because of ALT's own offset.
- **Re-runs.** The photometry and z_cal were re-run with the rule-chosen offsets; the numbers above are from that run.
