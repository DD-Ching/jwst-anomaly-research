# W1 field-wide row/column-segment veto for the CANUCS flanking fields (plan, frozen 2026-10-18 before any mask)

Cloud routine, claim PR #325. Handoff of `results/w1_streak_rowdiag/RESULT.md`: the D-097 trail finder misses the
detector-wide horizontal row striping through the macs0417-ncf S_max aperture, at any threshold. The post hoc V1r
of `w1_peak_vet.py` (stripe rows of a 30″ cutout) removed that peak (S 1.35), but it was local and post hoc.
This plan makes it a field-wide veto, applied identically to the screen and the injections.

**Hypothesis (H).** As `results/w1_volume/PLAN.md`: W1 lenses give a radial aperture-mass peak.
**Ordinary explanation tested.** Low-level stripes (1/f-type row or column offsets) bias the moments of the sources
on them along the stripe, giving a coherent shear pattern. **What changes the plan:** see "Outcomes".

**Veto (frozen; every number an ASSUMPTION; function `w1_streak_mask.segment_stripes`, unit-tested before the run).**
On the block-4 median image the trail finder already reads (`stream_binned`, 0.125″ pixels):
1. Cut every row into consecutive segments of `L` pixels (the last partial segment dropped). A segment needs
   ≥ 50 % finite pixels; its value is the median of its finite pixels.
2. Detrend each segment column (the segments at the same x) along y with a running median over 41 rows (≈ 5″),
   edges by nearest-value padding; residual = value − running median.
3. σ = 1.4826 × MAD of all finite residuals in the field. A segment is a stripe if |residual| > 5σ.
4. The same on columns (the image transposed): rows and columns are handled identically.
5. The stripe pixels are OR-ed onto the trail mask's 0.25″ grid (a mask pixel is set if any of its 2 × 2 block-4
   pixels is a stripe pixel). The union of trail mask and stripe mask is the field's mask; `inject_shear.py
   --streak-masks` uses it unchanged (1″ veto, real rows and painted images alike).

Variants, tried in this order, the first that passes K1r and the cost check adopted, no other tried:
- A: L = 80 (10″, the aperture);
- B: L = 40 (5″).

**Known case K1r (computable).** The union mask, with the 1″ veto, covers the macs0417-ncf S_max position
(`screen.s_max_radec` of `results/w1_shear_ncf/summaries/z0.3_macs0417-ncf.json`).
**Cost check (computable).** In every one of the five NCF, the screen with the union mask keeps ≥ 75 % of the rows
the D-097-masked screen used: kept = `screen.n_used_before_veto − screen.n_veto_used`; kept in the new z0.3
summary ≥ 0.75 × kept in `results/w1_shear_ncf/summaries/z0.3_<field>.json` (function `cost_ok`; the z0.3 run is made first for the check; the screen does not
depend on Z).
Neither variant passes both → the route fails, recorded as a failed approach; the headline stays at fourteen fields.

**Run (if a variant is adopted).** All five NCF re-run with the union masks (`inject_shear.py --streak-masks`, the
frozen `results/w1_shear_ncf` grid: Z ∈ {0.1, 0.3, 0.6, 1.0, 1.5}, masses 2 × 10¹², 8 × 10¹², 2 × 10¹³ M☉, 100
lenses per cell, seed 50, R 0.425). **Combined limit:** `w1_volume.py` over the eleven stored cluster summaries
(D-097 masks, unchanged) plus every NCF whose new screen does not pass the E/B rule. An NCF whose new screen passes
stays out of the limit and is looked at on a sheet; it goes to `w1_peak_vet.py` / `/vet-candidate`.

**Outcomes.**
- Both flagged NCF admitted and n₉₅ at 2 × 10¹³ M☉ within 2× of Takahashi & Asada → blank deep fields next
  (`results/w1_shear_ncf/PLAN.md`).
- macs0417-ncf still passes with the veto → the stripe explanation is weakened; vet it (`/vet-candidate`).
- A previously null NCF passes with the veto → flagged and vetted like any other.

**Reported regardless.** Per NCF and variant: stripe segments (rows, columns), mask pixels, share of used rows
vetoed; for the run: S_max, p_rot, max |S_×|, recoveries per cell. Every mask sheet is looked at before any
conclusion.
