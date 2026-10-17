# W1 shear screen: what sets Abell 2744's B-mode floor? (diagnostic, 2026-10-17)

**Why.** TASKS "JWST focus" 1 (W1): injection recovery needs S > the field's max |S_×| (D-053), and Abell 2744's
4.21 at (88, −42) is the highest of the five fields. The claim was a streak veto (S1063's top peaks sit on streaks);
the sheet below showed no streak here, so the run became a diagnostic (**deviation** from the claimed scope; no plan
was frozen and nothing is adopted).

Labels: every number `derived` from the committed inputs (F200W GLIMPSE catalogue sha256 `9953c869…`, DJA v7.2 zout
`43662686…`, Bergamini+23 model); all tolerances are ASSUMPTIONs. Script: `scripts/w1_floor_diag.py`.

## Results (`diag_abell2744.json`, `null_vs_bmax_abell2744.json`)
- **Streaks: not seen.** `sheet_abell2744_bmode.png` (F200W, 28″ box, S3 byte-range read; looked at): the aperture
  holds a bright elliptical, three edge-on bright galaxies and a compact bright galaxy; the top contributors are
  small sources inside or next to their light. No long linear feature crosses the box.
- **The peak is carried by many sources:** 44 in the aperture; the top 1 / 3 / 5 give 0.12 / 0.31 / 0.47 of S_×.
  Aperture mean ε 0.044 against a field mean |ε| 0.331, so a uniform PSF-anisotropy offset is not indicated (a
  pattern that varies across the aperture is not tested).
- **It exceeds shape noise (indicative):** 1 of 200 rotation-null maxima of S reaches 4.21 (99th percentile 4.02).
  This compares a two-sided max |S_×| with the one-sided max S null, so p ≈ 0.01 is indicative only.
- **Collinear chains (streak proxy) are common by chance at these tolerances:** 21 % of used sources have ≥ 2
  aligned, collinear neighbours within 20″ (10°/10°); 50 % in the peak aperture. Not used as a veto.
- **Bright-neighbour rule (source within k semimajor σ of a neighbour ≥ Δm brighter), screen only, no injections:**

  | Δm, k | used (of 4,951) | b_max | S_max | share of peak S_× from flagged rows | flagged share of field |
  |---|---|---|---|---|---|
  | 2, 2 | 4,333 | 4.19 | 3.39 | 0.00 | 0.12 |
  | 2, 3 | 3,325 | 3.86 | 3.40 | 0.21 | 0.33 |
  | 2, 5 | 1,270 | 3.89 | 2.97 | 0.81 | 0.74 |
  | 3, 3 | 3,436 | 3.86 | 3.40 | 0.18 | 0.31 |
  | 3, 5 | 1,361 | 3.47 | 2.81 | 0.76 | 0.73 |

  Flagged rows are not over-represented at the peak (share ≤ the field fraction in every row), so this rule does not
  single out the floor's sources; its b_max drops come with a 31–74 % source loss and are not separated from the
  effect of fewer sources per aperture (not tested).

## Conclusion (narrow)
Abell 2744's B-mode floor is not a streak and not one or a few outlier shapes; it is a coherent cross pattern over
the 44 aperture sources around a group of bright galaxies, above the rotation null (indicative). None of the tested vetoes is
adopted. W1 limits are unchanged.

## Next
Test a spatially varying systematic directly: a PSF model from stars across the mosaic (ε_PSF at each source), or a
bright-galaxy light-gradient term (local background slope at each source from the i2d), then injections under a
frozen adoption rule (lower limit at each mass, as in `results/w1_edge_cut/PLAN.md`).
