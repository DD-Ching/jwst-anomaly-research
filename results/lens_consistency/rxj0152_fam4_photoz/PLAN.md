# RX J0152 family 4: RELICS photo-z of the three chance-level matches (plan, frozen before any catalogue row is read)

Context: D-034 addendum bg (`../rxj0152_fam4_photometry_v2/RESULT.md`) found three colour-and-flux matches
(segments 112, 127, 131) to family 4 (SLICE model z 1.69, photometric) near the predicted third image, at chance
level (λ 2.27, P 0.40), and named photo-z of those segments as the reopen input. DJA v7 has no RX J0152 field
(checked 2026-10-26 in the public v7 index); the RELICS HST BPZ catalogue does
(`hlsp_relics_hst_acs-wfc3ir_cl0152m13_multi_v1_cat.txt`, HTTP 206 on a header read).

**Hypotheses.** H_img: one of 112 / 127 / 131 is family 4's third image (the models are right; no residual).
H_field: all three are unrelated field or cluster galaxies, leaving the predicted image unidentified (an input to
H_miss, the lens-model residual; not a flag, because gate V of bg failed).
**Ordinary explanations of a mismatch:** BPZ catastrophic outliers (Balmer / Lyman break confusion), HST blending at
JWST-resolved positions, astrometric offset between RELICS and VENUS frames, faint sources missing from RELICS.

## Method (frozen)
1. Read the RELICS catalogue (column names from its `# N name` header, as `scripts/arcs_spt0615.py::read_bpz`).
   Keep rows with `stel` < 0.8. z interval = [`zbmin`, `zbmax`] (the catalogue's 95 % interval; ASSUMPTION that
   it is the right width for a consistency test). Reliability = `odds` ≥ 0.9 (ASSUMPTION).
2. Match radius 0.6″ (ASSUMPTION; RELICS ACS/WFC3-IR resolution), nearest row, no offset correction. The median
   offset of the matched control images is recorded as a diagnostic only.
3. **Gate P (photo-z control):** the 18 SLICE images of families 1–6 (`../slice_relics_rxj0152/slice_images.ecsv`;
   family 1 spectroscopic z 3.93, others SLICE model z). Passes when ≥ 8 images match and ≥ 60 % of the matched ones
   have their family z inside [zbmin, zbmax] (ASSUMPTION thresholds).
4. Per segment (112, 127, 131; positions from `../rxj0152_fam4_photometry_v2/result.json`):
   - `consistent`: matched and 1.69 ∈ [zbmin, zbmax];
   - `excluded`: matched, 1.69 ∉ [zbmin, zbmax] and odds ≥ 0.9;
   - `undetermined`: unmatched, a star, or excluded with odds < 0.9.
   4.1 / 4.2's own RELICS rows are reported (derived) beside the SLICE z.

## Outcome (frozen)
- Gate P fails → **U**: no reading.
- Otherwise **C** if any segment is `consistent` (H_img stays open for that segment; next step a vet / colour
  comparison with 4.1 at matched resolution); **E** if all three are `excluded` (H_field: the predicted third image
  is not among the region's colour-matched sources; still no flag, because bg's injection gate failed); **U**
  otherwise.

**What changes the plan.** E → the third-image question narrows to a depth test at the predicted positions (a new
frozen injection gate). C → the line closes as model-consistent unless the vet fails. U → paused as now.
Labels: catalogue values `observed` (RELICS), photo-z `derived` (RELICS BPZ), family z `model_prediction` /
spectroscopic, classes `derived`.
