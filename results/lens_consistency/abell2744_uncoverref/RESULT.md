# Abell 2744: UNCOVER-reference extra-image test, stage 1 (result)

Plan: `PLAN.md` (frozen and committed before the run). Script: `scripts/a2744_uncoverref.py`. Numbers:
`summary.json` (`stage1_gate_K`), per-image table `imageplane_abell2744_abell2744-uncover2-maps.ecsv`, image list
`uncover2_images.ecsv` (`observed`, the release's table with `z_spec` from `Redshift_type`). Image-plane positions are
`model_prediction`; offsets are `derived`.

**Frozen outcome: gate K fails; the UNCOVER v2.0 maps are not used for the extra-image stages (PLAN.md rule 1).**
On the 172 spectroscopic images, image-plane rms 4.90″ (gate ≤ 1.0″, ASSUMPTION), 94.2 % within 3.0″ (gate ≥ 80 %).
Run time 1 min 47 s.

**Post hoc (not part of the frozen rule).**
- *What drives the rms.* Three families get a single predicted image for a catalogued triple: 34.1 (a/b 32.2″ /
  31.0″), 80.1 (24.8″ / 22.6″), 74 (19.6″ / 18.7″); then 4.1e 9.7″. Without the 10 images not reproduced within 3″
  (`drop_unreproduced`) rms is 2.22″ on 162 (still > 1.0″, refitted 80.4a lands 25.3″ off); the median offset of all
  172 is 0.66″, and the rms over the images within 3″ is 0.93″.
- *Likely cause (ordinary, untested here).* `lensmodel.imageplane_residuals` takes the family source as the **mean**
  back-traced position; `a2744_uncover.py` (D-041) uses the magnification-weighted least-squares source and reproduced
  all 181 own-list images at rms 0.84″. With a mean source placed outside the caustic, the grid search returns one
  image, which matches the single-image families above. The plan listed the estimator as an ordinary explanation; it
  was not swapped after the run, because that would be a choice made after seeing the outcome.
- `extra_predicted` lists 33 families (e.g. system 2's seven knots, 4 catalogued vs 5 predicted); not read further,
  since the gate failed.

**Reference catalogue (rule 3, availability only, no flux read).** From the cloud session, Google Drive links on
`https://jwst-uncover.github.io/DR3.html` (accessed 2026-10-23) answer with byte ranges: DR3 LW SUPER photometric
catalogue `UNCOVER_DR3_LW_SUPER_catalog.fits` (83.7 MB, 74,020 rows, 156 columns incl. `f_f200w`, `e_f200w`,
`a_image`, `b_image`, `theta_J2000`, `use_phot`, `flag_nearbcg`, `z_spec`, `z_phot`), the LW detection segmentation map
(7.0 MB gzip) and the bCG-subtracted F200W mosaic `uncover_v5.2.0_abell2744clu_f200w_bcgs_out_sci.fits.gz` (1.84 GB
gzip, so no byte-range cutout; a stage-2 local-noise rule would need the whole file or another noise source). Not
pinned (no sha256 yet).

**Next.** This model is stopped for the extra-image test under the frozen estimator. Options, each a new frozen plan:
(a) gate K with the WLS source estimator of D-041, frozen as a deviation (like BUFFALO K2) and labelled as chosen
after this failure; (b) gate K with Bergamini+23 (`abell2744-bergamini23`, Lenstool) on its own spectroscopic list.

## Stage 1b: gate K with the WLS source (deviation after the stage-1 failure)

Plan: `PLAN_1b.md` (frozen and committed before the run; rule 1 committed as the tested function
`a2744_uncoverref.gate_k_wls` before the run). Numbers: `summary.json` (`stage1b_gate_K_wls`), per-image table
`imageplane_wls_abell2744_abell2744-uncover2-maps.ecsv`. Run time 29 s (not stored).

**Frozen outcome: gate K-WLS passes.** On the same 172 spectroscopic images with the same 0.25″ grid, swapping only
the family source (mean → magnification-weighted least squares, D-041): rms 0.835″ (≤ 1.0″, ASSUMPTION), 100 %
within 3.0″ (≥ 80 %). **This is a deviation chosen after stage 1 failed**, and it validates the maps on their own
constraints only (the images are the model's inputs), as BUFFALO K2 (D-034 addendum h). Model sha256 as stage 1.

**Reported (rule 2, not rules).** No image beyond 3″ (largest offset 2.31″, median 0.60″; from the per-image table).
`fewer_predicted`: none, so the three single-image families of stage 1 (34.1, 80.1, 74) are reproduced, which
supports the estimator reading of stage 1's RESULT (no other setting changed). 24 images share their nearest
predicted image with another image of the family (`shared_match`). 31 families have one more predicted image than
catalogued (`extra_predicted`; 3 vs 2 or 5 vs 4). These are the inputs of the extra-image stages; they are not
read here and are not anomalies by themselves (an unlisted counter-image can be faint, demagnified or blended).

**Next.** Freeze stage 2 with the WLS source: predicted / control tables, the UNCOVER DR3 SUPER catalogue as reference
(pin by sha256, freeze the column mapping), gate C2 with the frozen shared-pair rule (D-034 addendum i), then flags.
The 1.84 GB bCG-subtracted mosaic cannot be cut by byte range, so the local-noise rule needs another noise source
(e.g. catalogue errors) or a stated download reason.

## Stage 2a: position-angle convention of the DR3 SUPER reference

Plan: `PLAN_2a.md` (frozen and committed with its tested rule functions before the run). Numbers: `summary.json`
(`stage2a_theta_convention`, with the 500 segment ids used). Run time 36 s (not stored).

**Frozen outcome: H_a adopted.** `theta_J2000` is radians counter-clockwise from +x (SEP's `theta`): median residual
2.7° (≤ 10°), 98.4 % within 20° (≥ 70 %). The README says "degrees clockwise from the first image axis (along y)"; read as clockwise from +y (H_b) it fails:
46.6°, 23.2 %. **Post hoc** (review): read as clockwise from +x it also fails (43.4°, 24.0 %;
`posthoc_cw_from_x`). Sample: 2,188
rows pass rule 1 before the segment cut; the 500 with the largest `a_image` and ≥ 50 segment pixels are used.

**Reported.** Catalogue pixel grid: 0.0400″ per pixel, +y at PA 0.006°, +x at PA −90.0° (north-up, east left). So sky PA
= deg(θ) + 90 (mod 180). A contact sheet of six of the 500 segments (ids 54653, 51026, 26828, 25022, 37572, 33881) with both
readings drawn was looked at in the session (image not stored): the H_a line follows each segment's long axis.
**Post hoc** (review): segment centroid vs catalogue `x`, `y` (1-based pixels assumed), median 4.9 px, 90th
percentile 17.4 px (`posthoc_centroid_offset_pix`; per source in `stage2a_sources.ecsv`); id = segment label is not
tested as a rule, but a mismatch would scatter the H_a residuals, which have median 2.7°.

**Limits.** Validates the column against the release's own segmap (unweighted segment moments), not against an
independent image. No minimum sample size was frozen (here n = 500). The full segmap is read into memory (fine in this
session; a cutout reader would be lighter). The README's flux definition (`f_f200w` total, 10 nJy) is taken as stated, not tested here.

**Next.** Stage 2b, its own frozen plan: predicted / control tables with the WLS source over the GO 2561 tiles; gate
C2 with the shared-pair rule; local noise from GO 2561 byte-range cutouts (decision in `PLAN_2a.md`).

## Stage 2b: extra-image test with the WLS source, DR3 SUPER reference and local noise

Plan: `PLAN_2b.md` (frozen and pushed with its tested functions before the run). Numbers: `summary.json` `stage2b`;
tables `predicted_abell2744.ecsv`, `control_abell2744.ecsv` (frozen), `without_contain/`, `frame/`. Run time 21 min.

- **Frame (rule 1, reported, not used):** 172 spectroscopic images, WLS source; image classes 161 observed, 20
  candidate, 6 missing, 10 demagnified, 5 no_flux_ref, 1 faint; 9 images unpredicted. Pipeline forced test on GO 2561
  `o001_t003` F200W: gate C 60 / 60 testable recovered or confused; no pipeline flag.
- **Gate C2 (rule 5, shared pairs) passes:** all controls 30 / 34 = 0.882, arc-referenced 11 / 13 = 0.846 (≥ 8, ≥ 0.80,
  ≥ 3). Five shared-pair groups merged (sizes 4, 4, 4, 3, 2). Plain `gate_c2` also passes (arc-referenced 18 / 20).
- **Flags (rule 6): none.** Target classes after rule 4: 4 recovered, 1 confused, 2 ambiguous, 11
  undetectable_local, 2 undetectable, 11 no_reference, 1 inconsistent_reference.
- **Reported, not rules:** without the local-noise rule C2 passes and three targets would be `absent`: 3.2
  (|μ| 48.2; F_pred / σ_loc 2.2 → undetectable_local), 4.2 (|μ| 9.4; 9.7 σ → ambiguous) and 78.3 (|μ| 6.7; 8.7 σ →
  ambiguous). Without the containing-source rule C2 fails (arc-referenced 7 / 13).
- **Contact sheet looked at** (F200W, 6″ stamps, 1″ search circle; image not stored): 3.2 and 4.2 lie in the smooth
  BCG / ICL light (the mosaic is not bCG-subtracted, so σ_loc is large there, as PLAN_2a expected); 4.2 is the D-041
  4.2c position. 78.3 lies in a clean field with no obvious source inside the circle except a faint knot near its
  edge. Not interpreted: 78.3 is system 78 (= 700.1, `../a2744_s700_vet`), already undetermined there.

**Outcome.** The DR3 SUPER reference with pipeline-mosaic local noise passes gate C2 on Abell 2744 (exploratory: the
local-noise and shared-pair rules were set on earlier fields, not tuned here). No frozen flag. The two `ambiguous`
targets sit just under the 10 σ absent threshold, so the result depends on that ASSUMPTION.

**Limits.** Single band (F200W); pipeline mosaic not bCG-subtracted (completeness near members lost); F_pred is a total
flux against an aperture σ (lenient for extended images); no bCG-subtracted re-check of 3.2 / 4.2.

**Next.** A frozen follow-up of 78.3 and 4.2 on a bCG-subtracted image or a second band (F277W / F356W, same
rules), or `/vet-candidate` only if a frozen rule flags them; not a re-read of this stage's thresholds.

## Stage 2c: per-band re-test in F277W and F356W

Plan: `PLAN_2c.md` (frozen and pushed with `outcome_2c` and its unit test before the run). Numbers: `summary.json`
`stage2c`; tables `f277w/`, `f356w/` (`predicted_` / `control_abell2744.ecsv`). Frame: the committed stage-2b tables.

- **Gate C2 (shared pairs) passes in both bands:** F277W all 28 / 32, arc-referenced 11 / 13; F356W all 36 / 42,
  arc-referenced 14 / 16.
- **Frozen flags (rule 5):** F277W `4.2@3.579980,-30.408545`; F356W the same 4.2 position and
  `78.3@3.575198,-30.355756` (78.3c below). **Outcome (rule 6):** `flag` for 4.2 and 78.3c; `present` for 26.1,
  26.2 (×2) and 67.2; every other target `undetermined` (mostly `undetectable_local` / `no_reference`).
- 4.2: `absent` in F277W (F_pred / σ_loc 15.8) and F356W (12.5); stage 2b F200W 9.7 (`ambiguous`).
- 78.3c: F277W 9.4 σ (`ambiguous`), F356W 13.2 σ (`absent`); stage 2b F200W 8.7 σ. Brightest catalogue source in the
  1″ circle: 9.6 % of F_pred in F356W (`max_ratio`).
- 3.2 (both positions) stays `undetectable_local` in both bands.
- **Reported, not rules:** without the local-noise rule both bands flag 3.2, 4.2 and 78.3c.

**Contact sheet** (`sheet_2c.png`: F200W / F277W / F356W, 6″ stamps, 1″ circles; made by a session script, not
stored; looked at):
- 4.2 lies on the smooth halo of a bright member with no source in the circle; its reference 26191 is a compact knot
  with extended structure.
- 78.3c lies in a clean field; a compact source sits inside the circle near its edge (catalogue 51679, 0.73″,
  F356W 9.6 % of F_pred).
- The reference 54653 of system 78 is a long thin arc (`a_image` 115.9 px = 4.6″, `b_image` 10.2 px) next to a
  bright member.

**Ordinary explanations already on record; vetting not finished:**
- **4.2 = 4.2c, vetted in D-036** (`docs/candidates/abell2744-family4-c.md`: μ systematics of the family 4 c images
  next to member 34423). In this run the sibling control 4.1c (same reference 26191) is also `absent` in all three
  bands (F_pred / σ_loc 11.2 / 11.6 / 16.5; `control_abell2744.ecsv`), so the deficit is shared by the catalogued c
  image. Not a new flag.
- **78.3c is not yet vetted.** Its reference 54653 is the arc that is also the counterpart of 78.1a and 78.2a, so
  F_pred scales a shared arc. Yet the same reference predicts the controls 78.1c and 78.2c, about 1.5″ away (separation not stored), which are
  `recovered` at 10.7 / 10.9 σ in F356W. Known ordinary routes from `../a2744_s700_vet` (system 78 = 700.1): a
  member at 2.1″ and a group halo at 4.4″ make μ there depend on the local mass; 78.3 has a single catalogued image
  (78.3a), so its source rests on one image. Exotic readings are not considered.

**Outcome.** Exploratory, as stage 2b. A second band turns the stage-2b `ambiguous` 4.2 and 78.3c into frozen
flags. 4.2 repeats the vetted 4.2c. 78.3c is an unvetted flag, not a candidate.

**Next.** `/vet-candidate` on 78.3c, cheapest test first:
1. model dependence (D-036 rule: untestable if μ moves > 2× or an independent model disagrees). Bergamini+23
   reproduces 700.1 at no z (`../a2744_s700_zscan`, outcome B), so its members cannot be perturbed for this system as
   they are. The cheapest independent check is whether the CATS v4.1 maps (D-035) predict 78.3c from 78.3a;
   then μ(78.3c) under ±30 % of the nearby member and group components;
2. 78.3a's knot identity (a single-image system);
3. the flux of the 51679 source.
