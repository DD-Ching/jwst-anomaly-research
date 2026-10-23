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
