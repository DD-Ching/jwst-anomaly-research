# MACS0647: DJA NIRSpec redshifts for the Meena+2023 families (plan, frozen before the table is read)

Why: the held-out test (`../macs0647_clash_heldout/RESULT.md`, D-034 addendum w) fitted each family's z, so it could
not test the redshift structure; it named spectroscopic redshifts as the input for a stronger test. DJA
`dja_msaexp_emission_lines_v4.6.csv.gz` (2026-09-15; the newest on the bucket, v4.7–v5.1 HEAD 403 on 2026-10-24) is
the public graded JWST spectroscopic table already used in `../macs1149_s16_specz/dja_v46/`.

- **Hypothesis under test:** at their spectroscopic z, the CLASH LTM-Gauss maps (the usable model of D-034 addendum
  w) reproduce every certain image of the families with a spec-z within 3″.
- **Ordinary explanations for a miss:** model error (HST-era model, no JWST constraints), a wrong spec-z match (a
  neighbour or cluster member within the radius), a misidentified family member.
- **Data:** the v4.6 table, streamed and filtered on the fly (not stored; cloud disk rule): rows within 0.5″
  (ASSUMPTION) of any Meena+2023 image (`../macs0647_inputs/meena23_images.ecsv`, uncertain rows included for the
  match only).
- **Rules (functions in `scripts/macs0647_specz.py`, tests before the run):**
  - S1 spec-z per image: DJA `grade` 3 (ASSUMPTION, as before); several grade-3 rows → their median if they agree
    within 0.01 (1 + z), else the image is `ambiguous` and gets no z.
  - S2 spec-z per family: the median of its images' S1 z if all agree within 0.01 (1 + z) (ASSUMPTION); else the
    family is `inconsistent` (recorded, not used: a mismatch needs vetting first, e.g. a neighbour).
  - S3 control gate: the pre-JWST families (1–12) with an S2 z, run at that z, must have ≥ 4 certain images and
    ≥ 80 % of them within 3″ of a predicted image (ASSUMPTIONs, as gate K0). If S3 fails or has < 4 images, the
    held-out reading is `not decided`.
  - S4 reading: a certain image of a held-out family (13–23) with an S2 z that has no predicted image within 3″ at
    that z → `flag` (to `/vet-candidate`); otherwise `no_flag`. Only when S3 passes.
- **Outcomes:** no family with a spec-z → `no input` (the route waits for new spectra); S3 fails → `not decided`;
  S3 passes and no flag → the held-out redshift structure is consistent for the families tested (state which);
  flags → vet.
- Post hoc, not rules: z_fit (D-034 addendum w) against each spec-z.
