# Band-restricted MIRI confirmation (follow-up to `RESULT.md`)

Status: plan, written before measurement (2026-11-07, about 01:08Z; see the claim commit).

- **Question.** `RESULT.md` adopted NIRCam confirmation for the 163 single-band-only MIRI sources as a group
  (61 matched, chance 0.051), but the F1500W/F1800W match rates (0.071, 0.086) sit near the group's chance rate,
  and the per-band chance rate was not computed. Which MIRI detection bands can be confirmed on their own?
- **Hypothesis.** F770W and F1000W pass a per-band gate; F1500W and F1800W fail it.
- **Ordinary explanations for a match.** Chance coincidence with the dense DJA catalog (per band, from shifted
  positions); the larger PSF at long wavelengths (0.59″ FWHM at F1800W) shifting centroids beyond 0.3″, which
  would lower the real rate there.
- **Rule (ASSUMPTION, frozen here).** Same radius (0.3″), S/N floor (5), shifts (10, 5–10″) and gate as `PLAN.md`,
  applied per detection band of the single-band-only group: band b is confirmable when n_b ≥ 15 and
  chance_b < 0.2 × real_b. No other band selection, radius change or post hoc split.
- **What changes the plan.** If every band passes, the group rule stands as is. If F1500W/F1800W fail, the
  pipeline confirmation (when wired) is restricted to the passing bands. If F770W fails, the group adoption in
  `RESULT.md` is withdrawn.
