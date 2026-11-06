# SMACS 0723 jwst 3.0.0 with SIMBAD (CfA mirror): the star rule changes nothing; the six records keep their ranks

Plan: `PLAN.md` (pushed before the run). Run `20261106T171253Z-28af4751`, `configs/smacs_jwst3_simbad.yaml`
(as `smacs_jwst3_gaia_only.yaml` plus SIMBAD in `classify` and `crossmatch`). Tables: `candidates.csv`,
`summary.json` (`scripts/smacs_jwst3_compare.py --out`).

## Tool check first (observed)
- CDS still resets connections (SIMBAD TAP, XMatch); `simbad.harvard.edu/simbad/sim-tap` answers.
- Known case: on the five crossmatch fixture targets, the mirror returns the same 3 SIMBAD objects (main_id, otype,
  all types) as the CDS recording of 2026-10-07 (row order differs). The run used backend `simbad: mirror`
  (crossmatch meta), after `xmatch` and `tap` failed.

## Result (derived)
- Quality-passing star basis: 36 SIMBAD/Gaia astrometry, 97 stellar locus, 6 Gaia position only, the same counts as
  the Gaia-only run (#743). The 2938 ranked galaxies, and every record's rank and score, are identical to #743
  (`candidates.csv` equal column for column): 2804 6, 829 7, 1032 8, 438 11, 1243 15, 2915 27; 1096 551, 2559 687;
  2925 no counterpart. **Hypothesis (records within ±3 of #743): holds, with Δ = 0 for all six.**
- **Why nothing moves:** none of the six Gaia-position-only stars (f200w_3141, 3209, 3469, 3794, 3999; f277w_816)
  has any SIMBAD object within 2″ (mirror query, observed), so the "nearer SIMBAD galaxy wins" rule cannot apply.
  Post hoc ordinary check: their five Gaia DR3 sources (G 20.6–21.0) all have
  `classprob_dsc_combmod_star` ≥ 0.9993 and are not in `in_galaxy_candidates` (observed catalogue values; DSC is
  less reliable at G ≳ 20, so this supports, not proves, "star"). #743's bound on the galaxy top 25 (at most these four
  entering at ranks 9–21) is therefore no longer needed: they stay stars under both available rules.
- **Lens flags on the galaxy top 20:** SIMBAD adds no lens-related type. NED `G_Lens` flags the records 829
  (f277w_743) and 1032 (f200w_921), as before. New SIMBAD identifications: 1243 = `[RBI2023] 18` (G, 0.06″) and
  f200w_3774 = `[VBG2023] SMACS 511` (G, 0.02″): catalogued galaxies, a known-population label only.
- Contact sheet seen: identical top 20 to #743; f200w_3177 (#10) still sits in a bright star's spike and halo
  (likely artefact, not screened by D-019); f356w_1098 has no F200W flux (F356W detection), as recorded.

## Reading
The SIMBAD-dependent part of the SMACS jwst 3.0.0 re-run is closed: no recorded verdict changes, the Gaia-only
star separation was not biasing the galaxy stratum, and no new lens-related identification enters the top 20.
Nothing here is a candidate; none of the new top-25 entries is vetted.

## Limitations
- The mirror's content date is not stated by the service; it matched CDS on the 3-object fixture only.
- Gaia DSC probabilities are a post hoc check (ASSUMPTION: ≥ 0.999 counts as stellar support).
