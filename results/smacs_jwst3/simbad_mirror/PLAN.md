# SMACS 0723 jwst 3.0.0: SIMBAD star rule and lens flags via the CfA SIMBAD mirror — plan (before the run)

Follow-up to `../galaxy_stratum/RESULT.md` (#743), which separated stars with Gaia only because CDS
(simbad.cds.unistra.fr, cdsxmatch.u-strasbg.fr) resets every connection from the cloud session.

- **Observation (2026-11-06 17:15Z):** the CfA SIMBAD mirror `simbad.cfa.harvard.edu` answers SIMBAD TAP
  (`SELECT ... FROM basic`, 200, SMACS-field objects returned); CDS still resets.
- **Tool change:** a `simbad` backend `mirror` (SIMBAD TAP basic + alltypes on the CfA mirror), chain
  `xmatch → tap → mirror`. Validated first on the recorded `tap` fixture targets: the mirror must return the same
  main_ids / otypes for them (known-case check) before it is used on SMACS.
- **Hypothesis:** with the SIMBAD rule (a nearer SIMBAD galaxy wins over a Gaia position-only match) and
  `is_lens_related` flags, the six surviving records stay in the galaxy top 50 and keep ranks within ±3 of the
  Gaia-only run (#743); the four `gaia_position`-only stars that ranked 9–21 post hoc are settled one way or the other.
- **Ordinary explanations to keep in mind:** compact cluster members labelled stars; SIMBAD lens-related otypes
  (arcs, multiple images) flagging catalogued lensed images, which is a known-population explanation, not a hit.
- **What changes the plan:** the mirror disagreeing with CDS-recorded rows on the fixture (then it is not used), or
  a record leaving the galaxy top 50 (then the earlier verdicts are re-opened in the RESULT).
