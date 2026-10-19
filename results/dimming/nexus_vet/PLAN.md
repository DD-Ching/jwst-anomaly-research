# NEXUS W3 dimming survivors: vetting plan (frozen before any read, 2026-10-19 11:20Z)

Input: the 5 artefact-test survivors and 10 flagged controls of `../nexus_calib/` (#369).

**Hypothesis under test:** none of the 5 survivors is a W3 (wormhole-lens) dimming, so the three-field headline in
`../nexus_calib/headline.json` holds as written.

**Ordinary explanations, tested in order (docs/methodology.md):**
1. Single-band single-epoch points (d53823, d60156, d79928): cosmic ray / snowball / detector residual. Test: is the
   point present in each of the individual `_cal` exposures (or the per-exposure DQ / CR flags) of that epoch? A
   real astrophysical point appears in every dither at the same sky position; a CR appears in one.
2. d73390: extended-nucleus aperture / PSF-mismatch variability. Test: aperture-size dependence of the F200W
   variation (0.15″ vs 0.3″ vs 0.5″) and the same measurement on nearby galaxies of similar brightness.
3. d79956: ordinary transient. Test: TNS / SIMBAD / NED at the position; light-curve shape (rise-only = not a dip).
4. 10 flagged controls: are they the same kind of single-point residual (check 1 on a sample)?

**Decision rule (ASSUMPTION):** a survivor is "explained" if one of the tests above gives a specific ordinary cause;
"unexplained" otherwise. Any unexplained survivor that shows a dimming in both bands goes to `/vet-candidate` in full
and the headline is quoted as conditional. What changes the plan: an unexplained two-band dimming.
