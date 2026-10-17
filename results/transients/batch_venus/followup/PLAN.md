# h002 identity, h003 limited re-check, wording (plan, written before any new query)

Owner instruction 2026-10-17 (after #296). Budget 90 min: 07:11Z → **08:41Z**. h002 has priority. No new sky
regions, no batch re-run, no threshold changes, no new solver.

## Already done (not repeated)
- **#296:** search, sheet, forced photometry in every covering mosaic (`boxes/review/lightcurves.ecsv`) and the first
  vetting record.
- **#299** (open, another session; branch `claude/h002-host-fix`):
  - `[FPP2024] NS_0481` is **3.77″** from h002. The #296 record read NED's arcmin column as arcsec.
  - z = 1.3236 ± 0.0002 is spectroscopic (NIRSpec MSA, DDT 4446; Frye et al. 2024 Table 2, arXiv:2309.07326).
  - SN H0pe's images are 25–38″ from h002 (literature check, not stored).
  - TNS: HTTP 403, not run.
  - This unit builds on #299 and does not edit its branch.

## h002
- **Known facts (observed):**
  - Forced-photometry detections at MJD 60033 (PEARLS 1176, 8 bands, S/N 14–41) and MJD 60056 (DDT 4446,
    F090W / F150W / F200W).
  - No detection at MJD 60815 (4744) or MJD 61125 (VENUS).
  - The F150W and F444W fading peaks are 0.03″ apart.
  - No Gaia or SIMBAD source within 3″.
- **To verify:**
  1. Reported transients: TNS (retry once; if blocked, "could not query"); the G165, PEARLS and SN H0pe papers and
     their tables (Frye et al. 2024; Pierel et al. 2024; Chen et al. 2024; Pascale et al. 2025; lens-model papers);
     public ATel / AstroNote / GCN text search where reachable.
  2. SN H0pe image positions (2a, 2b, 2c) with their source, compared with h002 by position. h002 is detected in
     both 2023 epochs while H0pe's images are seen in the same epochs, so the dates cannot separate them; position
     decides.
  3. NS_0481: trace the source table, redshift type and object type. With 3.77″ separation it is a neighbour, not a
     confirmed host; whether any galaxy lies under h002 is checked on the existing stamps.
  4. Light-curve evidence: what the multi-band, multi-epoch fluxes and upper limits support (rise or decline, colour,
     timescale, absolute magnitude only if a redshift were known), and what is missing. No SN template fitting.
- **Decision rules (fixed now):**
  - **Matches a reported transient:** a published or TNS position within 0.3″, or within the paper's stated
    uncertainty, at a consistent epoch.
  - **Matches only an ordinary catalogue object:** a catalogue object within 0.3″ that is not a transient.
  - **Identity not found:** no match from the sources that could be queried. Each source not reachable is listed as
    "could not query".
  - **Physical type confirmed:** only with a spectrum or published classification of the transient itself. A light
    curve alone gives "consistent with" wording.

## h003 (MACS0416 box 2, 64.031622 −24.075255)
1. From the committed `lightcurves.ecsv`, identify which ≥ 3σ differences are in the same band and the same aperture
   (0.15″, fixed sky position) between programmes, and whether the two search bands' S peaks coincide on the sky
   (#296: 0.49″ apart).
2. One fixed check with existing tools, if it fits the time. Programme-to-programme photometric scatter of comparable
   sources: `hit_lightcurve_controls.py pick` (arc or extended-light controls in the same box, same flux range) plus
   `compare` on the band pair with the largest |z|. The control rule is fixed before measuring:
   - h003 is "explained by programme-to-program scatter" if its control tail fraction is ≥ 0.05, as in #286 h043;
   - otherwise it stays undetermined.
   - **Fixed parameters (added 07:20Z, before the controls run; step 1 done from the committed table):**
     `pick-arc` on `jw01176-o211_t009_nircam_clear-f356w_i2d.fits` at the h003 position, `--half-arcsec 16
     --match 1.5 --step-arcsec 0.2 --min-sep-arcsec 0.6 --n 80 --seed 1` (the Abell 2744 h182 settings); controls and
     h003 measured together by `hit_lightcurve.py` (0.15″); the decisive pair is the largest-|z| same-band pair,
     F356W `jw01176-o211` → `jw01208-o004` (z = −9.09). Other pairs are reported, not decisive. Fewer than 20 controls
     measured in both mosaics of the decisive pair ⇒ undetermined.
3. If the controls cannot run (too few controls, tool limits), it stays undetermined and the minimal missing input
   is recorded.

## Wording
- Summaries, the candidate record and the state files are aligned to the owner's evidence levels.
- The review queue's remaining count (0) is separated from the unresolved science questions (h002 identity and type,
  h003).
- Original records are kept; corrections are added as dated sections.
