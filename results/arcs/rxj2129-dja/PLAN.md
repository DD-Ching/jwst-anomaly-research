# Plan: `arcs` on RXJ2129 with the DJA v7.0 catalogue (frozen before any `arcs` count on it; 2026-10-16)

**Why.** The RXJ2129 run on the program 2767 F200W pipeline catalogue (#283, `results/arcs/rxj2129/`) gave **no
result**: 12 strong-shear background rows, 8 aligned (67 %) under both models, so both controls failed. The DJA
grizli v7.0 photometric catalogue of the same field (`rxj2129-grizli-v7.0-fix_phot_apcorr.fits`, already in the
photo-z tarball used by #283) has 158 rows passing e ≥ 0.5, a ≥ 2 px, S/N ≥ 10 within 50″ (counted while choosing
this step; no shear, model or class was computed), against 52 selected rows from the pipeline catalogue.

**Not independent.** Same field and same two models as #283, whose post hoc anti rows (447, 464) were looked at.
A result here is a re-test of RXJ2129 with a different shape source, never a fifth independent field. Rows that
match #283's anti rows within 0.3″ are listed.

**Shape source (changes from #283).** SEP moments on the DJA `ir` detection image (HST WFC3/IR + JWST NIRCam
combined, 0.04″ px, north-up): `a_image`, `b_image`, `theta_image`. The PSF is broader than F200W alone, which
rounds small sources; this can only lower the aligned fraction (ASSUMPTION).

**Conversion to the pipeline columns (`scripts/arcs_rxj2129_dja.py`, unit-tested):**
- `sky_orientation` = PA east of north of the major axis. For a north-up, east-left image, PA = θ + 90° (mod 180)
  with θ = `theta_image` in degrees counter-clockwise from +x. **Validated before the run** on sources matched
  (≤ 0.2″) to the #283 pipeline catalogue with e ≥ 0.5 in both and S/N ≥ 20: the median |ΔPA| must be ≤ 10°,
  else the run stops and the convention is recorded as failed.
- `ellipticity` = 1 − b/a; `semimajor_sigma` = a × (0.04 / 0.031) so the D-024 rule a ≥ 2 px keeps its
  F200W-pixel meaning (0.062″); `isophotal_flux(_err)` = `flux_iso`, `fluxerr_iso`; `isophotal_area` = `area_iso`;
  `isophotal_abmag` = 28.9 − 2.5 log₁₀ `flux_iso` (catalogue ZP); `label` = `id`.
- `is_extended` = True for every row (no star flag in the catalogue; ASSUMPTION: the e ≥ 0.5 cut removes most
  stars; anti rows are looked at by eye).

**Everything else as #283 (frozen there):** D-024 defaults (ASSUMPTIONs), the same Lenstool (primary) and GLAFIC
(second) maps by zip range reads and `validate`, the same DJA photo-z, and the same outcome rule
`arcs_macs0647.outcome(primary, second)`: the primary decides if its control passes (≥ 70 % of strong-shear
background rows aligned), else the second, else **no result**; the deciding run **replicates** the S1063 pattern if
P(≥ n_anti) < 0.05 under the shape-noise expectation, else **does not replicate**.

**What changes the plan:** replication sends the anti rows to the S1063 vetting chain (ICL subtraction, measurement
noise, clump-pair census) before any reading, and the tally records it as a non-independent re-test; non-replication
or no result is recorded as such for RXJ2129.

**Data:** the 24 MB photo-z tarball (sha256 `d5b19aa8…ef83`; catalogue member sha256 `6278025d…a7fe`), the 1 MB
pipeline catalogue (validation only) and the map members by range reads, all in scratch space, deleted after use.
