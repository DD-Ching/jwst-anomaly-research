# MACS1423 batch-2 box-2 h000: HST archive light curve (plan, fixed before any HST flux is measured)

Owner instruction 2026-10-18 (after #332). Budget 09:28Z → **11:28Z** (session clock), records included.

## Target (locked)
- **B2-MACSJ1423+2404-box2-h000 (215.949099, +24.086948)**; not the batch-3/4 h000 entries in other fields.
- Known (`../batch_venus2/boxes/review/lightcurves.ecsv`, 0.15″, host light included): JWST 1208 (MJD 60001.6,
  2023-02-26) point source in 8 bands, e.g. F150W 551.9 ± 4.4 nJy; VENUS 6882 (MJD 61104.7, 2026-03-05) F150W
  275.1 ± 8.4, no point source. Host candidate GLASS 00487 (z 0.525, 0.72″).

## HST products (metadata, MAST query 09:35Z; `hst_obs.ecsv`)
- Pre-detection: ACS F435W and F606W, 16667, 2022-02-05/12 (10 visits) and 2022-08-26/30 (5 visits); CLASH 12790 and
  others 2006–2014 (many bands).
- Post-detection: ACS F435W 16667 2023-08-25 (+180 d); 17504 WFC3/IR F110W + WFC3/UVIS F606W 2023-12-12/14 (4 visits,
  +289–291 d); 17504 F110W + UVIS F814W 2024-04-10/11 (2 visits, +409 d).
- Visits within 7 days form one epoch; dithers and visits inside an epoch are not independent epochs (their scatter
  is reported as a consistency check).

## Measurement (reused `hit_lightcurve.aperture_photometry`, the #323 HST convention)
- Pipeline association drizzled products (`_drc` for ACS/UVIS, `_drz` for IR), read by S3 byte range (box ±8″);
  e⁻/s × PHOTFNU → nJy. Aperture 0.2″, annulus 0.6–1.0″, no aperture correction.
- **Local error:** robust std of the same aperture at 300 random positions in the ±8″ box, ≥ 1.5″ from the target
  (stated per visit; checked against the visit-to-visit scatter of the target within an epoch).
- **Host:** constant per instrument + filter; only differences between epochs of the **same instrument and filter**
  are interpreted (F435W ACS 2022-02 → 2022-08 → 2023-08; F606W ACS 2022-02 → 2022-08; F110W WFC3/IR CLASH 2013 →
  2023-12 → 2024-04). Cross-instrument pairs (F606W ACS → UVIS, F814W ACS 2006 → UVIS 2024) are reported with a
  zero point from field sources only and labelled limited.
- **Field-source zero point:** compact field sources in the box (local maxima with S/N ≥ 20 in the reference
  epoch, ≥ 1″ from the target, de-duplicated at 0.5″, at most 30): median and robust scatter of their flux ratio per
  epoch pair; ≥ 8 needed, else the pair is reported without a zero-point correction and labelled.
- **Rule:** Δ = F(epoch) − m · F(reference), σ = sqrt(err_epoch² + (m·err_ref)² + (s·F_ref)²) with m, s the
  field-source median and scatter. |Δ|/σ ≥ 5 ⇒ change detected; otherwise a 3σ limit on |Δ|. Stamps are looked at.
- Interpretation, fixed now (no type is forced): a pre-detection change (2022) means the source was variable or
  present before 2022-08; none means it was not there at the measured level. A post-detection non-detection is a
  limit on the decline; a detection is a late-time point (rate of decline). Colour from one filter pair at one epoch
  only. None of this decides lensing vs intrinsic.
- Cap: if the S3 reads stall, a 60 s socket timeout and one retry; stop measurement at 10:45Z.
