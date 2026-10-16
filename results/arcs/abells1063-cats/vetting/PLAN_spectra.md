# Plan: can public spectroscopy settle 1316 and 1557? (owner request 2026-10-16; ≤ 1 h; written before any cutout is measured)

**Question.** Does existing public spectroscopy give a realistic chance of a redshift for 1316 or 1557? This is a
feasibility check, not a new pipeline.

**Order.** (1) catalogued redshifts or 1-D spectra; (2) server-side cutouts of public MUSE cubes (ESO SODA, each
< 200 MB); (3) whole cubes only as a cost report for the owner.

**Line windows (MUSE 4750–9350 Å; photo-z 16–84 % range from `../arcs.ecsv`, widened by ±0.05):**
- 1557, z 1.27–1.43: [O II] 3727/3729 at 8460–9060 Å (the only strong line in range); [Ne III] 3869 is weak.
- 1316, z 1.87–2.36: C III] 1909 at 5480–6410 Å; C IV 1549 at 4750–5210 Å only for z ≳ 2.07; [O II] and Lyα are
  outside MUSE. A line elsewhere would mean the photo-z is wrong (Lyα for z ≳ 2.9; [O II], Hβ, [O III] or Hα at
  low z).

**Measurement (fixed here).** Circular aperture of radius 0.5 × PSF FWHM (ObsCore `s_resolution`) and at least
0.4″ around the catalogue position. Sum DATA; variance = sum of STAT. A 6-channel (7.5 Å) running boxcar, S/N per
window position. Outside the target windows, the same statistic gives the empirical S/N spread (sky residuals); the
noise is rescaled by the robust (MAD) S/N spread of the whole spectrum if that exceeds 1.
- **Line detection:** boxcar S/N ≥ 5 after rescaling, inside a target window, and not on a channel where the
  rescaled S/N spread is driven by a sky line (|S/N| ≥ 5 at the same wavelength in a 2″-offset sky aperture).
- **Limit:** otherwise 5 × the rescaled boxcar σ (erg s⁻¹ cm⁻² per 7.5 Å), median over the target window.
- **Expected flux (model_prediction, ASSUMPTIONs):** continuum from the F090W isophotal AB magnitude (pipeline
  catalogue, `observed`), rest EW 30 / 100 / 300 Å.

**What changes the reading (stated before measuring).**
- Line at a low redshift (z ≤ 0.4, i.e. a member or foreground) → the source is not lensed, so its orientation is
  random; this supports an ordinary explanation of the anti axis.
- Line in the photo-z window → the background redshift is confirmed. Only the photo-z failure explanation is
  removed; the clump-orientation and shape explanations stay open. It does not support exotic lensing.
- No line → undetermined. It supports nothing, neither ordinary nor exotic.
- MUSE seeing (0.6–1.5″) cannot separate the 0.15–0.21″ clumps; no per-clump redshift is attempted.
