# MACS0416 arc-difference hit o211/o213 h108 — vetting record

## 2026-10-15 (cloud routine 01VE5srt)
- **source_uid / run.** `h108` of `arc_diff.py` o211→o213 (PEARLS jw01176 t009, 20″ box at 64.0377, −24.0687;
  `--kernel-px 4 --kernel-stamps`; commit 97ef7ab, #242). Position (pipeline frame) 64.03682, −24.06626.
  Hit: F200W + F277W, max |S| 5.7, sign +, not near a bright core (`model_prediction`/`derived`; thresholds ASSUMPTIONs).
- **Hypothesis before measuring:** a point source on the faint thin arc that brightened between o211 and o213 (a
  caustic-crossing star would be the ordinary population). It would show as one PSF-sized spot, with the same sign in
  every band, and its fluxes would add up across the four epochs. A subtraction residual of the arc would instead
  appear on the neighbouring knots of the same arc too, and its sign could differ between bands.
  Plan change: a point-like, additive, band-consistent excess → full /vet-candidate with a lens-model check;
  otherwise close it as an artifact.
- **Epochs (observed, primary header DATE-BEG):** o211 2022-10-07 (PA_APER 293°), o212 2022-12-29 (33°), CANUCS
  jw01208 o004 2023-01-11 (48.5°; 8 exposures), o213 2023-02-10 (71°). There is no jw01176 o214 (the #242 handoff
  meant the third epoch, o213). CANUCS o004 is an independent fourth epoch between o212 and o213.
- **Five pairs** (8″ box centred on h108, same settings), forced aperture fluxes (nJy, r = 0.08″ SW / 0.12″ F277W,
  F356W / 0.15″ F444W; d = later − earlier; ± = robust spread over apertures with similar epoch-1 arc flux; derived):

  | Pair | F150W | F200W | F277W | F356W | F444W | arc-coherent (top band) |
  |---|---|---|---|---|---|---|
  | o211→o212 | −1.4 ± 2.1 | +1.6 ± 1.9 | +5.6 ± 2.1 | +9.9 ± 1.8 | +6.3 ± 2.4 | no (F356W 0.33) |
  | o211→c004 | +1.9 ± 1.8 | +7.5 ± 1.7 | +10.7 ± 1.7 | +13.4 ± 1.5 | +13.7 ± 2.1 | yes (F356W 0.67) |
  | o211→o213 | +4.1 ± 2.2 | +4.9 ± 1.9 | +12.8 ± 2.1 | +5.0 ± 1.9 | +11.1 ± 2.6 | yes (F277W 1.00) |
  | o212→o213 | +5.5 ± 2.1 | +2.4 ± 1.9 | +7.1 ± 2.1 | −4.3 ± 1.8 | +4.1 ± 2.7 | yes (F356W 1.00) |
  | c004→o213 | +2.1 ± 2.0 | −2.7 ± 1.8 | −2.3 ± 1.8 | −4.5 ± 1.6 | −2.9 ± 2.5 | yes (F356W 1.00) |

  The fluxes add up across pairs within the errors (for example F356W: o211→o212 + o212→o213 = +5.6 vs direct +5.0).
  Additivity holds for any per-epoch flux difference, so it does not separate a real source from an epoch-level
  residual. The largest change (o211→c004, ≈ 13 nJy ≈ 28.6 AB in F277W–F444W) would be a typical caustic-star amplitude.
- **Arc coherence (the deciding test):** the 3–4 knots of the same arc 0.25–0.7″ from h108 carry residuals of the
  same sign and similar S in the same band and pair. Examples: c004→o213 F356W −5.5, −5.3, −4.2 (h108 −3.7);
  o212→o213 F277W +4.2, +6.6, +4.8 with F356W −4.1, −3.8, −5.2 (F277W rises while F356W falls); o211→o213 F277W
  +7.5, +7.7 (h108 +8.0). A single star cannot vary coherently over ≳ 0.7″ of arc with opposite signs in adjacent
  bands. This is a few-per-cent, band-dependent residual of the arc segment's extended light (PSF or background
  matching; the epochs' PAs differ by 15–138°).
- **Validation of the coherence test on known transients (derived):** the published Yan et al. 2023 arc transients
  in their discovery pair (o211→o212 for D21, o211→o213 for D31) are **never** flagged arc-coherent (0 of 9
  testable; D23-S3 and Dc2-W4 have no discovery pair here). Flags occur only in the other pair, where the transient
  is weak (|S| 2–4) and the Warhol arc carries other transients. Caveat: on an arc with several real transients,
  coherence can veto a weak one.

| Test | Tier | Outcome | Evidence | Label |
|---|---|---|---|---|
| image edge / low exposure | A | does not explain | WHT and exposure count at h108 normal in every epoch (4, 4, 8, 4 exposures) | observed |
| bright-core rosette | A | does not explain | `near_bright_core` false; the bright galaxy ≈ 1.2″ away leaves a rosette in the S maps that does not reach h108 (contact sheet) | derived |
| cosmic ray / snowball | A | does not explain | ≥ 3 bands and ≥ 3 epochs involved | derived |
| arc-segment matching residual | A | **explains** | arc coherence in 4 of 5 pairs, including both pairs without o211; band-dependent signs | derived |
| SIMBAD / NED / Gaia cross-match | B/C | not run (moot after the tier-A explanation) | — | — |
| point transient (caustic star) | C | not supported | no excess above the arc residual of the neighbouring knots; cannot exclude a ≲ 10 nJy point source hidden in it | derived |

- **Verdict:** `artifact: arc-segment-coherent PSF/background matching residual` (not a point transient).
  Not a candidate. Unexplained-candidate count from MACS0416 arc differencing stays 0.
- **Regenerate:** `arc_diff.py` with the bands of `results/transients/macs0416_arcdiff/vet_h108/summary_<pair>_h108box.json`
  ("command" field; `--maps` to a scratch path), then `scripts/arc_diff_vet.py --maps <pair> <npz> ... --ra 64.03682099747297
  --dec -24.0662625071866 --out .../h108.json`. Yan transients: same, on the 20″ box maps, one JSON per transient.
- **Open questions:** the arc-segment residual itself (why a few per cent, band-dependent) is not removed; a local
  kernel (per arc segment) might. Contact sheets were looked at (1.5″ and 0.6″ stamps, five pairs, five bands).
