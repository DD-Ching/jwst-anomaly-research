# RELICS Lenstool v1 amplitude offset: is it a z_ref normalisation? (2026-10-25)

Frozen plan `PLAN.md` (pushed in 909bd8dd before any header was read). Numbers come from `summary.json`
(`h_u`, `h_c`, `posthoc_makers`), written by `scripts/relics_lenstool_zref.py`. Factors are the stored rule-7
medians (derived); header values are observed (FITS primary headers of the published maps, byte-range reads).
No imaging was read.

## Frozen result
**H_U (one RELICS-wide z_ref = 2.09): rejected.** Prediction f_pred = D_LS/D_S(z_lens, 2.09); rule |f_obs −
f_pred| ≤ 0.10 (ASSUMPTION).

| Cluster | z_lens | f_obs | f_pred | diff | z_ref implied | Consistent |
|---|---|---|---|---|---|---|
| RXC J0232 | 0.284 | 1.002 | 0.787 | +0.215 | none (f > max) | no |
| RXC J2211 | 0.397 | 1.017 | 0.710 | +0.307 | none | no |
| MACS J0553 | 0.412 | 0.665 | 0.700 | −0.035 | 1.73 | yes |
| RX J0152 | 0.827 | 1.045 | 0.462 | +0.584 | none | no |
| Abell 3192 | 0.425 | 0.692 | 0.692 | 0.000 | 2.09 (defines z_ref) | yes |

Three of five clusters reject it; the two offset clusters imply different z_ref (1.73 and 2.09).

**H_C (a normalisation keyword in the map headers): not supported.** The `kappa` and `x-arcsec-deflect` primary
headers of all five maps carry 44 cards each and none names a source redshift, D_LS/D_S or a normalisation (the
keyword/comment pattern in the script matched nothing; COMMENT cards are only the FITS-standard text). Offset and
non-offset maps have the same keyword set.

So by the plan, **H_M stays open**: the offset belongs to the models as released (their mass scale relative to the
SLICE redshifts), not to a documented product convention. The released κ and μ(z = 2) maps agree with the deflection
maps read at D_LS/D_S = 1 (gate P medians 1.000 / 1.006 for MACS J0553, 1.000 / 1.000 for Abell 3192, earlier runs'
`validate.json`), so a deflection map written at another source z would have to be shared by the μ maps too.

## Post hoc (not in the plan)
- **Map maker.** The header keywords that differ between maps (`posthoc_makers`, all 10 registered RELICS Lenstool v1
  maps plus SMACS 0723) include `MAPMKR_F/L`. The two offset maps (MACS J0553, Abell 3192) and SMACS 0723 (offset per
  addendum aw, |α| 1.14–1.56 × GLAFIC, no per-family factor) carry the same map maker and DATE 2017-12-03…06. The three
  near-1 maps carry three other makers (RXC J0232 2018-03-06, RX J0152 2018-01-04, RXC J2211 2017-08-28). That
  association is 3 / 3 vs 0 / 3 on clusters chosen before the headers were read, but defined after; it is not a
  cause. Untested maps from the same maker: RXC J0600 and MACS J0417 (2017-12-05 / 06), RXC J0018 (2020-07-06).
- **Prediction for the next test (hypothesis H_K):** if the offset follows that 2017-12 map batch, the RXC J0600 and
  MACS J0417 Lenstool v1 maps need a factor ≈ 0.65–0.72 against z_spec families or a second model, and RXC J0018 /
  SPT0615 / WHL0137 do not. A clean failure of either half rejects H_K.

## Deviation
- PLAN.md said "primary and first extension headers". Every map stores its image in the primary HDU (NAXIS 2,
  `EXTEND` set but no extension was read), so only primary headers were read.

## Reading
- The Abell 3192 z_ref ≈ 2 coincidence (addendum bh) does not generalise: no single normalisation fits the five
  clusters, and the maps document none.
- The addendum-bd per-cluster amplitude check stays the only rule before reading a RELICS Lenstool v1 map; no global
  correction is applied.
