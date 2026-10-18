# W1 shear screen at lower masses: filter pilot result (2026-10-18)

Plan: `PLAN.md` (frozen 11:20Z before any injection). Runs: `inject_shear.py --lens-z 0.3 --mass 2e12 1e12 5e11
--n-inject 100 --seed 50`, no streak masks, three filter variants on `ceers-o001` and `goodsn-o004`
(`summaries/<variant>_<field>.json`). Labels: injected lenses `simulated`; recoveries `derived`; R = 0.414 and every
threshold an ASSUMPTION. About 80 s per job, six jobs in parallel on four cores.

## Recoveries (out of 100 per cell; θ_E at z_s = 2: 3.63″ / 2.57″ / 1.82″)
| Variant | field | 2 × 10¹² | 1 × 10¹² | 5 × 10¹¹ | N (1e12 + 2e12) |
|---|---|---|---|---|---|
| V0 Schirmer 10″ (frozen screen) | ceers-o001 | 11 | 4 | 1 | **15** |
| V0 | goodsn-o004 | 18 | 3 | 0 | **21** |
| V1 point-mass 10″ | ceers-o001 | 2 | 3 | 1 | 5 |
| V1 | goodsn-o004 | 8 | 1 | 0 | 9 |
| V2 Schirmer 5″ | ceers-o001 | 8 | 6 | 1 | 14 |
| V2 | goodsn-o004 | 11 | 1 | 0 | 12 |

**Frozen outcome: V0 stays.** No variant exceeds V0's N by ≥ 50 % in both fields (V1 is 0.33× / 0.43×, V2 0.93× /
0.57×). The stop condition does not fire: V2 recovers 6 / 100 at 1 × 10¹² M☉ in ceers-o001 (≥ 5). Under the plan the
lower-mass step is a V0 grid at 1 × 10¹² M☉ over the forty fields.

**Reproduction check.** The V0 real-data screens reproduce the committed Z = 0.3 summaries: S_max 3.9934
(ceers-o001) and 4.3458 (goodsn-o004) equal `results/w1_shear_ceers/summaries/z0.3_ceers-o001.json` and
`results/w1_shear_goodsn/summaries/z0.3_goodsn-o004.json` (that pointing's streak veto removed 0 rows). Recoveries
at 2 × 10¹² M☉ differ (ceers-o001 11 here vs 14 committed) because the pilot paints without the streak masks.

**New looks (recorded, not vetted, in no limit).** V1 and V2 screens of the real data: no E/B rule pass in either
field (V1 S_max 3.84 / 3.36, p 0.015 / 0.48; V2 3.44 / 3.69, p 0.155 / 0.055). V0 goodsn-o004 passes the rule as it
did in the GOODS-N run (the mosaic-edge flag explained there by V1 ∪ V2, `results/w1_shear_goodsn/RESULT.md`).

## Reading
- Re-weighting the aperture does not open lower masses: the point-mass filter loses most of the 2 × 10¹² M☉
  recoveries, and the 5″ Schirmer aperture trades 2 × 10¹² for at most a few 1 × 10¹² recoveries. With N ≈ 5–21
  (Poisson ±2–5) only V1's loss is clearly outside the counting noise.
- Hypothesis (not tested): the inner, heavily weighted radii are where a negative point mass leaves no images
  (inside θ_E) and where few resolved sources fall, so the 1/θ² weight amplifies shape noise more than signal.
- **Expected value of the next step (derived, rough):** V0 recovers 3–4 / 100 at 1 × 10¹² M☉ against 11–18 at
  2 × 10¹², so a forty-field 1 × 10¹² grid would give a limit roughly 3–5× weaker than the 2 × 10¹² one (9.32× Takahashi
  & Asada), i.e. about 30–45× Takahashi & Asada. It would close the 1–2 × 10¹² gap of the comparison's mass range
  (|M| > 10¹²) but not reach it. T4 (cost per gain) argues for lower priority than other W1 work.

## Limitations
- Two pointings, one lens redshift (Z = 0.3), 100 lenses per cell; the selection used baseline efficiency.
  goodsn-o004 is a flagged pointing (out of the GOODS-N limit); its recovery threshold sits above a mosaic-edge
  artefact peak, so its efficiencies are a lower bound for a clean pointing of its depth.
- Three variants only; a point-mass filter with a larger inner radius (r_min ≈ θ_E) was not tried.
- No streak masks (applied to every variant alike).
