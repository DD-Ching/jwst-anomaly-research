# Clump-pair census on five other `arcs` fields (exploratory)

Plan: `PLAN.md` (committed before the run, 468f018). Script: `scripts/arcs_census_fields.py` (S1063 census code
unchanged; outcome rule `plan_outcome`, unit-tested) → `census_fields.{json,ecsv}`, `two_knot.png`. Labels: pixels
observed; classes derived; injections simulations; f̂ and P model_prediction. Depends on #273 (S1063 census).

## Census (frozen rule, own-flux gate per row; `census_fields.json` "per_field")
| Field | rows | gate passed | `two-knot` | f̂ | 95 % upper |
|---|---|---|---|---|---|
| Abell 2744 | 35 | 33 | 8 | 0.242 | 0.423 |
| Abell 370 | 56 | 34 | 9 | 0.265 | 0.444 |
| El Gordo | 12 | 12 | 2 | 0.167 | 0.484 |
| MACS0416 | 27 | 19 | 7 | 0.368 | 0.616 |
| MACS1149 | 24 | 19 | 6 | 0.316 | 0.566 |
| **Pooled** | 154 | 117 | 32 | **0.274** | **0.364** |
| S1063 (#273, for comparison) | 24 | 21 | 5 | 0.238 | 0.472 |

No row fell outside its mosaic. Every field has ≥ 12 gate-passed rows.

## S1063 under the five-field rate (frozen outcome; `census_fields.json` "s1063_p_ge_4_anti")
- Mixture null with the S1063 per-row p_cp, p_sn (`abells1063-cats/vetting/census.ecsv`): P(≥ 4 anti) =
  **0.031** at the pooled f̂ and 0.054 at its 95 % upper.
- **Frozen outcome: "S1063 excess remains at the five-field clump rate (exploratory)"** (0.031 < 0.05).
- Reading (narrow): an independent clump rate (0.27) is close to S1063's own (0.24), so it does not remove the
  excess; at the pooled upper (0.36) P is just above 0.05. The S1063-only upper (0.47) gave 0.093 (#273); the pooled
  rate narrows that range.

## Looked at (`two_knot.png`, all 35 `two-knot` rows incl. 3 gate-failed)
- Most are two compact knots or a knot inside or beside a larger body, as in S1063.
- Several are a compact source plus a neighbouring object, not one clumpy galaxy: by eye (not stored) Abell 370
  3979, 4401 and 4614 (a knot on a long streak or arc), and Abell 370 3846 (compact core in a halo). MACS0416 2411,
  2580 and 2980 show no clear pair at this stretch. Such rows raise f̂, so the pooled rate, if anything, overstates
  the clump-pair fraction; the outcome direction (excess remains) does not depend on them.

## Post hoc (not in the plan)
- The other fields' screen-anti rows: of the 7, 4 passed the gate and none of those is `two-knot` (Abell 370 340
  and 4540, El Gordo 742: `single`; MACS1149 4236 `two-knot` but gate-failed) (`census_fields.ecsv`). In S1063, 2
  of the 4 anti rows (1316, 1557) are `two-knot`.

## Limits
- Same point-pair ASSUMPTIONs as the S1063 census (isotropic separations, point clumps, ≤ 0.4″ pairs).
- The five fields differ from S1063 in depth, filter set of the catalogues, lens models and photo-z sources; the
  pooled f̂ assumes one clump rate for all of them.
- Injection seeds are seed + row index within each field (as frozen), so rows at the same index in different fields
  share a random stream on different images (code review; not re-run).
- 2004 (κ ≈ 1) is still counted in k = 4 with the mixture probability.

## Status
- 1316, 1557: open, unvetted, not candidates. The anti excess in S1063 stays at P ≈ 0.03 (exploratory) under the
  clump-aware null with an independently measured clump rate.
- Next: spectroscopic redshifts (owner decision on MUSE: wait), or the arcs test on a new field with a deep F200W
  mosaic to see whether S1063's anti excess recurs (an independent replication of the excess, not of f).
