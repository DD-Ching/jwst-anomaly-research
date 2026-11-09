# `arcs`: single-component anti rows after vetting (post hoc, exploratory)

Plan: `PLAN.md` (committed before the computation). Script: `scripts/arcs_anti_ledger.py`; class selection, row
dropping and the outcome rule unit-tested (`tests/test_arcs_anti_ledger.py`). Every number is in `summary.json`.
Input: the stored `../pooled-census/census.ecsv` only; no new data.

## Frozen outcome: **single-component excess not supported after vetting (post hoc)**
| Variant | Rows | Anti | Expected | P(≥ anti) |
|---|---|---|---|---|
| reproduce (addendum 9 post hoc split) | 54 | 3 | 0.631 | 0.0242 (stored value reproduced) |
| **primary**: drop 2194 and 581 | 52 | 1 | 0.572 | **0.440** |
| drop 2194 only (581 an arc, Lenstool reading) | 53 | 2 | 0.582 | 0.114 |
| obs_only: drop the two rows, keep all probabilities | 54 | 1 | 0.631 | 0.472 |

`p_anti_row` is the shape-noise null (σ_ε 0.26; model_prediction). ASSUMPTION threshold 0.05.

## Row ledger (derived)
| Row | p_anti_row | Status |
|---|---|---|
| WHL0137 2194 | 0.0484 | `model_sensitive` (addendum 13) |
| RXC J0600 581 | 0.0107 | `member (model input)` under GLAFIC v3; an arc under Lenstool (addendum 19) |
| MACS0308 2157 | 0.0008 | open: `robust` to Zitrin model error (14); MUSE redshift `undetermined` (16, 17) |

## Reading (narrowest statement)
- The post hoc single-component excess (P 0.024) does not survive the later vetting: with 2194 and 581 out, one anti
  row against 0.57 expected (P 0.44). Even keeping 581 as an arc, P is 0.11.
- 2157's own row probability is small (8.4 × 10⁻⁴), but it is one of 52 rows; the class-level tail is the frozen
  test. Looked at `../anti-bands/anti_bands.png`: 2157 is a bright (F200W 20.6) single galaxy with an extended
  envelope, axis stable across bands; nothing beyond an ordinary galaxy shape is visible.
- Together with addendum 9 (pooled variant B explained by clump orientation, P 0.484), the `arcs` replication fields
  hold no anti excess. 2157 stays an open single row, not a candidate.

## Limitations
- Post hoc: the class split was defined after the counts were seen.
- One-sided selection: only anti rows were vetted for membership or model sensitivity; aligned rows that are members
  would also leave the class and lower the expected count. The `obs_only` variant bounds this in the other direction
  only for the two dropped rows.
- 581's status depends on the model (GLAFIC member, Lenstool arc); the `drop_2194_only` row covers the Lenstool reading.
