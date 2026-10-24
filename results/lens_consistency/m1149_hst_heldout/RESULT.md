# HST-era MACS1149 models vs the Schuldt et al. 2024 new images: both fail gate H as written; neither miss is shared (D-034 addendum au)

Plan frozen in 2ab445c2; script and tests committed before the run. Numbers from `summary.json` (per-image rows in
`heldout_macs1149-cats.ecsv`, `heldout_macs1149-sharon.ecsv`) unless marked post hoc. Labels: Schuldt+2024 and model
catalogue positions and z are observed (as the authors list them); residuals are model predictions; thresholds are
ASSUMPTIONs (`../mars_hff_heldout/PLAN.md`).

## Frozen result
| Model | Seen list | Frame offset O | Seen / held-out | Types E / N / M | Outside map | Gate H rms | Within 3″ | H |
|---|---|---|---|---|---|---|---|---|
| CATS v4.1 | `arcs.txt` (145 rows) | −0.033″, 0″ (7 pairs, applied) | 7 / 11 | 0 / 11 / 0 | 0 | 2.46″ (11) | 82 % | **fail** (rms) |
| Sharon v4cor | gold list (23 rows) | not applied (0 pairs) | 0 / 18 (SXa untested) | 0 / 17 / 0 | 0 | 1.12″ (17) | 100 % | **fail** (rms) |

- **CATS v4.1 already contains families 201 and 301 and SXa** (within 0.06″ / 0.23″ after the frame offset, `summary.json` `seen_rows`; Schuldt's "new" is relative to
  Grillo et al. 2016, not to CATS), so CATS is tested on 202, 203, 204 and 401 only. Its fail is carried by family
  401 (z 5.98): 401.1 at 4.50″ and 401.3 at 6.00″; the other 9 images lie within 2.12″.
- **Sharon v4cor** reproduces every image within 3″ (max 2.40″, 301.3) and fails only on rms (1.12″ vs ≤ 1.0″).
  Its 200″ maps cover families 301 and 401, which fall outside the MARS map.

## Post hoc: same images across the four models (`posthoc.json`, from the stored tables)
| Images | CATS | Sharon | MARS | CANUCS v1 (saw them) |
|---|---|---|---|---|
| 202–204 (7) | 0.88″, max 1.18″ | 1.23″, max 2.08″ | 2.73″, max 3.72″ | 0.76″, max 1.33″ |
| 201–204 (MARS's 10) | — (201 seen) | 1.10″, max 2.08″ | 2.51″, max 3.72″ | 0.72″, max 1.33″ |
| 401 (4) | 3.92″, max 6.00″ | 0.86″, max 1.13″ | outside | 1.24″, max 1.77″ |

## Reading (narrow)
- No image fails in more than one held-out model: the 401 miss is CATS-only (Sharon and CANUCS put 401 within
  1.8″), and the MARS misses on 202 / 203 are absent from CATS and Sharon (≤ 1.45″). This favours model-specific
  incompleteness over a family, z or MUSE-position problem, as addendum at's post hoc CANUCS check suggested.
- On the same images MARS's rms is 2.3× Sharon's (201–204: 2.51″ vs 1.10″) and 3.1× CATS's (202–204: 2.73″ vs
  0.88″), so the MARS miss is mainly MARS-specific; but neither HST-era parametric model passes gate H as written either.
- Neither CATS nor Sharon is a held-out-validated second model on MACS1149 under the frozen rule; no two-model
  extra-image stage there. Sharon's near-pass (100 % within 3″) is not promoted post hoc.
- No flag (frozen: no flag reading). No image pixels read; no contact sheet (catalogue positions only).

## Limits
- All families type N (mean source includes the held-out images), a weak test; 401 is the one z ≈ 6 family with
  4 images, so CATS's rms rests on one family.
- Sharon's seen list is the 23-image gold table transcribed from its readme (Table 5); if Sharon's fit used images
  not in that table, a Schuldt image it saw would count as held out here.
- Rule 1 cannot measure a Sharon frame offset (0 shared images); a sub-arcsecond Sharon-to-HST offset would enter
  its 1.12″ rms directly (ASSUMPTION: none).
