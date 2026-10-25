# Held-out residuals across fields and models: synthesis (result)

Plan frozen in f9c48518, rule functions and tests in the next commit, both before the run
(`scripts/heldout_synthesis.py`, `tests/test_heldout_synthesis.py`). Numbers: `summary.json` (frozen rules),
`posthoc.json` (post hoc). Inputs: the 11 stored held-out tables of D-034 addenda ap, at, au, aw, ax, ay; no new data.
Labels: offsets are model_prediction residuals of observed positions; thresholds (3.0″ miss, 1.0″ fit, 95 %) are
ASSUMPTIONs.

## Frozen outcome
- **Flags (rule 4): none.** No held-out image is both missed by two primary models of its field and isolated in
  its family.
- **Shared misses (rule 3):** only SMACS 0723, images 6.1, 7.1, 17.1 (glafic v2 and zitrin-ltm-gauss v2). None is
  isolated: the other images of families 6, 7, 17 sit at 1.1–4.9″ in both models (`heldout_*.ecsv` of addendum aw),
  so these are family-wide model offsets, not single-image anomalies.
- **Isolated misses (rule 2): none** in any of the 11 tables.
- **Per-field floor (rule 5)**, 95th-percentile held-out offset of the lowest-rms primary model:

| Field | Model | n | rms | median | P95 = floor |
|---|---|---|---|---|---|
| MACS0416 | MARS | 77 | 0.514″ | 0.261″ | 1.31″ |
| MACS1149 | Sharon v4cor | 17 | 1.12″ | 0.865″ | 2.40″ |
| SMACS 0723 | glafic v2 | 44 | 2.23″ | 1.98″ | 3.39″ |
| Abell 2744 | Bergamini+23 | 43 | 1.78″ | 0.606″ | 4.07″ |
| Abell 370 | MARS | 40 | 14.81″ | 0.165″ | 5.76″ |
| El Gordo | Caminha+23 | 162 | 8.98″ | 0.938″ | 7.36″ |

  El Gordo's n covers all held-out images (geometric-z families included, as the plan's inputs are the full table);
  its spec-z-only gate value (1.64″ on 46) is in addendum ay. Abell 370 has one primary model, so its floor comes
  from a model that fails gate H.

## Post hoc (not in the frozen plan)
- Rule 2 only counted held-out family members. With the model's *seen* images of the family added as members (their
  stored fit residuals, where a stage stored them: Abell 2744 / Bergamini+23, El Gordo / Caminha+23, MACS1149 /
  CATS v4.1), El Gordo 6c (3.05″), 13c (6.26″) and 20c (5.10″) become isolated misses; Abell 2744 and MACS1149 give
  none. These are the JWST-proposed third-image candidates already discussed in addendum ay (Caminha's own predicted
  images there have m_pred 28.5–29.2); with one model on El Gordo they cannot be shared misses, so no flag. The
  addendum-ay contact sheet was looked at; this stage made no new images.

## Reading
- Under these rules the stored held-out tests contain no single-image position anomaly. The misses are either
  family-wide (SMACS 6, 7, 17; Abell 2744 family 80) or single-model third-image identifications (El Gordo).
- The floors are the stated limit of the image-position residual route on these fields: with today's public
  models an image offset below about 1.3″ (MACS0416), 2.4″ (MACS1149), 3.4″ (SMACS 0723), 4.1″ (Abell 2744),
  5.8″ (Abell 370) or 7.4″ (El Gordo) does not stand out from published-model error at about 5 % per image.

## Limitations
- One primary model on four of six fields, so rule 3 can only act on MACS1149 and SMACS 0723.
- Held-out type N images anchor their source on all their own images, so their offsets understate a true
  prediction error; type E (seen siblings) is the cleanest and is a minority.
- P95 on 10–17 images (MACS1149) is coarse.
