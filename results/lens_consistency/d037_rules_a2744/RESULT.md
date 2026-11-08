# D-037 rules as functions: second-field gate on Abell 2744 (result)

Plan: `PLAN.md` (frozen before the run). Script: `scripts/d037_rules_a2744.py`. Outputs: `summary.json` (frozen
route), `summary_post_hoc.json` (post hoc alias). Inputs: the stored 2026-10-14 D-041 run
(`results/lens_consistency/abell2744-uncover2/`) and Bergamini+23 Table 2 from vizier.cfa.harvard.edu
(J/ApJ/952/84/table2, sha256 `e79aa3c7…`, accessed 2026-11-08; not stored).

**Frozen outcome: route check fails; gate A not run as written.** Table 2 has 149 images, but 3 ids of the stored
`imageplane_bergamini23.ecsv` (46.1a–c) are absent: Table 2 names that system `JD1A`–`JD1C`. The six docs-table
positions agree with Table 2 within 0.024″ (limit 0.1″).

**Post hoc (deviation): alias JD1B/JD1A/JD1C → 46.1a/46.1b/46.1c.** The pairing is the one where each predicted
image's stored `sep_image_arcsec` (0.386 / 0.755 / 0.279″) equals its separation from the Table 2 position
(0.38 / 0.76 / 0.28″; session check, not stored). With it the route check passes and **gate A passes 6 / 6**:

| Model | Row | Want | Got | Inputs |
|---|---|---|---|---|
| UNCOVER | 3.1 `absent` | `model_topology` | `model_topology` (3.1a/3.1b) | no second-model image |
| UNCOVER | 34.1 `absent` | `model_topology` | `model_topology` (34.1a/34.1b) | Bergamini 2.74″, μ ratio 1.60 → not μ-dependent |
| Bergamini | 700.1 `absent` ×2 | none | none, none | pair 700.1a/b beyond 8″; UNCOVER image 3.56″ / 12.74″ (> 3.0″: μ rule not applicable) |
| Bergamini | 18.1 `confused` | none | none | no second-model image |
| Bergamini | 22.1 `confused` | none | none | UNCOVER image 28.8″ away |

Unread rows (no recorded hand reading; outside the gate): UNCOVER 3.3 (none; Bergamini 1.52″, μ ratio 1.09),
C200.3 (`model_topology`: C200.3a/C200.3b share one UNCOVER prediction), 301.1 (none). C200.3's label is a new
diagnostic reading, not checked by eye in this run.

**Limits.** The Bergamini model's `shared_match` is the recorded 2026-10-08 list, not a re-run (www.fe.infn.it reset
every connection, 16:12Z), unlike the MACS1149 route. Bergamini 700.1's first row has μ ratio 11.7 at 3.56″, just
outside the 3.0″ radius (ASSUMPTION): the μ label depends on that radius. Only one rule (repeated pair) is exercised
positively here; `mu_model_dependent` gives only `False` / not-applicable on this field, so its positive case rests
on MACS1149 1192 alone. Provenance: forced classes derived; μ model_prediction; positions observed.

**Reading.** Both functions agree with every hand reading on a second field (post hoc alias only); they can be
reported as a `d037_label` diagnostic column. No flag, no candidate; nothing here bears on exotic lensing.
