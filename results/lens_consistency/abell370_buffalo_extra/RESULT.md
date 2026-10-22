# Abell 370 BUFFALO model A extra-image test: result

Plan: `PLAN.md` (frozen and pushed before the run). Numbers: `summary.json`, tables in this directory.

- **Gate K** (rule 1): 104 BUFFALO gold images, rms 0.910″, 99.0 % within 3″ → passes (repeats gate K2 of
  `../abell370_buffalo`).
- **Stage 2** (reported, not used): pipeline forced-test gate C 26 / 27 testable targets recovered; no pipeline flag.
- **Gate C2** (rule 4, local-noise and containing-source rules): all controls 35 / 42 = 0.833; arc-referenced controls
  19 / 25 = 0.76 < 0.80 → **fails**. As frozen, the reference route is not validated on Abell 370 and no flag is read
  (the flag list is empty by rule).
- Reported, not rules: without the containing-source rule C2 also fails (arc-referenced 18 / 25 = 0.72); without the
  local-noise rule, arc-referenced 19 / 25.
- Post hoc, not used: the 7 absent controls are all counted under "sources too faint" in `posthoc_absent_controls`;
  cause not tested.
- **Deviation:** the CANUCS files were downloaded through the MAST download API (`mast:HLSP/canucs/dr1/a370/clu/...`)
  and the BUFFALO maps through the same API, because archive.stsci.edu reset connections from the cloud session; the
  bytes are verified against the sha256 pinned before the run. The script's URLs name archive.stsci.edu.
- Outcome: as on MACS1423 (C2 failed on arc-referenced controls), the CANUCS-reference route does not pass on Abell 370.
