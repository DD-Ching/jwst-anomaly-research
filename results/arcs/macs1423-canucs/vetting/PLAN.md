# Vetting the MACS1423 `arcs` anti rows, step 1: ICL-subtracted shapes and measurement noise (plan, 2026-10-16)

Written and committed before either script runs. Parent: `../RESULT.md` (frozen outcome "replicates" at screen level).

- **Code:** the unmodified S1063 vetting scripts of #273 at commit `43fbbbd` (`scripts/arcs_bgsub_shapes.py`,
  `scripts/arcs_vet_noise.py`), run with `--model macs1423-canucs-maps --arcs-dir results/arcs/macs1423-canucs`
  on `s3://stpubdata/jwst/public/jw01208/L3/t/o010/jw01208-o010_t005_nircam_clear-f200w_i2d.fits` (byte ranges).
  The two files are copied unchanged next to this branch's `scripts/` (which holds the model entry; #273's
  `lens_consistency.py` predates it). Every parameter is the script's default (S1063 values), except
  `--targets 1818 2197 2322 4234`.
- **Hypothesis to falsify (ordinary):** the anti rows have true axes aligned with the prediction, and their anti
  axes come from ICL gradients (bgsub step) or photon/ICL noise in the shape measurement (noise step).
- **Frozen rules (the scripts' own functions, unit-tested on #273):**
  - bgsub known-case check: as in #273 (pipeline-aligned rows stay aligned under the re-measured axis);
  - noise gate `gate_ok`: bright (m 24) injections come back aligned in ≥ 95 % for every row, else the noise test
    is invalid and supports nothing;
  - per target `noise_class`: `noise-explicable` if P_anti ≥ 0.05, else `not noise-explicable`;
  - sample `excess_p`: Poisson-binomial P(≥ n_obs) with the rows' P_anti, n_obs = rows `anti` under
    `bgsub_shapes.ecsv`; `excess explained by measurement noise` if ≥ 0.05, else `excess remains`.
- **What changes the plan:** sample excess explained → the MACS1423 replication is closed as a measurement-noise
  effect, and the same calibration becomes the next test of S1063 itself; excess remains → the remaining S1063
  chain (deblend, clump-pair census, second band) runs on the surviving rows next.
