# Lenstool in the cloud, and the port vs Lenstool on CANUCS MACS0416 (result, 2026-10-26)

Plan `PLAN.md` (committed 09:12Z; gate R amended in a commit at 09:18Z, tool committed 09:19Z, before any port value was computed, see "Deviations"). Tool
`scripts/lenstool_cloud_check.py`, rules tested in `tests/test_lenstool_cloud_check.py`, both committed before the
run. Numbers: `result.json`; Lenstool's output `image.all.gz`. Every image position is a `model_prediction` of the
pinned `macs0416-canucs` `best.par` (D-044); the catalogued images are `observed`.

## Result
- **Install works.** `micromamba create -p <scratch>/lt -c conda-forge lenstool` installs Lenstool 8.15.5 (conda-forge
  build `py314h03d985e_1`, with cfitsio 4.6.4, gsl 2.8, wcslib 8.2.2; build and dependency versions not stored)
  in 53 s (not stored) in a cloud session. No
  compiler, no host outside the allowlist; `conda.anaconda.org` and `micro.mamba.pm` answer.
- **Lenstool runs the pinned model unchanged:** `lenstool best.par -n` (its own `runmode image 1`, 4 OpenMP
  threads, not stored) writes `image.all` in 7 s (not stored).
- **Gate R passes** (`gate_R` true). All 303 catalogued images appear in `image.all`; all 912 Lenstool-predicted
  images have a port image within 0.1″ (one-to-one): median 0.00093″, 95th percentile 0.0052″, maximum 0.033″.
- **The D115b extra image is Lenstool's too:** for D115b.1 and D115b.2, both codes predict an image within 0.7″ of
  the #458 flag. So the flag is not an artefact of the port.
- **Reported, not gated:** the port finds 67 images with no Lenstool counterpart (`n_port_unmatched`); among them
  D115a.1, D115a.2 and D115b.1 (6 port vs 5 Lenstool images each, `per_image`; D115a.3 1 vs 1, D115b.2 5 vs 5). Not diagnosed (post hoc reading,
  untested: Lenstool's `image.all` may omit faint or central images).

## Reading
The D115b re-fit that the vetting record (`docs/candidates/macs0416-canucsref-D115b.md`) left to the owner machine can
now run in a cloud session, on a Lenstool that reproduces the port's image positions to ≤ 0.033″ on this model.

## Limits
One model, one runmode (image-plane prediction from each image's back-traced source); Lenstool's optimiser
(`runmode inverse`) is not tested here. The port was compared single-plane (as `best.par`); the D-046 two-plane
option has no Lenstool counterpart in this check.

## Deviations
- Gate R was re-stated after the Lenstool run and a look at `image.all`, before the port comparison was computed
  (PLAN "Amendment"): the original reference, `predicted_macs0416.ecsv`, uses the mean source per family and is not
  like-for-like with `image.all`.
- The first invocation was started from `scripts/` and failed at the JSON write (relative output path); it was re-run
  from the repository root.
- `/code-review` fixes after the run (outcome unchanged, `result.json` byte-identical on the re-run): the output path is
  anchored to the repository; matching prices pairs beyond 0.1″ out of the assignment (so far pairs count as unmatched);
  gate R also fails when a catalogued image is missing from `image.all` or has no finite back-trace (0 here); the
  worker pool is spawn-safe (initializer).

## Next
A frozen plan for the D115b Lenstool re-fit: `runmode inverse` on `input.par` with member 8515's σ (and the nearest
members') free and the rest fixed at `best.par`, image-plane χ², then read Δχ² and the flag image with the D115b
rules. Validate first: the same restricted run with nothing freed must return the `best.par` χ².
