# Lenstool in the cloud: install and reproduce the port on CANUCS MACS0416 (plan, committed 2026-10-26 09:12Z)

Why: the D115b record (`docs/candidates/macs0416-canucsref-D115b.md`) waits on a Lenstool re-fit with member 8515's
σ free, so far listed as "owner machine only" because Lenstool was not installable here (DECISIONS.md D-024 era note:
"C code that needs compiling, not pip-installable"). conda-forge now lists `lenstool` 8.15.5 (linux-64).

Hypothesis: the conda-forge build installs in a cloud session and reproduces the repository's dPIE port on the
pinned `macs0416-canucs` model. Ordinary failure modes: missing linux-64 build, proxy/TLS blocks, a version whose
`.par` syntax differs from the CANUCS file, or a numerical disagreement from a different cosmology or z handling.

Gate R (frozen): Lenstool's predicted image positions (`runmode` image-plane prediction from the catalogued
sources) for the D115a / D115b families agree with the port's `predicted_macs0416.ecsv` positions within 0.1″
(ASSUMPTION; the port's 0.1″ solver grid). Pass → a restricted re-fit (8515 σ and its neighbours free, cluster halos
fixed) becomes a cloud step. Fail → record the cause; the re-fit stays on the owner machine.

## Amendment (committed 2026-10-26 09:18Z, before any port value of this comparison was computed)
Deviation, post hoc with respect to the Lenstool output: the install (micromamba, conda-forge `lenstool` 8.15.5,
53 s) and one Lenstool run of the pinned `best.par` (its own `runmode image 1`, 7 s) were done first, and `image.all`
was looked at (it predicts, for each catalogued image, the images of that image's own back-traced source). The port's
`predicted_macs0416.ecsv` uses the mean source per family, so it is not a like-for-like reference. Gate R is therefore
re-stated on the same inputs, for **every** catalogued image (not only D115):
- For each catalogued image i: the port back-traces it (`backtrace_images`, the fiducial single-plane model, as
  `best.par`) and solves for all images on the 0.1″ grid (`find_images`); Lenstool's rows for i in `image.all` are
  the reference.
- Per image i, port and Lenstool predictions are matched one-to-one (Hungarian on distance); a pair matches when
  ≤ 0.1″ (ASSUMPTION).
- **Gate R passes** if (a) ≥ 95 % (ASSUMPTION) of all Lenstool-predicted images have a matched port image, and (b)
  for both D115b.1 and D115b.2 both codes predict an image within 0.7″ of the flag (64.048551, −24.064825).
  Both conditions are required. Function `gate_r` in `scripts/lenstool_cloud_check.py`, unit-tested.
- Reported, not gated: port images with no Lenstool match, the median and 95th percentile of matched distances.
