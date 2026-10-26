# Lenstool in the cloud: install and reproduce the port on CANUCS MACS0416 (plan, 2026-10-26 09:25Z)

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
