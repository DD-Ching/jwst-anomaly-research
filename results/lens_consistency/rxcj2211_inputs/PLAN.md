# RXCJ2211-0350 input check (plan, 2026-10-24)

Next step of the lens-consistency route (TASKS: "another field with public JWST-era maps"). A Zenodo search
(2026-10-24) found a VENUS GLAFIC v0.91 model of RXCJ2211-0350 (Zenodo 10.5281/zenodo.20312327, CC0, 43 MB zip:
kappa, gamma1, gamma2 maps and `full_output_v91.csv` with source- and image-plane outputs).

Input check only, as for MACS J1931 (addendum aa): list the files, read the README and headers, parse the image list,
and record which of (image list, model maps or predictions, reference catalogue) exist. No residual is read and no
prediction is compared with data in this step.

## Frozen before any comparison (2026-10-24, 13:25Z)

Already seen when this section was written (input listing only): the zip (md5 matches the record), the README
(zl 0.397, χ² 174.99 / 180, rms 0.61″; spec-z references for systems 1, 14, 32 only), the CSV header and row counts
(146 rows, 43 systems, 10 rows with z −1), and the kappa header (2000 × 2000, 0.1″ / px, RA---CAR, ZS = 3) with its
min / max / median. The release has **no deflection maps**: only κ, γ1, γ2 at z_s = 3, plus per-image deflections and
back-traced sources in `full_output_v91.csv`. No deflection or source value has been read.

Checks (script `scripts/rxcj2211_inputs.py`, numbers to `summary.json`):
- **S (CSV convention):** for each row with a deflection, the source equals the image minus the deflection. Read both
  sign conventions (x along +RA or −RA) and report the one with the smaller median mismatch; it passes if that median
  is ≤ 0.01″. ASSUMPTION: deflections are in arcsec at the image's own z.
- **D (deflection from κ):** α from the κ map by FFT convolution on a zero-padded (2×) grid, scaled from z_s = 3 to each
  image's z with D_ds / D_s (flat ΛCDM, Ωm 0.3, h 0.7, the README's cosmology), bilinear at the image position,
  compared with the CSV `deflect_*_best`. Passes if median |Δα| ≤ 0.10″ and 90th percentile ≤ 0.25″ (ASSUMPTION:
  small next to the 0.61″ model rms).
- **D-affine (only if D fails):** the κ map covers 200″, so mass outside it is missing. Add an affine term
  α += a + B θ (6 parameters) fitted on odd-numbered systems and tested on even-numbered ones; same thresholds on the
  test half. Passing makes the maps usable for gate K with that term; failing stops the field's map route.
- **R (reference catalogue):** the DJA v7 imaging index has a field named for RXCJ2211 (`rxj2211`, `rxcj2211`,
  `j2211`, `rj2211`).
Nothing here reads a gate-K residual or a predicted image.
