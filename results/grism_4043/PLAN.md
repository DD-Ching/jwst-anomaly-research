# Program 4043 F444W grism: coverage of SMACS candidate records and a line test for f200w_2915

Written before any grism flux was read (footprints and c1d source positions only).

## Step 1 (done before this plan): coverage
`scripts/grism_4043_coverage.py` places the nine `results/smacs_jwst3/candidates.csv` records and the two point-like
top-20 sources (`f200w_1816`, `f200w_384`) against the level-3 grism footprints (MAST ObsCore `s_region`) and the
pipeline c1d source list (jwst 3.0.0, `SCATFILE` = the 4043 F444W catalog). Only `f200w_2915` has a c1d spectrum
(0.03″ in both dispersions); `f200w_2925` matches at 0.31″, beyond the 0.15″ radius (ASSUMPTION, as in
`results/smacs_jwst3/RESULT.md`), and has no 2736 counterpart anyway.

## Step 2: line test for f200w_2915
- **Hypothesis under test:** the recorded verdict (`docs/candidates/jw02736-o001_t001_nircam_f200w_2915.md`): a
  quiescent galaxy at DJA z_phot 2.84 (95 %: 2.68–3.10).
- **Ordinary alternative:** a dusty star-forming galaxy or AGN at the same z (DJA A_V 1.8), which would show
  Paschen or He I emission. Lines in band for z 2.68–3.10: He I 1.0833 µm, Paγ 1.0941 µm, Paβ 1.2822 µm
  (rest, vacuum).
- **Statistic (frozen, `line_snr` in `scripts/grism_4043_lines.py`):** continuum = running median, 51 pixels;
  residual summed in a 5-pixel boxcar; S/N = boxcar sum / robust σ of boxcar sums over the spectrum
  (1.4826 × MAD, finite pixels with DQ = 0 and N_INPUT > 0), so the S/N is empirically normalised.
- **Detection rule (frozen, `decide`):** a line is detected when, in **both** grismr and grismc, the maximum S/N
  inside the window of that line for z 2.68–3.10 is ≥ 5 (ASSUMPTION) and the two peak wavelengths agree within
  3 nm (ASSUMPTION). Any detected line → verdict `emission_detected` (quiescent reading contradicted; re-vet);
  none → `no_line_detected` (consistent with the verdict, stated with the injection limit).
- **Validation gate (before reading the result):** inject Gaussian lines (σ 1.5 nm) of boxcar S/N 7 at 200 random
  z in 2.68–3.10 for each line into both spectra; the rule must recover ≥ 90 % (ASSUMPTION). Zero-amplitude
  injections give the false-positive rate over the same windows. If the gate fails, no verdict.
- **What changes the plan:** `emission_detected` re-opens the 2915 vetting record; `no_line_detected` closes the
  grism route for the SMACS candidate records (none of the others is covered).
