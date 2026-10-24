# RXCJ2211-0350 input check (2026-10-24)

Next field for the fitted-model route (TASKS: "another field with public JWST-era maps"), found by a Zenodo search
this run: the VENUS GLAFIC v0.91 model, Zenodo 10.5281/zenodo.20312327 (CC0, published 2026-05-20; the record cites
arXiv:2512.05180; the README's model reference is "in prep."). Checks frozen in `PLAN.md` (e9ccd27), functions and
tests in a separate commit before the run. Script `scripts/rxcj2211_inputs.py`, numbers from `summary.json`,
image list `rj2211_images.ecsv`, tests `tests/test_rxcj2211_inputs.py`.

| Input / check | Found | Status |
|---|---|---|
| Image list | `full_output_v91.csv`: 146 rows, 43 systems; 10 rows (systems 5, 7, 43) with z −1; README `ZS_REF` (spec-z reference) only for systems 1, 14, 32 | present; most z are the model's, not spectroscopic |
| Model fit | README: χ² 174.99 / 180 dof, rms 0.61″, z_l 0.397 | reported by the authors |
| Model maps | κ, γ1, γ2 only, 2000 × 2000 at 0.1″ (200″ field), RA---CAR, header ZS = 3; median κ on the map edge 0.208 | **no deflection maps released** |
| S: CSV source = image − deflection | median mismatch 0.00016″ with x along −RA (70.8″ with +RA) | **passes** (≤ 0.01″) |
| D: α from κ (FFT, zero-padded, scaled from z_s 3 by D_ds / D_s) vs the CSV deflections, 136 images | median 14.13″, 90th percentile 15.67″, max 20.32″ | **fails** (≤ 0.10″ / ≤ 0.25″) |
| D-affine: + 6-parameter affine term fitted on odd systems (70 rows), tested on even (66) | test median 2.66″, 90th percentile 4.35″ | **fails** |
| R: DJA v7 imaging field (35 listed) named rxj2211 / rxcj2211 / j2211 / rj2211 | none | **fails** |

**Reading.** As frozen, the κ map, read at its header z_s = 3, cannot stand in for the missing deflection maps, so
the map route for RXCJ2211 stops here, and gate C2 could not run anyway without a public deblended reference. No
gate-K residual and no predicted image were read. The plan declared only the header normalisation; unlike the MACS
J1931 unit check U it did not state both readings before the run, so this failure may come from an input convention
rather than from the model (see the post hoc line).

**Review fixes (/code-review, before the PR):** pixel positions from the header WCS (`astropy.wcs`) instead of a hand
offset, images outside the map refused (136 / 136 inside), ids read as strings, `sample`, `dls_ds`, `read_csv` and
the README parser unit-tested, the cosmology stored, and the post hoc block written only when D-affine fails. The
verdicts and numbers are unchanged to 0.001″.

**Post hoc (added after D and D-affine failed; no verdict depends on it):** one scale factor fitted to all 136 rows
is 0.763, close to D_ds / D_s at z_s 3 (0.7575); after it the median mismatch is 1.71″ (90th percentile 2.42″). So the
κ map may be normalised to D_ds / D_s = 1 rather than to its header's z_s = 3 (hypothesis, not tested), and a
remaining ~1.7″ would still need mass outside the 200″ map. Any retry is a new plan frozen before its run (scale term
plus affine, cross-validated as here), or the authors' deflection maps if they are released.

**Data:** the 43 MB zip was read in memory and in the session scratchpad only (md5 matches the record), then deleted;
nothing stored beyond the image list and the summary.
