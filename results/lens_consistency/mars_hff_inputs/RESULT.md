# MARS HFF free-form maps: C / U pass on six fields, gate K passes on four (D-034 addendum ao)

Plan frozen in 896c4c1 (before any map value), functions and tests in b28d5dd (before the run). Numbers from
`summary.json` (per-system rows there; image-plane tables `imageplane_<tag>.ecsv`). Labels: catalogue positions and
spectroscopic z are observed (as listed by the authors); deflections are model predictions; thresholds are
ASSUMPTIONs (PLAN.md).

## Search (`search.json`)
A Zenodo search on 2026-10-24 found no JWST-era cluster deflection-map record that this route has not used or
listed. The one unused cluster-map record with an image list is HST-era: Cha & Jee 2023 MARS (10.5281/zenodo.7804575,
2.7 MB zip, md5 checked, deleted after the run). Record 13843384 ("El Gordo … HST and JWST") has no files.

## Result (per field; spec images = gold rows without a MARS-fitted z, ≥ 2 per system)
| Field | Spec images / systems | Map pixel | U median A / P | Gate K rms (mean source) | Within 3″ | K | K-WLS (reported) |
|---|---|---|---|---|---|---|---|
| Abell 370 | 105 / 35 | 1.515″ | 0.053″ / 78.3″ | 0.850″ | 98.1 % | pass | 0.942″, 97.1 % |
| Abell 2744 | 120 / 26 | 1.010″ | 0.041″ / 0.760″ | 1.112″ | 96.7 % | **fail** | 1.051″, 96.7 % |
| Abell S1063 | 57 / 20 | 1.515″ | 0.102″ / 62.0″ | 8.020″ | 96.5 % | **fail** | 0.078″, 100 % |
| MACS0416 | 236 / 67 | 1.263″ | 0.039″ / 17.0″ | 0.513″ | 100 % | pass | 4.520″, 94.5 % |
| MACS0717 | 31 / 10 | 2.020″ | 0.061″ / 138.2″ | 0.099″ | 100 % | pass | 0.075″, 100 % |
| MACS1149 | 123 / 40 | 0.758″ | 0.084″ / 10.1″ | 0.742″ | 99.2 % | pass | 0.336″, 100 % |

- C passes everywhere (every spec image inside the map). U accepts only the arcsec reading (as the README says) on
  all six fields.
- **Record outcome: usable on 4 / 6 fields** (Abell 370, MACS0416, MACS0717, MACS1149).
- Abell S1063 fails on system 4 alone: under the mean source it gets one predicted image (rms 34.95″); with the WLS
  source the field gives 0.078″. Abell 2744 fails on rms (1.112″ > 1.0″); its largest per-system rms are systems 5
  (3.73″) and 40 (2.90″).
- The authors' median `Delta_rms` of the same images is 0.004–0.035″; our values are larger. Ordinary causes not
  separated here (not tested): bilinear interpolation of 0.76–2.02″ maps (the authors interpolate bicubically) and the
  mean-source estimator.
- MACS0416's K-WLS (4.52″) is far worse than its mean-source value (0.51″); not diagnosed.

## Not done / limits
- Code fix after the run (review finding): the ≥ 2-images rule is re-applied after C drops images outside the map. No
  effect on these numbers, since C kept every image on all six fields.
- No extra-image, flux or flag reading (frozen). These are HST-era fits to these same images, so the stage tests
  our loading and the maps, not lensing consistency.
- The catalogue labels MACS1149 "MACS J1149.5-2223" (the map header CRVAL2 is +22.40°); the script maps it to the `MACS J1149.5+2223`
  folder by name.

## Next
A held-out stage on the four usable fields: JWST-era spectroscopic images that are not in the 2023 MARS catalogue
(e.g. CANUCS DR1 MACS0416 / MACS1149 / Abell 370 lists already pinned here), with a frozen plan: does a model fitted
before JWST reproduce them?
