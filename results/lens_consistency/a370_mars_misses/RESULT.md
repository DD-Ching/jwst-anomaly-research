# Abell 370: the four MARS held-out misses are reproduced by BUFFALO model A (D-034 addendum aq)

Plan frozen in fae1d495; `classify` and its test in e58a959f; both before the run. Numbers from `summary.json`
(per-image rows in `residuals.ecsv`). Labels: CANUCS positions and listed z are observed; residuals are model
predictions; the 3.0″ radius is an ASSUMPTION (PLAN.md).

## Result (40 Abell 370 CANUCS images held out from MARS, same anchors as addendum ap)
| Class | Images |
|---|---|
| both models within 3″ | 36 |
| MARS only | 0 |
| BUFFALO only | 4 (3b.1 1.45″, 5.3 0.62″, 7.6 0.49″, 24e.2 2.18″) |
| neither | 0 |

- Every MARS miss of addendum ap is reproduced by BUFFALO model A, which gives 24e.2 three predicted images where
  MARS gave one, and 7.6 a residual of 0.49″ where MARS gave 20.91″.
- Reported, not a gate: BUFFALO on all 40 images, rms 1.170″, 100 % within 3″ (gate H's rms threshold 1.0″ is
  not met).

## Reading (narrow)
- The four large Abell 370 residuals of addendum ap are specific to the MARS maps; an independent model with
  0.12″ pixels reproduces them within 3″. They are closed as model-side misses (no flag, no candidate).
- Which MARS property causes them (1.515″ pixels near critical curves, free-form freedom away from the HST-era
  constraints) is not separated here.
- None of the four lies within 1″ of a BUFFALO gold row (PLAN.md); BUFFALO model A was fitted to its own gold list,
  so this is not a fully independent held-out test for BUFFALO's nearby systems (3, 5, 10, 24).
- No pixels read; positions only.

## Next
Nothing further on these four images. MARS stays validated on MACS0416 only (addendum ap).
