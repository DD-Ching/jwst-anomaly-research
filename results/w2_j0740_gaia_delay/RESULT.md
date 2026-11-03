# J0740+2926: can Gaia DR3 G epoch photometry identify a delayed replay of one source? (2026-11-03)

Owner-requested single-pair test, at most 60 minutes.

**Sources:**
- Pre-registration: `PREREG.md`, frozen at commit 4ef8f8c (06:11:57 UTC) before any light curve of the pair was
  fetched. Its text says "~06:15"; the commit time is the record.
- Numbers:
  - `result.json`: the frozen run;
  - `post_hoc.json`: labelled post hoc;
  - `sis.json`: a model prediction;
  - `literature.json`.
- Script `scripts/w2_j0740_gaia_delay.py`; contact sheet `contact_sheet.png` (looked at).

**Labels:**
- Gaia light curves: observed.
- Delays, D values and exclusions: derived.
- Controls: simulated.
- SIS delays: model prediction.

## Verdict

**Frozen rules:** "excludes the listed delay intervals". Neither named range is excluded: R_lens (B lags by +20 … +400 d)
and R_all (−400 … +400 d) both stay open.

**No support for a same-source delay:**
- C1 failed at every injected delay. Known delays were recovered within ±20 d in only 20–37.5 % of simulations; the
  pass mark is 68 %.
- C2 failed. All 100 simulated independent pairs passed the same-source acceptance (simulated false-pass
  fraction 1.0; not a sky rate).
- The two PyCS3 estimators disagree on the target: the spline gives +255 d, regdiff gives −88 d.

**What the data can and cannot do:**
- Gaia DR3's two G series cannot identify a delayed replay for this pair.
- They exclude scattered delays only under the simulation assumptions (§4). Across the error-scale uncertainty, only
  B leading A by 160–220 d stays excluded.
- For the delays an ordinary galaxy lens would give, the data are insufficient.

## 1. Published delay (literature, `literature.json`)
- **Measured delay: not found.**
  - J0740+2926 is not among the 73 systems in Gaia GraL X's delay compilation (arXiv:2603.12810v1).
  - No monitoring paper is among its SIMBAD references.
- **Model-predicted delay: not found.**
- **Published classes:**
  - NIQ (Lemon+2018: same redshift, no lens galaxy detected).
  - Binary quasar (Hennawi+2006). Kayo & Oguri 2012 give |Δv| = 230 km/s.
  - **"cont."** in Shajib+2021's Keck follow-up table (arXiv:2011.01971, App. A). That class means "revealed as
    non-lenses … either through our high-resolution imaging or spectroscopy". There is no object-specific note, so
    the basis is not stated.
- **This round's own SIS illustration** (`sis.json`, model prediction):
  - image radii 2.04″ / 0.60″ from the Gaia G ratio 3.42 [assumption: the G ratio is the macro-magnification ratio];
  - with the lens at z 0.2–0.7, B lags by 60–603 d.

## 2. Data (observed; Q1–Q3 as frozen)
| | A (878723302007625344) | B (878723302009111040) |
|---|---|---|
| G epochs kept / DR3-rejected | 22 / 2 | 25 / 0 |
| visit groups (≤ 1 d) | 16 | 17 |
| median error after Q2 | 0.017 mag | 0.040 mag |

- **Error scale (Q2):** s = 1.73, from χ²_ν = 3.01 about the means of 14 same-day groups (14 dof). Its 95 % range
  is 1.27–2.74 (post hoc, `post_hoc.json`).
- **Span and shared transits:** the span is 971.6 d, and 22 of B's transits have a partner transit of A.
- **Control pair (C3):** the two nearest Gaia DR3 `vari_agn` sources within 1° with G 18.0–19.8, run as a pair:
  - 878708630399208960, G 19.43, 0.127°;
  - 878728730846385664, G 18.60, 0.149°;
  - its own error scale is 1.32.

## 3. Controls (simulated on the actual epochs and errors; `result.json`)
**C1, recovery** (PyCS3 spline, 40 simulations per delay). Fraction recovered within ±20 d:

| Injected Δt (d) | −200 | −100 | +50 | +100 | +200 | +300 |
|---|---|---|---|---|---|---|
| Fraction within ±20 d | 0.275 | 0.375 | 0.225 | 0.225 | 0.200 | 0.200 |

None reaches the 68 % pass mark, so C1 fails.

**Thresholds** from D(Δt_true) over the 240 C1 simulations: T95 = 35.5, T99 = 48.0; the median is 10.4 (post hoc).

**C2, independent pairs** (N = 100):
- Every pair gave D(Δt̂) ≤ T95; the largest was 18.9 (post hoc). The simulated false-pass fraction is 1.0, so C2
  fails.
- The fitted delays pile up at the range edges: 66 of 100 have |Δt̂| ≥ 200 d (post hoc), where the curves
  overlap least.
- The exclusion scan rejects a median 36 % of the range for an independent pair (0–84 %; 40 pairs). This is the
  power of the exclusion test.

**C3, real independent pair:** D(Δt̂) = 35.54 against T95 = 35.47. It passes, but only narrowly.

**Reading:**
- With 16–17 visits per image, the spline finds a shift that fits two unrelated curves as well as, or better than,
  a true same-source pair at its true delay.
- D therefore cannot detect a replay here.

## 4. Target
**Delay estimates:**
- Spline Δt̂ = +255 d, with D = 3.5 (r² 21.8 joint against 18.3 separate).
- Regdiff Δt̂ = −88 d.
- Neither is evidence, because C1 and C2 failed.

**Exclusion:** D(Δt) > T99 on the 5-day grid, with the frozen zero-lag rule applied post hoc. Excluded intervals
(d): −400…−350, −330, −250…−135, −105…−65, −20…−5, +5…+10 and +320…+380.
- That is 40 % of the interpretable grid and 17 % of R_lens.

**Post hoc sensitivity to the error scale:**
- Rescaling all errors multiplies D by (1.73/s′)², with T99 kept.
- Post hoc: at the 95 % upper scale s′ = 2.74, only −240…−235, −220…−160 and +375 remain (10 % of the grid, 1 % of R_lens).
- Post hoc: at the lower scale s′ = 1.27, 58 % of the grid and 31 % of R_lens are excluded.

**Assumptions the exclusion rests on:**
- intrinsic variability is a damped random walk (τ = 300 d, σ = 0.10 mag);
- no microlensing beyond a linear trend in B over 972 d;
- common per-transit errors of 0.01 mag;
- 1 % A→B crosstalk.

**Feature to check against microlensing (post hoc):** A's last-season median is 0.21 mag fainter than its earlier
median; B's changes by −0.0668 mag (brighter). If A's change were microlensing rather than intrinsic, the exclusions that require it to recur
in B would not hold. This was not tested.

## 5. Deviations and notes
- **Implementation deviation (after the first run aborted, before any result):**
  - PyCS3's `opt_fine` default bokeps = 10 raised "Knot spacing min < epsilon" on these sparse curves.
  - The rerun uses bokeps = 2.0 (PyCS3's `fit` default). It skips failed starts and records failed fixed-delay fits
    as not evaluated; 0 were.
- **Iteration counts:** the PyCS3 iteration counts (3 per start, 10 for refinement, 5 for fixed delays) and C3's role
  order (the nearer source as A) were set in code before the run. PREREG does not name them.
- **Zero-lag rule:** the run's exclusion omitted PREREG's |Δt| < 5 d rule. `post_hoc.json` applies it. The only change
  is that −20…+10 becomes −20…−5 and +5…+10.
- **Post hoc additions:** the error-scale sensitivity, the C2 alias counts and the last-season changes.
- **Not done:** no second variability model, and no microlensing simulation beyond a linear trend.

## 6. What this does not say
- It neither confirms nor rules out lensing. Shajib+2021's class is the strongest published statement, and its basis is
  not given.
- A delay test for the lens range would need denser resolved photometry of both images, 2.6″ apart. None was
  examined here.
- J0740+2926 is not promoted: it stays a published NIQ / binary / contaminant system, not a lens candidate from this
  test.
