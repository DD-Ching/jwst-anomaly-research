# Acceleration V2 audit, first pass (2026-10-10)

The owner brief "Autonomous research acceleration V2" asks for a measured audit before any optimization, without
repeating what the repository already covers. This file records what was measured, which parts of the brief are
already in place (with their canonical locations), and the next units, ranked by expected information per hour.
Everything below is measured on this date unless labelled ASSUMPTION or estimate.

## Already covered (not repeated)
| Brief item | Where it lives |
|---|---|
| Profile before optimizing; background fetches; durable derived-data store | D-078 (S-1 to S-5) |
| Claim, heartbeat and takeover rules; one unit per session | D-078 S-5, docs/cloud-routine-prompt.md, research-cycle skill |
| Inject through the whole chain; limits only where injections recover; thresholds as ASSUMPTIONs | scripts/CLAUDE.md, w3-survey skill §5 |
| Pipeline overlap of I/O and compute; process pools sized to the cores | scripts/CLAUDE.md "Parallel topology" |
| Repository as the shared knowledge index (decisions, failed approaches, sources) | DECISIONS.md, CHANGELOG.md, SOURCES.md, skills' failed-approach lists |
| Theory program with its own budget; pause rule | D-069, D-075, D-076 (rounds paused until a new owner brief) |
| Excluded vs undetectable vs untested | docs/exotic_limits.md per search ("valid only over systems where injection proved...") |

## MOA-II W3 pipeline: measured baseline (gb12, CHAIN_VERSION 2026-10-10.1, 4 cores)
Run time 1,528 s. Pre-screen stream 342 s (73,386 light curves, ~215 per s, CPU 99 %, ~80 MB/s). 51 fits 203 s.
Vetting of 46 flags ~2 s. 2,200 injections ~950 s (≈ 60 % of the run).

**Where injected W3 signals are lost** (2,000 W3 injections with `--sampling lf`, `injections_gb12.ecsv`):

| stage | kept | of previous |
|---|---|---|
| Cut-0 (published selection, emulated) | 904 | 45 % |
| shape cut (`prescreen_z/s/repeat`) | 212 | 23 % |
| shared-epoch test | 158 | 75 % |
| fit flag (ΔBIC) | 156 | 99 % |
| vetting chain | 99 | 63 % |

By source magnitude I_s (kept at Cut-0 → shape → recovered): 14–17: 68 → 53 → 32; 17–18.5: 139 → 78 → 39;
18.5–19.5: 206 → 45 → 20; 19.5–20.5: 322 → 25 → 6; 20.5–21.5: 169 → 11 → 2. Most of the loss happens at Cut-0
and at the shape cut for faint sources (I_s > 19.5 holds 76 % of the LF-weighted injections and 8 % of the
recoveries), where the signal is near the noise. In those cases the loss is mostly physical, not a filter defect.
Two losses are design choices and can be tested:
1. **Shared-epoch test: 25 % of shape-passing injections.** It removes 99 % of the real shape passes
   (5,532 → 51), so it is the main background filter. Injected epochs are random, so the loss is coincidence
   between a real dip box and frames that many stars share. Candidate improvement (needs a dev/validation split):
   mask the shared epochs and re-scan, instead of vetoing the light curve.
2. **Vetting at short and long t_E.** By t_E (flagged → recovered): 3 d 23 → 9 (`jackknife_nights` 7: one or two
   nights carry the spikes when sampling is sparse); 300 d 22 → 9 (`feature_bracketed` 10: a 300-d umbra plus the
   20-epoch baseline on both sides rarely fits inside a season gap pattern). 10–100 d: 72–76 % kept. This answers
   the brief's question "can a real peak-dip-peak lose its peaks to sparse sampling and be rejected as an ordinary
   dip?": yes, measurably at t_E = 3 d through `jackknife_nights`. At other t_E the losses go to
   `repeated_deficit` (9) and `exotic_feature_sampled` (5).
   These tests exist because real false positives passed without them (w3-survey failed-approach rules). Any change
   needs the background flags (real, 46 per field) as the false-positive check. It also needs injections it was
   not tuned on.

## Redundant computation found and removed
- `_chunk_done` compared the whole `Params`, so the D-068 addendum's vetting-only parameters invalidated every
  streamed pre-screen chunk. gb12 re-streamed 342 s and gb7 487 s, and the rows came out byte-identical to the
  tracked chunks (checked for gb12 chunks 1, 4, 7; only the metadata changed). Fixed in this PR: chunk validity
  now depends only on the pre-screen parameters, the tracked-row constants and a `PRESCREEN_CODE` version.
  - Before: a field re-run re-streams its whole tar, about 340–490 s and 24–30 GB over HTTP.
  - After: 0 s for unchanged chunks.
  - Measured for the three remaining fields (gb15, gb17, gb18): all 31 tracked chunks pass the new check.
  - Scientific impact: none (identical rows). Validation: unit test for vetting-only changes, pre-screen changes,
    code bumps and legacy chunks; full offline suite.

## Next units, ranked (qualitative information per hour)
1. Remaining MOA fields gb15, gb17, gb18 with the reuse fix (~20 min each instead of ~26), then `combine`.
   This tightens an existing limit at known cost.
2. Shared-epoch masking study. Split the injections into dev and validation sets. Compare (a) veto vs (b) mask
   and re-scan, on both injections and the real shape passes. Ship only if recovery rises at a non-increasing
   real-flag count. Up to +25 % efficiency at the pre-screen; the vetting cost of the extra flags is unknown.
3. Injection speed. Injections are about 60 % of a field run. Profile `run_inject` before choosing (the fits
   dominate per injection; early rejection of injections that fail Cut-0 or the shape cut is already in place).
4. Short-t_E jackknife study (t_E = 3 d keeps 9 / 23): only with the background flags as the false-positive check.

Not done here, and why: GPU (the fits are small scalar optimisations with branchy vetting; no profile suggests a
GPU-shaped kernel); adaptive grids for the W3 limit (the 5 × 2 grid is set by the published-limit format, and
cells are cheap compared with the streaming).
