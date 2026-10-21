# Parallel round 2026-10-21: X4/R1, S1 CHIME, W2 HST (coordinator plan, frozen before any line runs)

Owner request (interactive, 2026-10-21): a manual, bounded reopening of three lines, run in parallel by up to three
worker subagents under one coordinator. Routines and schedules are not changed. This file is the round's
arrangement only; it does not change the permanent collaboration policy in CLAUDE.md.

## Budgets (owner)
- Wall clock 6 h for the round: 00:05–06:05 UTC. Downloads 8 GB in total. No paid resources, accounts or new
  schedules.
- A (X4/R1): compute ≤ 5 h. B (S1 CHIME/FRB Catalog 2): ≤ 3 h. C (W2 HST image test): ≤ 3 h. Any line may stop early.

## Machine (measured 2026-10-21 00:06 UTC)
- 4 CPU cores, 15 GiB RAM, no swap, about 30 GB writable disk.
- Allocation, so heavy jobs do not slow each other:
  - A: at most 3 concurrent CPU-bound processes (one fit per process, `OMP_NUM_THREADS=1`), ≤ 8 GB RAM; a 4th
    process only after B and C have finished.
  - B and C: mostly network I/O; each runs at most one CPU-bound process at a time, under `nice -n 10`, ≤ 3 GB RAM.
- Downloads: B ≤ 5 GB, C ≤ 2.5 GB, A none (the lalsuite wheel only). Data are streamed or deleted after use (cloud
  disk rule); each line logs URL, bytes and sha256 of every download.

## Hard stops (UTC)
- B and C: 03:30 (results saved, write-up done).
- A: computation stops at 05:15; write-up done by 05:35.
- Coordinator: state files, review and PR by 06:05.

## Lines, working and result directories
| Line | Unit | Scripts | Results | Scratch |
|---|---|---|---|---|
| A | X4 S1g: bounded P-search budget, validation, then R1 | `scripts/x4_s1g.py` | `results/x4/s1g_*` | `scratchpad/x4/` |
| B | S1 on CHIME/FRB Catalog 2 dynamic spectra, ≤ 30 events | `scripts/s1_chime*.py` | `results/s1_chime/` | `scratchpad/s1-chime/` |
| C | W2 PSF-subtracted HST test of blended or close lensed quasars, ≤ 10 targets + controls | `scripts/w2_hst_blend*.py` | `results/w2_hst_blend/` | `scratchpad/w2-hst/` |

- Each line writes and commits its own frozen plan (`PLAN.md` in its result directory, or the unit-file section for
  A) before it looks at any new result. Implementation bugs may be fixed and re-run; both versions are kept.
- The coordinator alone edits TASKS.md, CHANGELOG.md, DECISIONS.md, SOURCES.md, `docs/state/trial.md` and
  `docs/state/pending.md`; workers hand their text for those files back in their reports.

## Checked before starting (orientation, 00:05–00:28 UTC)
- `main` b26b62b. Open PRs: #418 (draft, MACS1149) and the needs-human PRs in `docs/state/pending.md`; none covers
  X4, CHIME Catalog 2 or a pixel-level HST test of the 307 blended or close quasar lenses. No `claude/*` branch does.
- X4: S1f was the last stage (R1 undetermined; S1e's P search failed V9 seed 11). Not redone: S1–S1f.
- S1: #167 (Fermi M2 re-screen) and #169 (TTE, E_peak) are done on unmerged PRs and are not redone.
- W2: D-056, D-060, D-064 and the D-064 HST addendum (`scripts/w12_niq_hst.py`, `w12_niq_epsf.py`) are reused, not
  redone.
- Kept paused: ASAS-SN, S6a expert review, E1 CHIME uptime, Gaia DR4, S2 (no real GW events in this round), every
  parked candidate. No fourth line is opened.

## Environment
- Project env `.venv` (Python 3.12, `-e ".[dev,cloud]"`).
- A's separate waveform env: lalsuite 7.26.16 from PyPI, which provides lal 7.7.1 and lalsimulation 6.2.1, the
  versions X4 used. The earlier record "lalsuite 7.7.1" names the lal version.

## Actual use (filled in by the coordinator at the end)
| Line | Wall clock (UTC) | Compute | Downloads | Outcome |
|---|---|---|---|---|
| A | 00:24–03:06 (fits 00:49–03:03) | 6.8 CPU-h on 3 cores | none (PyPI wheel only) | P search passes 6 / 6 frozen checks; R1 seeds disagree, R1 undetermined |
| B | 00:22–00:52, review fixes to about 02:00 | about 5 CPU-min | 0.90 GB, streamed and deleted | frozen S1 method does not transfer to CHIME (recovery 6.9 %) |
| C | 00:22–01:16, review fixes to 01:52 | about 2 CPU-min, plus downloads | 2.28 GB, deleted | negative control S/N 13 ≥ 5σ: method validation fails, targets not interpreted |

- Round total: about 3.2 GB downloaded of the 8 GB budget. No 4th concurrent process was used.
- Two /code-review rounds on B and C changed no committed verdict (B: `chi2_union` δ choice; C: depth counted once per
  visit, which changed the selection but not the gate or the negative control). Both v1 result sets are kept.
- Owner clarifications during the round (2026-10-21): freed cores may only run planned fits earlier; B's conclusion is
  limited to "the frozen method did not transfer to CHIME"; a negative control at threshold fails C's method
  validation; R1 only after the whole validation set passes and within the frozen budget.
