# W3 ZTF: ZTF alone vs ZTF + ATLAS, paired injections on the 40 calibration stars (feasibility trial)

Owner brief 2026-10-22 after #461:
- 4 h wall clock, ≤ 2 GB, ≤ 4 cores. Session start 08:56 UTC.
- ATLAS through the documented shared Astro Data Lab token. At most 2 of our tasks queued or running; only our own tasks
  are touched. No account, no e-mail.
- A feasibility trial on the already-used calibration sample, not an independent validation. No search group, no
  search.
- X4, ASAS-SN and the other paused units are untouched.

## 0. Prior state (checked 08:56–09:00)
- `main` at b85ca54 (#461 merged). Open PRs: #442 (dependabot) and the `needs-human` backlog; none on W3 / ZTF or ATLAS.
- #461 kept only 498 in-gap epochs of 4 stars (`../atlas/analysis.json`). The full light curves were deleted, so they
  cannot be recomputed and all 40 stars are requested again. The 4 #461 stars are **development data**:
  601209200002983, 601213300003718, 632215200005622 and 632211200001819.
- The official pages (output description, FAQ, API guide, schema), re-read 08:57, still say nothing about the value
  9999. The NOIRLab notebook still allows the Data Lab token.

## 1. Acquisition (committed before any request; `targets_script.py.txt` → `targets.json`, `acquire_script.py.txt`)
- **Targets:** the 40 stars of `../asassn/selection.json`, one request each, reduced mode (total flux), MJD 58150–61000.
  That covers both areas' ZTF carrier ranges (58198.4–60969.3).
- **Priority order** (committed counts only): most A0-insufficient W3 events, then most W3 events, then most fades,
  then key. No star is swapped for ATLAS coverage or results.
- **Limits:**
  - ≤ 2 tasks in flight, polled every ≥ 30 s; 429 honours `Retry-After`.
  - Never resubmitted.
  - A task unfinished 45 min after submission is recorded as `timeout` and left alone.
- **Acquisition deadline:** 170 min from the start of the scheduler. Unfinished tasks are then recorded for
  continuation, not as "no data".
- **Kept for recomputation:** compact per-star measurement tables (< 1 MB each; the columns needed by the quality cuts,
  time and flux), task ids, time stamps and source sha256. The token is in no product.
- The method (quality, time, baselines, errors, injection, rules, denominators, success and stop) is fixed in §2–§5,
  committed before any ATLAS flux of this round is read.
