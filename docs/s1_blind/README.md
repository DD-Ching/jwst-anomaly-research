# S1 blind spots: single-pulse, faint and short bursts

Owner brief (2026-10-11): bring the 2,975 bursts that S1 excludes into the search with TTE. Plan, with decision
rules written before any score: `docs/s1_blind/plan.md`. Code: `scripts/s1_blind_spots.py`; results
`results/s1_blind/`; test `tests/test_s1_blind_spots.py`. No candidate search was run on these bursts. This round
asks whether one can be informative at all.

## 1. Why the 2,975 are excluded [derived]

| Reason (bcat light curves) | Bursts | TTE can help? |
|---|---|---|
| No resolved pulse (no peak at 5σ prominence): faint | 846 | Barely. TTE has the same photons; finer bins lower the S/N per bin. |
| One resolved pulse, T90 > 10 s | 1,197 | No. The 64 ms–1 s bins already resolve pulses lasting ≥ 10 s. |
| One resolved pulse, T90 2–10 s | 410 | Little. |
| One resolved pulse, T90 < 2 s (short) | 522 (699 short in total, against 42 eligible) | Possibly, since 64 ms is coarse for short bursts. |

## 2. Step A: can single-pulse bursts be searched at all? [obs + sim; `single_pulse.json`]

- **Setup:** the frozen M2 method (#167), applied to npulse = 1 bursts with peak S/N ≥ 10 (1,167 bursts).
  - Threshold on the calibration half's real-pair background (130,202 pairs, per-pair FAR 10⁻⁴): **0.9857**,
    against 0.9716 for multi-pulse bursts.
  - Injections: the frozen deformation set on the other half, with the copy required to have ≥ 1 pulse.
- **Unrelated pairs pass the consistency test at 39 %** (top 100: 47 %).
  - Check (owner review): among the calibration pairs that are actually above the threshold, 5 / 12 pass (42 %,
    Wilson 95 % about 19–68 %).
  - A full single-pulse screen (517,629 pairs, about 52 above threshold) would therefore give **about 20 false
    candidates (rough range 10–35)**. This is an extrapolation, not a measured count.
- **Deformed injections:** 46 / 960 recovered, **4.8 %** total. Undeformed: 4 / 120.
- **Decision (pre-declared):** both criteria fail (false candidates ≤ 3; deformed recovery ≥ 10 %).
  - Once stretch, lag and per-band scale are free, one smooth pulse is matched by many unrelated single pulses.
  - **Scope (owner review):** this holds for the current M2 statistic, two energy bands, this sample (npulse = 1,
    S/N ≥ 10) and this injection set. It is **not** a general statement that single-pulse bursts cannot be
    searched with any data product. A statistic that also constrains the spectrum might separate them (§5).

## 3. Step B: TTE pilot for short bursts [obs; `tte_pilot.json`, `tte_control.json`]

- **Data:** NaI TTE of the bcat detectors (HEASARC triggers directory), streamed one burst at a time into scratch
  and deleted after reduction.
- **Light curves:** two-band (50–300 keV; 8–50 + 300–900 keV) background-subtracted count rates, 2–16 ms bins.
  The short-burst window is padded by 0.5 T90 + 0.1 s (a pilot setting).
- **Pulse rule:** the same 5σ prominence rule as S1.

| Sample | n | ≥ 2 pulses, finest bins | ≥ 2 pulses, multi-resolution (2–64 ms) |
|---|---|---|---|
| Control: short bursts already resolved by bcat (≥ 2 pulses) | 15 | 9 (60 %) | **13 (87 %)** |
| Pilot: excluded short bursts with bcat peak S/N > 20 (seed 5, from 147) | 30 | 3 (10 %) | **5 (17 %)** |

- **Failed approach, kept on record:**
  - The first window padding (2 s, the S1 default) forced 32 ms bins for 0.1 s bursts.
  - The finest-bin count then lost 40 % of known multi-pulse structure in the control, because per-bin S/N drops.
  - Counting at several bin widths and keeping the maximum restores 87 %. This fix was made after seeing the
    control, and is labelled as such.
- **Decision (pre-declared ≥ 25 % gain):** 17 % falls short, so the full TTE run was paused **for resources**.
  - 5 / 30 has a Wilson 95 % interval of 7–34 %. The data do **not** show that the true gain is below 25 %.
  - The multi-resolution rule was changed after seeing the control.
  - After the round-2 robustness check (§5), 3 of the 5 survive: 3 / 30, 95 % interval 4–26 %.
  - *Withdrawn:* "these could only be paired with each other". The compatible-pair count is in §5.
  - **Cost** (measured): 497 MB and 88 s for 30 bursts, so about 0.6 h and about 10 GB streamed for all 699.
    It is cheap, but the yield is low.

## 4. Answer

- **Single-pulse bursts (2,129):** not searchable for deformed twins with **this** test (M2 shape match, two
  bands). Unrelated pairs pass too often (39–42 %). With a spectral constraint the answer could differ.
- **Faint bursts with no resolved pulse (846):** TTE does not add photons, so they stay outside a shape-based
  search.
- **Short bursts (699):** TTE resolves robust extra pulses in about 10 % (4–26 %) of the bright ones (§5). A full
  run is affordable; §5 weighs what it would add.
- These are blind spots of shape matching itself. Reaching them needs **another observable**, not more of the same
  light-curve data. Examples:
  - spectral-lag and E_peak(t) relations (CSPEC/TTE spectral fits) as the matching quantity;
  - a different messenger.
  This is a direction, not a result.

## 5. Round 2 (2026-10-12): does TTE time and energy information improve discrimination?

**Files:**
- `docs/s1_blind/plan_round2.md`: rules written before the checks.
- `scripts/s1_tte.py`, `tests/test_s1_tte.py`.
- Results: `results/s1_blind/tte/` (`manifest.json`, `robustness.json`, `epeak_test.json`, `candidates.json`, and
  reduced 8 ms light curves `lc_*.npz`).
- Downloads this round: **965 MB**, under the 1 GB cap. The NaI TTE of 5 new events, 8 pilot single-pulse bursts
  and 34 of the 42 eligible short bursts (a random subset that fits the cap) is cached in scratch and not
  committed. File names, versions and sha256 prefixes are in the manifest.

**5.1 Are the five new multi-pulse events real?** [obs]
- **Method:** 12 one-factor variants per event (bin offsets, background windows and model, detector subset).
  Exposure is checked against the TTE GTIs.
- **Failed approach, fixed:** the alternative background window started at −40 s, but older TTE starts at −26 s or
  −31 s. The no-data interval was counted as zero background. The background fit now uses only GTI-covered 1-s
  bins, and the exposure rule requires ≥ 10 s of covered background per side. This changed after seeing the first
  result.
- **Negative control:** 130 simulated smooth single pulses with independent Poisson noise, matched in peak and
  background rate to the 5 events and 8 pilot bursts, run through the full multi-resolution, all-offset counter.
  **0 / 130 are counted as ≥ 2 pulses** (95 % upper limit 2.9 %). The 5σ threshold stands.
- **Positive control:** 30 / 34 eligible short bursts recovered (88 %, 73–95 %). This is a property of this
  sample, not a general efficiency.

| Event | T90 (s) | nominal pulses | variants with ≥ 2 | robust |
|---|---|---|---|---|
| bn210528586 | 1.02 | 3 | 11 / 12 | **yes** |
| bn081113230 | 0.58 | 2 | 10 / 12 | **yes** |
| bn250330911 | 0.45 | 3 | 12 / 12 | **yes** |
| bn240905186 | 1.28 | 2 | 6 / 12 | no (depends on bin offset) |
| bn101129726 | 0.58 | 2 | 6 / 12 | no (depends on bin offset) |

**5.2 Joint time–energy constraint** [obs; `epeak_test.json`]
- **Model:** one factor s, with t → s t and E → E/s, so the fluence E_peak of a true copy is E_peak,A / s.
  - E_peak is the GBM catalogue value (fermigbrst best-fitting Band or Comptonized model). The GBM team folded
    it with the detector response and background.
  - s comes from the time-domain M2 refined stretch.
  - A systematic of 0.10 in ln E_peak per burst is assumed.
  - Pass if |z| ≤ 2.
- **Coverage:** 2,803 of 4,392 catalogue bursts have such an E_peak. About half of the screened pairs can be tested
  (1,034 / 2,000 and 25 / 64).

| Sample | pass rate |
|---|---|
| Unrelated, 2,000 most M2-similar calibration pairs (1,034 testable) | **20 %** (18–23 %) |
| Unrelated, the 64 full-catalogue pairs above the M2 threshold (25 testable) | 36 % (20–56 %) |
| Simulated true pairs, true systematic 0.10 (= assumed) | 96 % (95–97 %) |
| Simulated true pairs, true systematic 0.20 | 77 % (75–80 %) |

- **Reading:** on testable pairs the energy constraint cuts unrelated look-alikes by about 3–5× while keeping
  77–96 % of true copies.
  - **The pre-declared criteria are met only if the real systematic is about 0.10.** The 0.10 case is circular,
    since the simulation uses the assumed value. Hard-to-soft evolution is present in the real null pairs. The
    test still removes most of them, because the E_peak ratio must equal the time stretch.
  - **Not shown:**
    - the true E_peak systematic between unrelated detector geometries, which should be checked on bursts with
      several independent E_peak measurements;
    - whether the E_peak test is independent of the χ² shape test;
    - any false-alarm rate below about 1 / 1,000, given the sample size.
- **Time-resolved E_peak(t) from TTE is not done.** It needs response matrices (RSP2) per detector and a spectral
  fitter (the NASA Fermi GBM Data Tools) for a few bright bursts. The catalogue values tested here are fluence
  averages.

**5.3 Limited candidate check for the 3 robust events** [obs; `candidates.json`]
- **Compatible partners** in the whole catalogue (position-inconsistent, T90 ratio within the stretch range
  [0.5, 2], any data product): 397, 292 and 316. Of these, 37, 22 and 15 are bcat-multi-pulse, and 29, 17 and 11
  have TTE downloaded this round.
- **Scoring:** 59 new pairs at one common TTE bin width, each built from its own events. Highest M2 = 0.914,
  against the 0.9716 threshold (frozen on bcat data; not calibrated for TTE).
  - The TTE background of 372 unrelated short-burst pairs has a 99th percentile of 0.889.
  - **No candidate.** Nothing passes the M2 threshold, χ² and E_peak together.
  - Example of the energy test at work: bn210528586 / bn130701761 passes the shape χ² (p = 0.018) but has E_peak
    z = −14.9.
- **Not searched:** the 268–382 compatible partners per event without TTE. A common representation would need
  their TTE too.

**5.4 Answer and decision**
- **Does the new information separate same-source from unrelated signals better?**
  - TTE **time** structure: only marginally. It adds about 3 robust multi-pulse events per 30 bright short bursts.
  - **Energy:** yes on testable pairs, conditionally. Requiring the E_peak ratio to match the time stretch
    removes about 80 % of unrelated look-alikes, but the gain depends on an unvalidated systematic and covers only
    about half the pairs.
- **Real candidates:** none.
- **Expand TTE to the remaining bright short bursts?** **Not recommended now.**
  - Expected yield: about 15 robust events (95 % interval about 5–38) out of 147.
  - Each would bring about 15–37 multi-pulse partners that also need TTE. That is a few hundred pairs, about
    0.05 % of the 747,709 already screened.
  - Cost: about 2.5–4 GB streamed, plus a TTE-specific false-alarm calibration this sample cannot supply.
- **Higher value:** apply the catalogue E_peak constraint to the existing 747,709-pair screen and its 64 flagged
  pairs, which costs no download, after validating the E_peak systematic.
