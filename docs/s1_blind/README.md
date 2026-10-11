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
- **Unrelated pairs pass the consistency test at 39 %** (top 100: 47 %). A full single-pulse screen (517,629
  pairs, about 52 above threshold) would give **about 20 false candidates**.
- **Deformed injections:** 46 / 960 recovered, **4.8 %** total. Undeformed: 4 / 120.
- **Decision (pre-declared):** both criteria fail (false candidates ≤ 3; deformed recovery ≥ 10 %).
  - Once stretch, lag and per-band scale are free, one smooth pulse is matched by many unrelated single pulses.
  - The blind spot is **structural, not a data-access problem**. More time resolution cannot add shape
    information that a single smooth pulse does not have.

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
- **Decision (pre-declared ≥ 25 % gain):** 17 % falls short, so the full TTE run was not started.
  - Extrapolation: about 25 of the 147 bright short bursts would become multi-pulse, and fewer of the fainter ones.
  - At ms resolution these could only be paired with each other (|Δk| ≤ 2 versus long bursts). That adds a few
    hundred pairs at most.
  - **Cost** (measured): 497 MB and 88 s for 30 bursts, so about 0.6 h and about 10 GB streamed for all 699.
    It is cheap, but the yield is low.

## 4. Answer

- **Single-pulse bursts (2,129):** not searchable for deformed twins with this test. Unrelated pairs pass too
  often (39 %), whatever the data product.
- **Faint bursts with no resolved pulse (846):** TTE does not add photons, so they stay outside a shape-based
  search.
- **Short bursts (699):** TTE resolves extra pulses in only about 17 % of the bright ones (calibrated against a
  control at 87 %). A full run is affordable but would add few usable pairs. It was not started, per the
  pre-declared rule.
- These are blind spots of shape matching itself. Reaching them needs **another observable**, not more of the same
  light-curve data. Examples:
  - spectral-lag and E_peak(t) relations (CSPEC/TTE spectral fits) as the matching quantity;
  - a different messenger.
  This is a direction, not a result.
