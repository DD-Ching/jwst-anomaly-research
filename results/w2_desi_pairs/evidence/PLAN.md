# W2 DESI pairs, follow-up: identity, literature evidence and follow-up value of the 13 rejected pairs (bounded audit)

Owner brief 2026-10-23 after #472 (start 01:48 UTC):
- **Budget:** 2 h. No new schedule, no outside contact, no account.
- **Scope:** reuses #472's `pairs.json`, results and sources. No spectrum is downloaded again, and the failed
  classifier is not rerun. #472's numbers, frozen threshold and failed verdict stay as they are.
- **Not done:** no new classifier, statistic, threshold or PSF tool.

## Steps
1. **Evidence table for all 13 pairs** (`evidence.json`, `EVIDENCE.md`):
   - names, aliases and both components' coordinates and DESI TARGETIDs (zpix metadata only), plus any pairing
     ambiguity;
   - the original catalogue class (pinned D-064 inputs) and later published classes, each with its basis and limits;
   - whether a source used the same DESI DR1 spectra;
   - independent imaging, resolved spectra, time delays or lens-galaxy evidence;
   - the observation still missing that could change the call.
2. **Priority checks** (the others are not skipped):
   - J2244−0550: Lemon et al. 2023 §4.2.16 calls it a likely lensed quasar; trace the "UQP" label to its source and
     version.
   - J0130+0725: Chen et al. (arXiv:2607.00093) used DESI DR1; keep "dual favoured, lensing not excluded".
   - J0728+2607: Lemon et al. NOT spectra.
   - The two pairs that overlap the Hennawi controls, and the "Binary QSO" labels of J1511 and J1552: trace each
     label to its original evidence.
3. **Over-reading corrections in #472** (wording only):
   - D ≤ T means "no difference above threshold", not "same source";
   - 9 / 9 is not an accuracy for "same source".
4. **Up to 3 genuinely open pairs:** chosen by which independent evidence exists and what is missing, never by #472's
   failed score.
   - For each, check official product catalogues and published measurements for public data not yet used that could
     fill the gap: product id, epoch, band, resolution and access.
   - Not found means "not confirmed", never "does not exist". If only the same data could be re-analysed, say that
     no new observational information exists.
5. **Interpretation rules:**
   - no lens galaxy seen ≠ no ordinary lens;
   - similar spectra ≠ one source;
   - different spectra need time delay, variability, microlensing and blending considered before "binary";
   - conflicting or hedged literature stays uncertain.
