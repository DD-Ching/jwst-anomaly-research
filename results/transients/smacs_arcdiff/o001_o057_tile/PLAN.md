# SMACS ERO o001 → VENUS o057: tile-median registration check (plan, written before the run)

Pending item from #288 (`results/transients/fade_calib/RESULT.md`, "Effect on other searches"): this pair's
F150W whole-box registration is (2.11, −1.82) px and was never compared with a tile median
(`arc_diff.py --reg-tile-arcsec 6`, the El Gordo fix).

- **Hypothesis:** the whole-box shift is pulled by a few tiles (bright cores, the Atalanta host light), as in
  El Gordo; the tile median differs by ≥ 0.5 px (ASSUMPTION threshold).
- **Ordinary outcome:** the two agree within 0.5 px; the #237/RESULT.md screen stands as it is.
- **Run:** the stored pair-1 command unchanged except `--reg-tile-arcsec 6` and `--out` here.
- **What changes the plan:** shift difference ≥ 0.5 px → compare the hit lists (frozen `classify` rule in
  `scripts/smacs_arcdiff_images.py`): Atalanta A/B must stay `known` (positive control); any new `vet` row is
  vetted; any lost row is reported. Difference < 0.5 px → record and close the item.
