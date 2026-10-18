# UDS x3y0: L2 (`_cal`) exposure check at the S peak (plan, 2026-10-18, written before any pixel is read)

**Question.** Is the sharp horizontal line through the x3y0 S peak (34.50675, −5.25976; F200W) present in only one
(or a few, consecutive) of the o001 `_cal` exposures that cover the peak?

**Hypothesis (ordinary first).** A satellite / scattered-light trail: present in one exposure (or one detector
read sequence), absent in the dithers of the same pointing. Alternatives: a detector artefact fixed in detector
coordinates (present in every exposure but moving on the sky with the dither), or an astrophysical / sky feature
(present in every exposure at the same sky position).

**Test (frozen).** For every o001 F200W `_cal` exposure whose footprint contains the peak: byte-range read the SCI
cutout (±10″) around the peak, convert the line's sky track (row band through the peak in the c1002 frame,
± 3 px) to each exposure's pixels with its WCS, and compare the median band surface brightness with two parallel
control bands ±3″ away. A line is "in" an exposure if band − control > 5 × the robust scatter of control rows.
- in exactly one exposure, or only in exposures from one dither position → **trail**: x3y0 explained (post hoc,
  labelled; not the frozen V1 ∪ V2 rule);
- in every exposure at the same sky position → not a trail; x3y0 stays a flag needing the full vet;
- fixed in detector pixels → detector artefact (explained).
Thresholds are ASSUMPTIONs.
