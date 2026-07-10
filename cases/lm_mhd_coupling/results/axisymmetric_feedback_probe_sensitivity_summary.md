# Axisymmetric feedback probe-sensitivity scan

This red-team diagnostic varies probe placement and circular-loop
regularization for a fixed `CouplingSnapshot`. It does not run MUG or
TokaMaker.

Default-like row nearest offset `0.1 m` and regularization `0.03 m`:
- side = `inboard`.
- max |delta B|/|B_ref| = 1.968e-03.
- classification = `two-way-candidate`.

Extremes across the scan:
- strongest feedback = 2.047e-03 (`inboard`, offset 0.05 m, eps 0.01 m).
- weakest feedback = 2.361e-04 (`outboard`, offset 0.5 m, eps 0.1 m).
- fraction of rows at or above `1e-3` = 0.500.
- fraction of rows at or above `1e-2` = 0.000.

Generated files:
- `axisymmetric_feedback_probe_sensitivity.csv`
- `axisymmetric_feedback_probe_sensitivity_heatmap.png`
- `axisymmetric_feedback_probe_sensitivity_curves.png`

Red-team interpretation:
- This scan tests sensitivity of the reduced circular-loop metric, not
  plasma equilibrium sensitivity.
- Rows close to the current annulus are intentionally near-field-biased.
- If only near-field rows cross `1e-3`, that weakens the two-way case.
  If many offset/regularization rows cross `1e-3`, that strengthens the
  priority for coupled sensitivity studies.
