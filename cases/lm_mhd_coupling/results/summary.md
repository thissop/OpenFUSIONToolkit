# LM-MHD coupling regime scan

This is an executable triage artifact, not a MUG/TokaMaker result. It scans representative
liquid-metal properties and reports dimensionless groups plus an intentionally conservative
sheet-current feedback estimate.

Assumptions:
- Length scale L = 0.05 m.
- Conducting sheet thickness for feedback estimate = 0.05 m.
- Reference field for feedback ratio = 5 T.
- Closure factor multiplying sigma U B = 0.1. A value below one
  acknowledges that electric-potential closure can reduce or redirect the current.
- Classification thresholds: Rm >= 0.1 => question inductionless; feedback >= 1e-3 =>
  two-way-candidate; otherwise N >= 1 => drag-important.

| material | Ha range | Re range | Rm range | N range | feedback range | dominant labels |
|---|---:|---:|---:|---:|---:|---|
| PbLi | 5.41e+02-1.08e+04 | 2.78e+01-5.56e+05 | 5.03e-06-1.01e-01 | 5.26e-01-4.21e+06 | 2.51e-08-1.01e-02 | drag-important:6486, full-induction-candidate:80, low:23, two-way-candidate:611 |
| Li | 2.11e+03-4.22e+04 | 5.56e+00-1.11e+05 | 2.01e-05-4.02e-01 | 4.00e+01-3.20e+08 | 1.01e-07-4.02e-02 | drag-important:5521, full-induction-candidate:1040, two-way-candidate:639 |

Red-team interpretation:
- Large N with small Rm means one-way imposed-field LM-MHD can still have very strong
  Lorentz braking and unresolved Hartmann-layer pressure loss. This is the regime where
  effective drag/SM82-style closures are most valuable for full-device meshes.
- A feedback ratio above 1e-3 is a prompt to run a coupled field sensitivity study, not
  proof that the plasma equilibrium will change by that amount.
- The sheet-current estimate is an upper-bound scale. A real duct can have smaller or
  differently directed currents after solving electric potential and wall closure.

Generated files:
- `coupling_regime_scan.csv`
- `pbli_regime_heatmap.png`
- `li_regime_heatmap.png`
- `two_way_threshold_lines.png`
