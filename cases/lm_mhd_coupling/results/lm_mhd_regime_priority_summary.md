# LM-MHD B-U regime priority map

This reduced map scans representative liquid-metal properties over
magnetic-field and velocity ranges. It is a triage calculation, not a
resolved MUG or TokaMaker result.

Assumptions:
- length = 5.000e-02 m
- conducting thickness = 2.000e-02 m
- feedback reference field = 5.000e+00 T
- sheet-current closure factor = 2.500e-01

Representative row nearest `U = 0.2 m/s`, `B = 5 T`:

| material | Ha | Re | N | Rm | feedback | class |
|---|---:|---:|---:|---:|---:|---|
| PbLi | `5.465e+03` | `5.430e+04` | `5.499e+02` | `9.827e-03` | `4.966e-04` | `drag-important` |
| Li | `2.131e+04` | `1.086e+04` | `4.180e+04` | `3.931e-02` | `1.986e-03` | `two-way-candidate` |

Classification fractions over the scanned B-U grid:

| material | low | drag-important | two-way-candidate | full-induction-candidate |
|---|---:|---:|---:|---:|
| PbLi | `0.000` | `0.902` | `0.098` | `0.000` |
| Li | `0.000` | `0.655` | `0.195` | `0.150` |

Generated files:
- `lm_mhd_regime_priority_map.csv`
- `lm_mhd_regime_priority_summary.csv`
- `lm_mhd_regime_priority_metrics.png`
- `lm_mhd_regime_priority_classification.png`

Red-team interpretation:
- `N` flags Lorentz-force/drag importance; it says nothing by itself
  about magnetic feedback into the plasma.
- `Rm` flags when inductionless assumptions need scrutiny.
- The feedback panel is a sheet-current scale with a closure factor,
  not a solved current-closure or equilibrium perturbation.
- These maps prioritize what to simulate next; they do not validate
  MUG or claim blanket design performance.
