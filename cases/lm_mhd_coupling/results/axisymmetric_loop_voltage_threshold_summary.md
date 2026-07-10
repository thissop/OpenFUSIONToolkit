# Axisymmetric loop-voltage/conductivity threshold scan

This reduced manufactured scan estimates when loop-voltage-driven
axisymmetric blanket currents cross magnetic-feedback planning
thresholds. It does not run MUG or TokaMaker.

PbLi-like row nearest `sigma0 = 8e5 S/m`:
- two-way threshold (`1e-3`) V_loop = 1.046e-01 V.
- strong threshold (`1e-2`) V_loop = not reached.
- max feedback at max scanned voltage = 9.534e-03.

Threshold table:

| sigma0 [S/m] | V_loop at 1e-3 [V] | V_loop at 1e-2 [V] | max feedback |
|---:|---:|---:|---:|
| `2.000e+05` | `4.194e-01` | not reached | `2.383e-03` |
| `5.000e+05` | `1.676e-01` | not reached | `5.958e-03` |
| `8.000e+05` | `1.046e-01` | not reached | `9.534e-03` |
| `1.200e+06` | `6.967e-02` | `6.992e-01` | `1.430e-02` |
| `2.000e+06` | `4.133e-02` | `4.194e-01` | `2.383e-02` |

Generated files:
- `axisymmetric_loop_voltage_threshold_scan.csv`
- `axisymmetric_loop_voltage_thresholds.csv`
- `axisymmetric_loop_voltage_threshold_heatmap.png`
- `axisymmetric_loop_voltage_thresholds.png`

Red-team interpretation:
- Thresholds are interpolation over a coarse manufactured scan.
- The current closure, geometry, and field feedback are reduced models.
- Use these numbers only to prioritize coupled sensitivity studies, not
  as blanket design predictions.
