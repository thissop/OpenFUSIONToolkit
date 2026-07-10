# Axisymmetric Q1 FE inductionless potential MMS

This standalone Python reference assembles the weak form with the
axisymmetric R metric. It does not run MUG or TokaMaker.

| n | max error | observed order |
|---:|---:|---:|
| 9 | `1.943e-05` |  |
| 17 | `4.886e-06` | `1.992` |
| 25 | `2.175e-06` | `1.995` |
| 33 | `1.225e-06` | `1.997` |
| 49 | `5.448e-07` | `1.998` |

Finest-grid nodal postprocessing:
- max axisymmetric |div J| = 3.123e-07.
- max wall-normal current = 1.356e-06.

Red-team interpretation:
- This checks only the scalar axisymmetric R-weighted weak form.
- It does not include toroidal electric field, tensor conductivity,
  TokaMaker feedback, or MUG field extraction.
- Current diagnostics are structured-grid nodal postprocessing, not
  conservative FE flux recovery.
