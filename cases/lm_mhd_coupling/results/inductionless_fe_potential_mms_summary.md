# Q1 finite-element inductionless potential MMS

This standalone Python reference assembles the weak form
`int sigma grad(w).grad(phi) = int sigma grad(w).(u x B)` with
all-insulating natural current closure and a nodal mean-value gauge.
It does not run MUG.

Convergence:

| n | max error | observed order |
|---:|---:|---:|
| 9 | `9.333e-06` |  |
| 17 | `2.338e-06` | `1.997` |
| 25 | `1.040e-06` | `1.999` |
| 33 | `5.849e-07` | `1.999` |
| 49 | `2.600e-07` | `1.999` |

Finest-grid diagnostic:
- max interior |div J| from nodal postprocessing = 7.546e-13.
- max wall-normal current from nodal postprocessing = 3.617e-07.
- max boundary-face |J.n| from FE flux recovery = 5.827e-04.
- The recovered boundary-face current decreases approximately first order
  over this sweep, as expected for simple Q1 element-gradient recovery.

Red-team interpretation:
- This checks weak-form sign, gauge, and convergence on a structured Q1 mesh.
- The boundary-face flux is recovered from FE element gradients at Gauss
  points, but no interface flux recovery is implemented yet.
- The nodal current and wall-current diagnostics are still postprocessing
  checks rather than the finite-element residual operator.
- No Fortran field structure, conducting-wall Robin term, tensor conductivity,
  or axisymmetric R metric is implemented here.
