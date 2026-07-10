# Q1 FE two-region interface-flux MMS

This standalone Python reference aligns a conductivity jump with Q1
element faces, solves the insulating weak form, and recovers normal
current on both sides of the material interface. It does not run MUG.

Potential and interface-current recovery are at roundoff because the
piecewise-linear exact potential is in the conforming Q1 space and the
interface is aligned to element faces.

| n | max potential error | max interface target error | max interface jump |
|---:|---:|---:|---:|
| 17 | `5.551e-16` | `4.270e-15` | `3.109e-15` |
| 33 | `2.554e-15` | `1.271e-14` | `7.661e-15` |
| 49 | `1.238e-14` | `6.240e-14` | `8.882e-15` |

Finest-grid interface current:
- max interface target-current error = 6.240e-14.
- max left/right interface-current jump = 8.882e-15.

Red-team interpretation:
- Interface is vertical, conforming, and exactly aligned to element faces.
- Conductivity is scalar and piecewise constant.
- This checks interface flux recovery in a standalone Python Q1 prototype,
  not in MUG and not in an axisymmetric metric.
