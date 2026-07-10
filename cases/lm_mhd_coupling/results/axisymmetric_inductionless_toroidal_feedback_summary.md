# Axisymmetric inductionless toroidal-feedback diagnostic

This manufactured diagnostic reconstructs axisymmetric inductionless
currents with an explicit toroidal electric field and estimates magnetic
feedback from the resulting `J_phi` through the existing
`CouplingSnapshot` circular-loop machinery. It does not run MUG or
TokaMaker.

Selected-case metrics:
- E_phi = 2.000e-02 V/m.
- total toroidal current = 1.599e+04 A.
- max |delta B|/|B_ref| = 1.656e-03.
- feedback classification = `two-way-candidate`.
- max |J_phi| = 1.696e+04 A/m^2.
- max |J x B| = 3.917e+03 N/m^3.
- max Joule heating = 3.524e+02 W/m^3.
- integrated Joule power = 3.022e+03 W.
- integrated mechanical power f.u = 2.969e+00 W.
- integrated electric power J.E = 3.025e+03 W.
- integrated power-balance residual = -2.659e-15 W.
- max pointwise power-balance residual = 1.705e-13 W/m^3.
- net axisymmetric torque = 9.311e+00 N m.
- max interior axisymmetric |div J| = 1.858e-10.
- max wall-normal current = 0.000e+00.
- max integrated boundary current = 0.000e+00 A.
- net integrated boundary current = 0.000e+00 A.

Generated files:
- `axisymmetric_inductionless_toroidal_feedback_scan.csv`
- `axisymmetric_inductionless_toroidal_feedback_maps.png`
- `axisymmetric_inductionless_toroidal_feedback_scan.png`
- `axisymmetric_inductionless_toroidal_feedback_snapshot.h5`

Red-team interpretation:
- Manufactured structured-grid reference only.
- The electric potential is prescribed, not solved from MUG or TokaMaker.
- Circular-loop feedback is a triage estimate, not a Grad-Shafranov
  equilibrium update.
- The toroidal electric field is imposed explicitly; future coupled work
  must derive it from loop voltage or boundary conditions.
- Poloidal current is manufactured to satisfy the axisymmetric divergence
  operator; this is not a resolved blanket-flow prediction.
