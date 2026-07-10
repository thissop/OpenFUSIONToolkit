# Axisymmetric loop-voltage feedback diagnostic

This manufactured diagnostic derives `E_phi(R) = V_loop/(2*pi*R)`,
reconstructs inductionless current, and estimates toroidal-current
magnetic feedback. It does not run MUG or TokaMaker.

Selected-case metrics:
- V_loop = 2.000e-01 V.
- E_phi range = [1.592e-02, 3.183e-02] V/m.
- total toroidal current = 1.756e+04 A.
- max |delta B|/|B_ref| = 1.909e-03.
- feedback classification = `two-way-candidate`.
- max |J_phi| = 2.469e+04 A/m^2.
- max |J x B| = 4.977e+03 N/m^3.
- max Joule heating = 7.972e+02 W/m^3.
- integrated Joule power = 3.513e+03 W.
- integrated mechanical power f.u = -4.551e-01 W.
- integrated electric power J.E = 3.513e+03 W.
- integrated power-balance residual = -1.153e-14 W.
- max pointwise power-balance residual = 3.411e-13 W/m^3.
- net axisymmetric torque = 9.311e+00 N m.
- max interior axisymmetric |div J| = 1.858e-10.
- max wall-normal current = 0.000e+00.
- max integrated boundary current = 0.000e+00 A.
- net integrated boundary current = 0.000e+00 A.

Generated files:
- `axisymmetric_loop_voltage_feedback_scan.csv`
- `axisymmetric_loop_voltage_feedback_maps.png`
- `axisymmetric_loop_voltage_feedback_scan.png`
- `axisymmetric_loop_voltage_feedback_snapshot.h5`

Red-team interpretation:
- Manufactured structured-grid reference only.
- Loop voltage sets the toroidal electric field, but the potential and
  velocity fields are still prescribed.
- Circular-loop feedback is a triage estimate, not a Grad-Shafranov
  equilibrium update.
