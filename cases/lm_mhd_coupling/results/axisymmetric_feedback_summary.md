# Axisymmetric current-feedback diagnostic

This is an analytic/reduced-model coupling diagnostic, not a MUG or TokaMaker result.
A manufactured toroidal liquid-metal current patch is converted to circular-loop
current elements and evaluated at plasma-side probe points with the standard
axisymmetric Biot-Savart formula.

Baseline result:
- Peak J_phi = 1.000e+05 A/m^2.
- Integrated |J_phi| dA = 2.163e+04 A.
- Max probe |delta B| = 1.288e-02 T.
- Max |delta B|/|B_ref| = 2.576e-03.
- Feedback classification = two-way-candidate.

Generated files:
- `axisymmetric_feedback_snapshot.h5`: HDF5 roundtrip snapshot for the baseline current patch.
- `axisymmetric_feedback_scan.csv`: current-amplitude sweep and feedback metrics.
- `axisymmetric_feedback_field_map.png`: current patch, probe feedback map, and midplane line.
- `axisymmetric_feedback_scaling.png`: linearity and threshold plot.

Red-team interpretation:
- The field estimate includes geometric 1/R loop effects, but no vessel, wall-current,
  iron, plasma-response, or TokaMaker equilibrium solve.
- The sampled currents are manufactured and do not satisfy an electric-potential closure.
- Probes are placed outside the current patch to avoid unresolved self-field singularities.
- Crossing 1e-3 in |delta B|/|B_ref| is treated as a prompt for a coupled sensitivity
  study, not as a prediction of the plasma-equilibrium displacement.
