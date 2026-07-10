# CouplingSnapshot audit

This offline diagnostic reads a `CouplingSnapshot` HDF5 file and
computes reduced current-feedback and current-closure checks. It does
not run MUG or TokaMaker.

Headline metrics:
- snapshot = `cases/lm_mhd_coupling/results/axisymmetric_loop_voltage_feedback_snapshot.h5`
- points = 5265
- total toroidal current = 1.756e+04 A
- absolute toroidal current = 1.756e+04 A
- max |delta B|/|B_ref| = 1.909e-03
- classification = `two-way-candidate`
- max |J_phi| = 2.469e+04 A/m^2
- structured grid = `true`
- grid = 81 x 65
- max cell-volume relative error = 0.000e+00
- max interior axisymmetric |div J| = 1.858e-10
- max wall-normal current = 0.000e+00 A/m^2
- max integrated boundary current = 0.000e+00 A
- net integrated boundary current = 0.000e+00 A

Generated files:
- `coupling_snapshot_audit_metrics.csv`
- `coupling_snapshot_audit_current.png`

Red-team interpretation:
- Magnetic feedback is a circular-loop triage estimate.
- Boundary-current checks are only available for complete structured
  R-Z snapshots.
- A clean audit does not validate MUG physics; it verifies this HDF5
  handoff file against reduced postprocessing conventions.
