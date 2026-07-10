# LM-MHD CouplingSnapshot audit tool

Prepared on branch `lm-mhd-upgrades`, 2026-07-02.

## Scope

`cases/lm_mhd_coupling/coupling_snapshot_audit.py` is an offline audit tool for
LM-MHD current handoff files using the existing `CouplingSnapshot` HDF5 schema.
It is not a MUG run, not a TokaMaker run, and not a coupled equilibrium update.
Its job is to make future handoff files easier to red-team before they are used
in a two-way plasma/liquid-metal loop.

## What It Checks

For any valid snapshot with `cell_volume`, the tool computes:

- total and absolute toroidal current using `J_phi dA`,
- current centroid in the poloidal plane,
- circular-loop magnetic feedback at configured probe points,
- feedback classification using the existing planning thresholds,
- max `|J_phi|`, max `|J_poloidal|`, and max `|J|`.

If the snapshot points form a complete structured `R-Z` grid, it also computes:

- max relative mismatch between `cell_volume` and `2*pi*R*dR*dZ`,
- max interior axisymmetric `|div J|`,
- max pointwise wall-normal current,
- signed outward boundary-current integrals,
- net integrated boundary current.

The structured-grid extraction is implemented as the tested helper
`structured_snapshot_current_grid(snapshot)`. It sorts shuffled snapshot points
by increasing `Z`, then increasing `R`, and rejects incomplete rectangular
point sets.

## Default Run

Command:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 cases/lm_mhd_coupling/coupling_snapshot_audit.py
```

Default input:

```text
cases/lm_mhd_coupling/results/axisymmetric_loop_voltage_feedback_snapshot.h5
```

Generated artifacts:

- `cases/lm_mhd_coupling/results/coupling_snapshot_audit_metrics.csv`
- `cases/lm_mhd_coupling/results/coupling_snapshot_audit_current.png`
- `cases/lm_mhd_coupling/results/coupling_snapshot_audit_summary.md`

Measured default metrics:

| quantity | value |
|---|---:|
| points | `5265` |
| structured grid | `81 x 65` |
| total toroidal current | `1.756e4 A` |
| max `|delta B|/|B_ref|` | `1.909e-3` |
| feedback class | `two-way-candidate` |
| max cell-volume relative error | `0.000e0` |
| max interior axisymmetric `|div J|` | `1.858e-10` |
| max integrated boundary current | `0.000e0 A` |
| net integrated boundary current | `0.000e0 A` |

## Red-Team Notes

- The magnetic-feedback calculation is still a circular-loop reduced estimate,
  not a TokaMaker Green's function solve.
- The boundary-current closure audit is only available when the snapshot points
  form a complete structured `R-Z` grid.
- Passing the audit means the HDF5 file is internally consistent with the
  reduced postprocessing conventions. It does not validate MUG physics or prove
  that the represented current field is physically correct.
- A nonzero cell-volume mismatch would directly contaminate `J_phi dA`
  feedback integrals, so that check should be treated as a hard handoff audit.
- Future MUG-exported snapshots should include `cell_volume`; otherwise the
  axisymmetric `J_phi dA` feedback integral is under-specified.

## Verification

Focused coupling tests after adding `structured_snapshot_current_grid`:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src/python pytest -q \
  src/tests/physics/test_lm_mhd_coupling.py
```

Result: `8 passed`.

The full focused LM-MHD Python suite after the audit continuation:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src/python pytest -q \
  src/tests/physics/test_lm_mhd_coupling.py \
  src/tests/physics/test_lm_mhd_closures.py \
  src/tests/physics/test_lm_mhd_feedback.py \
  src/tests/physics/test_lm_mhd_inductionless.py \
  src/tests/physics/test_lm_mhd_full_induction_bridge.py
```

Current expanded focused result after later power-identity and full-induction
bridge additions: `53 passed`.
