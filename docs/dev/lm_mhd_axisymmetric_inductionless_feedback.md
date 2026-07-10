# Axisymmetric inductionless toroidal-feedback diagnostic

Prepared on branch `lm-mhd-upgrades`, 2026-07-02.

## Scope

This note records a manufactured Python reference diagnostic for the current
path that matters most for two-way plasma/liquid-metal coupling: axisymmetric
toroidal current in the liquid metal and its magnetic feedback. It is not a
MUG solve, not a TokaMaker solve, and not a coupled equilibrium update.

## Added reference utilities

The inductionless helper module now includes axisymmetric structured-grid
diagnostics:

- `axisymmetric_current_from_flux_function(psi, r, z)` with
  `J_R = (1/R) dpsi/dZ`, `J_Z = -(1/R) dpsi/dR`.
- `axisymmetric_charge_conservation_residual(J, r, z)` for
  `(1/R) d(R J_R)/dR + dJ_Z/dZ`.
- `axisymmetric_wall_normal_current_extrema(J)`.
- `axisymmetric_boundary_current_integrals(J, r, z)` and
  `axisymmetric_net_boundary_current(J, r, z)` for signed outward current
  leakage in amperes using the axisymmetric surface metrics.
- `axisymmetric_cell_volumes(r, z)` for `2*pi*R*dR*dZ` control-volume weights.
- `mechanical_power_density(f, u)` for local `f.u`.
- `axisymmetric_torque_density(f, r, z)` for torque-about-axis density
  `R f_phi`.
- `axisymmetric_volume_integral(q, r, z)`,
  `axisymmetric_integrated_power(q, r, z)`, and
  `axisymmetric_integrated_torque(f, r, z)`.
- `reconstruct_axisymmetric_inductionless_current(...)`, using
  `J = sigma(-grad(phi) + u x B + E_phi e_phi)`.

The explicit `E_phi` term is important because loop voltage or toroidal
electric field can drive `J_phi`, and `J_phi` is the component that produces
poloidal magnetic-field perturbations in an axisymmetric feedback diagnostic.

## Diagnostic script

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 cases/lm_mhd_coupling/axisymmetric_inductionless_toroidal_feedback.py
```

Generated artifacts:

- `cases/lm_mhd_coupling/results/axisymmetric_inductionless_toroidal_feedback_scan.csv`
- `cases/lm_mhd_coupling/results/axisymmetric_inductionless_toroidal_feedback_maps.png`
- `cases/lm_mhd_coupling/results/axisymmetric_inductionless_toroidal_feedback_scan.png`
- `cases/lm_mhd_coupling/results/axisymmetric_inductionless_toroidal_feedback_snapshot.h5`
- `cases/lm_mhd_coupling/results/axisymmetric_inductionless_toroidal_feedback_summary.md`

## Headline result

The default selected case uses `E_phi = 0.02 V/m`, representative
`sigma0 = 8e5 S/m`, and `B_phi = 5 T`.

| quantity | value |
|---|---:|
| total toroidal current | `1.599e4 A` |
| max `|delta B|/|B_ref|` | `1.656e-3` |
| feedback class | `two-way-candidate` |
| max `|J_phi|` | `1.696e4 A/m^2` |
| max `|J x B|` | `3.917e3 N/m^3` |
| max Joule heating | `3.524e2 W/m^3` |
| integrated Joule power | `3.022e3 W` |
| integrated mechanical power `f.u` | `2.969e0 W` |
| net axisymmetric torque | `9.311e0 N m` |
| max interior axisymmetric `|div J|` | `1.858e-10` |
| max wall-normal current | `0.000e0` |
| max integrated boundary current | `0.000e0 A` |
| net integrated boundary current | `0.000e0 A` |

The scan shows the intended regime transition:

| `E_phi` [V/m] | max `|delta B|/|B_ref|` | class |
|---:|---:|---|
| `0.000` | `2.238e-5` | `weak-feedback` |
| `0.002` | `1.682e-4` | `weak-feedback` |
| `0.005` | `4.161e-4` | `weak-feedback` |
| `0.010` | `8.294e-4` | `weak-feedback` |
| `0.020` | `1.656e-3` | `two-way-candidate` |
| `0.050` | `4.136e-3` | `two-way-candidate` |

## Loop-voltage continuation

The script `cases/lm_mhd_coupling/axisymmetric_loop_voltage_feedback.py`
derives the toroidal electric field from a loop voltage:

```text
E_phi(R) = V_loop / (2*pi*R)
```

This is a more tokamak-facing parameterization than a constant `E_phi`. The
default selected case uses `V_loop = 0.2 V`, which gives
`E_phi = 1.592e-2` to `3.183e-2 V/m` over the diagnostic annulus.

| quantity | value |
|---|---:|
| total toroidal current | `1.756e4 A` |
| max `|delta B|/|B_ref|` | `1.909e-3` |
| feedback class | `two-way-candidate` |
| max `|J_phi|` | `2.469e4 A/m^2` |
| max `|J x B|` | `4.977e3 N/m^3` |
| max Joule heating | `7.972e2 W/m^3` |
| integrated Joule power | `3.513e3 W` |
| integrated mechanical power `f.u` | `-4.551e-1 W` |
| net axisymmetric torque | `9.311e0 N m` |
| max interior axisymmetric `|div J|` | `1.858e-10` |
| max wall-normal current | `0.000e0` |
| max integrated boundary current | `0.000e0 A` |
| net integrated boundary current | `0.000e0 A` |

Generated artifacts:

- `cases/lm_mhd_coupling/results/axisymmetric_loop_voltage_feedback_scan.csv`
- `cases/lm_mhd_coupling/results/axisymmetric_loop_voltage_feedback_maps.png`
- `cases/lm_mhd_coupling/results/axisymmetric_loop_voltage_feedback_scan.png`
- `cases/lm_mhd_coupling/results/axisymmetric_loop_voltage_feedback_snapshot.h5`
- `cases/lm_mhd_coupling/results/axisymmetric_loop_voltage_feedback_summary.md`

## Threshold continuation

The script `cases/lm_mhd_coupling/axisymmetric_loop_voltage_threshold_scan.py`
scans loop voltage and conductivity to estimate where the manufactured
feedback metric crosses the planning thresholds.

For the PbLi-like row nearest `sigma0 = 8e5 S/m`:

| threshold | estimated loop voltage |
|---|---:|
| two-way candidate, `|delta B|/B_ref = 1e-3` | `1.046e-1 V` |
| strong candidate, `|delta B|/B_ref = 1e-2` | not reached in scan |
| max feedback at `V_loop = 1 V` | `9.534e-3` |

Generated artifacts:

- `cases/lm_mhd_coupling/results/axisymmetric_loop_voltage_threshold_scan.csv`
- `cases/lm_mhd_coupling/results/axisymmetric_loop_voltage_thresholds.csv`
- `cases/lm_mhd_coupling/results/axisymmetric_loop_voltage_threshold_heatmap.png`
- `cases/lm_mhd_coupling/results/axisymmetric_loop_voltage_thresholds.png`
- `cases/lm_mhd_coupling/results/axisymmetric_loop_voltage_threshold_summary.md`

## Red-team interpretation

- This is a manufactured structured-grid diagnostic, not a blanket prediction.
- The scalar potential is prescribed; it is not solved by MUG.
- The constant-`E_phi` diagnostic imposes the toroidal electric field directly;
  the continuation derives it from loop voltage, but still prescribes the
  potential and velocity fields.
- The magnetic feedback estimate uses circular-loop superposition; it is a
  triage metric, not a Grad-Shafranov equilibrium response.
- The global power and torque integrals are only as meaningful as the
  manufactured fields. The selected constant-`E_phi` and loop-voltage cases
  share the same net torque because the manufactured poloidal-current piece
  dominates that component, while the integrated mechanical power changes sign.
- The new boundary-current integral is a conservative closure audit. It is
  separate from the structured-gradient `div J` residual and the pointwise
  wall-normal current extrema.
- The result supports the priority ranking: two-way feedback is not always
  dominant, but plausible toroidal-current drives can push the blanket into the
  `two-way-candidate` regime, so any serious coupled LM/plasma workflow needs
  this diagnostic path.

## Verification

Focused test command after adding the axisymmetric utilities:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src/python pytest -q src/tests/physics/test_lm_mhd_inductionless.py
```

Earlier result before global integral utilities: `22 passed`.
Current inductionless-only result after adding global integral utilities:
`26 passed`.
Current expanded focused LM-MHD Python result after adding power-identity and
full-induction bridge checks: `53 passed`.

The diagnostic script was rerun, HDF5 snapshot roundtrip shape was checked, and
the saved PNGs for the constant-`E_phi`, loop-voltage, and threshold scans were
confirmed nonblank with PIL.
