# LM-MHD inductionless reference infrastructure

This note records the third autonomous LM-MHD upgrade slice from 2026-07-01.
It addresses the roadmap's highest-impact long-term upgrade, the low magnetic
Reynolds number electric-potential mode, but only at the Python reference and
diagnostic level. No Fortran solver physics was changed.

## What changed

Added `OpenFUSIONToolkit.LM_MHD.inductionless` with small tested utilities:

- `axisymmetric_cell_volumes(r, z)`
- `axisymmetric_charge_conservation_residual(current_density, r, z)`
- `axisymmetric_current_from_flux_function(flux_function, r, z)`
- `axisymmetric_wall_normal_current_extrema(current_density)`
- `ohms_law_current_density(sigma, velocity, magnetic_field, electric_potential_gradient)`
- `lorentz_force_density(current_density, magnetic_field)`
- `joule_heating_density(current_density, sigma)`
- `structured_gradient_2d(scalar, x, z)`
- `motional_electric_field(velocity, magnetic_field)`
- `potential_source_from_motional_emf(motional_emf, x, z)`
- `reconstruct_axisymmetric_inductionless_current(potential, velocity, magnetic_field, sigma, r, z, toroidal_electric_field)`
- `weighted_potential_source_from_motional_emf(conductivity, motional_emf, x, z)`
- `insulating_wall_normal_gradient(motional_emf, x, z)`
- `insulating_wall_normal_flux(conductivity, motional_emf, x, z)`
- `reconstruct_inductionless_current_2d(potential, velocity, magnetic_field, sigma, x, z)`
- `divergence_free_current_from_streamfunction(streamfunction, x, z)`
- `charge_conservation_residual_2d(current_density, x, z)`
- `wall_normal_current_extrema(current_density)`

The sign convention pinned down here is:

```text
J = sigma (-grad(phi) + u x B)
f = J x B
q = |J|^2 / sigma
```

That convention is the one future electric-potential work should red-team
against before adding a new finite-element field or Poisson solve.

## Manufactured diagnostic

The script `cases/lm_mhd_coupling/inductionless_mms_diagnostic.py` creates a
rectangular manufactured current using a streamfunction:

```text
J_x = d psi / dz
J_z = -d psi / dx
```

This makes `div J = 0` by construction under the same structured-grid
operators used by the diagnostic. A velocity field is then chosen so that, with
uniform `B_y` and zero electric potential, `sigma (u x B)` reconstructs that
current.

Generated artifacts in `cases/lm_mhd_coupling/results/`:

- `inductionless_mms_metrics.csv`
- `inductionless_mms_fields.png`
- `inductionless_mms_midplane.png`
- `inductionless_mms_summary.md`
- `inductionless_poisson_mms_metrics.csv`
- `inductionless_poisson_mms_solution.png`
- `inductionless_poisson_mms_convergence.png`
- `inductionless_poisson_mms_summary.md`
- `inductionless_neumann_mms_metrics.csv`
- `inductionless_neumann_mms_solution.png`
- `inductionless_neumann_mms_convergence.png`
- `inductionless_neumann_mms_summary.md`
- `inductionless_neumann_wall_mms_metrics.csv`
- `inductionless_neumann_wall_mms_solution.png`
- `inductionless_neumann_wall_mms_summary.md`
- `inductionless_full_mms_metrics.csv`
- `inductionless_full_mms_potential.png`
- `inductionless_full_mms_current.png`
- `inductionless_full_mms_midplane.png`
- `inductionless_full_mms_summary.md`
- `inductionless_variable_sigma_mms_metrics.csv`
- `inductionless_variable_sigma_mms_solution.png`
- `inductionless_variable_sigma_mms_convergence.png`
- `inductionless_variable_sigma_mms_summary.md`
- `inductionless_two_region_sigma_mms_metrics.csv`
- `inductionless_two_region_sigma_mms_solution.png`
- `inductionless_two_region_sigma_mms_flux.png`
- `inductionless_two_region_sigma_mms_summary.md`
- `inductionless_variable_sigma_neumann_mms_metrics.csv`
- `inductionless_variable_sigma_neumann_mms_solution.png`
- `inductionless_variable_sigma_neumann_mms_convergence.png`
- `inductionless_variable_sigma_neumann_mms_summary.md`
- `inductionless_variable_sigma_full_mms_metrics.csv`
- `inductionless_variable_sigma_full_mms_potential.png`
- `inductionless_variable_sigma_full_mms_current.png`
- `inductionless_variable_sigma_full_mms_summary.md`
- `inductionless_two_region_full_mms_metrics.csv`
- `inductionless_two_region_full_mms_potential.png`
- `inductionless_two_region_full_mms_current.png`
- `inductionless_two_region_full_mms_interface.png`
- `inductionless_two_region_full_mms_summary.md`
- `inductionless_postprocess_demo_input.csv`
- `inductionless_postprocess_demo_metrics.csv`
- `inductionless_postprocess_demo_current.png`
- `inductionless_postprocess_demo_summary.md`
- `inductionless_fe_potential_mms_metrics.csv`
- `inductionless_fe_potential_mms_solution.png`
- `inductionless_fe_potential_mms_convergence.png`
- `inductionless_fe_potential_mms_summary.md`
- `inductionless_fe_two_region_flux_mms_metrics.csv`
- `inductionless_fe_two_region_flux_mms_solution.png`
- `inductionless_fe_two_region_flux_mms_interface.png`
- `inductionless_fe_two_region_flux_mms_summary.md`
- `inductionless_axisymmetric_fe_potential_mms_metrics.csv`
- `inductionless_axisymmetric_fe_potential_mms_solution.png`
- `inductionless_axisymmetric_fe_potential_mms_convergence.png`
- `inductionless_axisymmetric_fe_potential_mms_summary.md`
- `mug_restart_manifest.csv`
- `mug_restart_field_ranges.png`
- `mug_restart_manifest_summary.md`
- `axisymmetric_inductionless_toroidal_feedback_scan.csv`
- `axisymmetric_inductionless_toroidal_feedback_maps.png`
- `axisymmetric_inductionless_toroidal_feedback_scan.png`
- `axisymmetric_inductionless_toroidal_feedback_snapshot.h5`
- `axisymmetric_inductionless_toroidal_feedback_summary.md`
- `axisymmetric_loop_voltage_feedback_scan.csv`
- `axisymmetric_loop_voltage_feedback_maps.png`
- `axisymmetric_loop_voltage_feedback_scan.png`
- `axisymmetric_loop_voltage_feedback_snapshot.h5`
- `axisymmetric_loop_voltage_feedback_summary.md`
- `axisymmetric_loop_voltage_threshold_scan.csv`
- `axisymmetric_loop_voltage_thresholds.csv`
- `axisymmetric_loop_voltage_threshold_heatmap.png`
- `axisymmetric_loop_voltage_thresholds.png`
- `axisymmetric_loop_voltage_threshold_summary.md`
- `coupling_snapshot_audit_metrics.csv`
- `coupling_snapshot_audit_current.png`
- `coupling_snapshot_audit_summary.md`
- `axisymmetric_feedback_probe_sensitivity.csv`
- `axisymmetric_feedback_probe_sensitivity_heatmap.png`
- `axisymmetric_feedback_probe_sensitivity_curves.png`
- `axisymmetric_feedback_probe_sensitivity_summary.md`
- `lm_mhd_regime_priority_map.csv`
- `lm_mhd_regime_priority_summary.csv`
- `lm_mhd_regime_priority_metrics.png`
- `lm_mhd_regime_priority_classification.png`
- `lm_mhd_regime_priority_summary.md`
- `lm_mhd_artifact_manifest.csv`
- `lm_mhd_artifact_manifest_summary.md`
- `full_induction_low_rm_bridge_scan.csv`
- `full_induction_low_rm_bridge_profiles.png`
- `full_induction_low_rm_bridge_scaling.png`
- `full_induction_low_rm_bridge_transient.png`
- `full_induction_low_rm_bridge_summary.md`

The current axisymmetric utility continuation also adds global integral
diagnostics:

- `mechanical_power_density(force_density, velocity)`
- `axisymmetric_torque_density(force_density, r, z)`
- `axisymmetric_volume_integral(density, r, z)`
- `axisymmetric_integrated_power(power_density, r, z)`
- `axisymmetric_integrated_torque(force_density, r, z)`
- `axisymmetric_boundary_current_integrals(current_density, r, z)`
- `axisymmetric_net_boundary_current(current_density, r, z)`

The selected constant-`E_phi` case now reports integrated Joule power
`3.022e3 W`, integrated mechanical power `2.969e0 W`, and net axisymmetric
torque `9.311e0 N m`. The selected loop-voltage case reports integrated Joule
power `3.513e3 W`, integrated mechanical power `-4.551e-1 W`, and the same
net axisymmetric torque `9.311e0 N m` for this manufactured field set. These
are reduced diagnostic integrals, not MUG force-balance or blanket-design
predictions.

Both selected cases also report zero signed integrated boundary-current
leakage in the regenerated CSV artifacts. This is a conservative current
closure audit using axisymmetric surface metrics, separate from the
structured-gradient `div J` residual.

The handoff-facing script `cases/lm_mhd_coupling/coupling_snapshot_audit.py`
now reads a `CouplingSnapshot` HDF5 file, computes the reduced magnetic
feedback metrics, infers a structured `R-Z` grid when possible, and applies the
same conservative boundary-current audit. On the default loop-voltage snapshot
it recovers `max |delta B|/|B_ref| = 1.909e-3`, classifies the file as
`two-way-candidate`, detects an `81 x 65` grid, and reports zero integrated
boundary leakage plus zero cell-volume relative error against
`2*pi*R*dR*dZ`.

The probe-sensitivity script
`cases/lm_mhd_coupling/axisymmetric_feedback_probe_sensitivity.py` red-teams
that same reduced feedback metric. Across inboard/outboard probe lines,
offsets, and loop regularizations, 50% of rows cross the `1e-3`
two-way-candidate threshold and no rows cross `1e-2`. The strongest row is
`2.047e-3`; the weakest row is `2.361e-4`.

The broader regime-priority script
`cases/lm_mhd_coupling/lm_mhd_regime_priority_map.py` scans representative
PbLi and Li over B-U space. With the default length, sheet thickness, and
closure factor, the representative PbLi point near `U = 0.2 m/s`, `B = 5 T`
is `drag-important` with feedback `4.966e-4`; the Li point is
`two-way-candidate` with feedback `1.986e-3`.

The artifact manifest script
`cases/lm_mhd_coupling/lm_mhd_artifact_manifest.py` records hashes, CSV row
counts, and PNG sanity metadata for the current coupling-diagnostic result set.
The latest manifest contains 32 entries: 9 CSV files, 14 PNG files, 2 HDF5
files, and 7 Markdown summaries; all 14 PNG files were nonblank.

The full-induction bridge script
`cases/lm_mhd_coupling/full_induction_low_rm_bridge.py` directly addresses the
model-choice question. In a prescribed 1D channel, the steady full-induction
Ampere current `(1/mu0) db_x/dz` matches the net-current-zero inductionless
Ohm's-law current. For the selected PbLi-like case (`U=0.2 m/s`,
`sigma=8e5 S/m`, `a=0.05 m`, `B0=5 T`), `Rm = 1.005e-2`,
`max |b_x|/B0 = 1.290e-3`, and the current bridge relative error is
`1.875e-5`.

Headline metrics from the default run:

| quantity | value |
|---|---:|
| sigma | `8.000e5 S/m` |
| B_y | `5.000 T` |
| max `|J|` | `1.000e5 A/m^2` |
| max `|u|` | `2.501e-2 m/s` |
| max interior `|div J|` | `2.328e-8 A/m^3` |
| max wall-normal current | `1.225e-11 A/m^2` |

## Poisson reference solve

The same module now includes `solve_potential_dirichlet_2d`, a uniform-grid
finite-difference solve for

```text
laplacian(phi) = source
```

with Dirichlet boundary values. The diagnostic
`cases/lm_mhd_coupling/inductionless_poisson_mms.py` verifies it against
`phi = sin(pi x) sin(pi z)`, for which `source = -2 pi^2 phi`. The observed
max-error convergence is second order:

| n | h | max error | observed order |
|---:|---:|---:|---:|
| 17 | `6.250e-2` | `3.219e-3` | |
| 25 | `4.167e-2` | `1.429e-3` | `2.003` |
| 33 | `3.125e-2` | `8.036e-4` | `2.001` |
| 49 | `2.083e-2` | `3.571e-4` | `2.001` |
| 65 | `1.562e-2` | `2.008e-4` | `2.000` |

This does not yet represent insulating MHD walls. It is a scalar solve target
for the eventual electric-potential equation and a useful way to keep sign and
operator conventions from drifting while the Fortran design is still pending.

## Homogeneous Neumann reference solve

The module also includes `solve_potential_neumann_2d`, a uniform-grid
finite-difference solve for the same scalar Poisson equation with homogeneous
Neumann walls and a mean-value gauge. The diagnostic
`cases/lm_mhd_coupling/inductionless_neumann_mms.py` verifies it against
`phi = cos(pi x) cos(pi z)`, which has zero normal derivative on all four
boundaries. The source has zero grid mean to roundoff, satisfying the Neumann
compatibility condition.

The observed max-error convergence is again second order:

| n | h | max mean-free error | observed order |
|---:|---:|---:|---:|
| 17 | `6.250e-2` | `3.219e-3` | |
| 25 | `4.167e-2` | `1.429e-3` | `2.003` |
| 33 | `3.125e-2` | `8.036e-4` | `2.001` |
| 49 | `2.083e-2` | `3.571e-4` | `2.001` |
| 65 | `1.562e-2` | `2.008e-4` | `2.000` |

This is closer to insulating-wall electric-potential checks than the Dirichlet
solve, but homogeneous wall data alone is not enough for a moving liquid metal.

## Nonzero Neumann wall-data reference

The Neumann reference now also accepts prescribed outward-normal derivative
data on `x_min`, `x_max`, `z_min`, and `z_max`. The diagnostic
`cases/lm_mhd_coupling/inductionless_neumann_wall_mms.py` verifies the wall-data
plumbing with a quadratic manufactured solution:

```text
phi = (x - 1/2)^2 + (z - 1/2)^2 - mean(phi)
laplacian(phi) = 4
d phi / dn = 1 on all four walls
```

The default `n = 65` run recovers the manufactured solution to roundoff:

| quantity | value |
|---|---:|
| max error | `3.220e-15` |
| rms error | `1.533e-15` |
| solved mean(phi) | `-1.245e-16` |

This is the scalar boundary-condition mechanism needed before testing the
insulating inductionless condition `d phi / dn = (u x B) . n`.

## Full constant-conductivity inductionless MMS

The diagnostic `cases/lm_mhd_coupling/inductionless_full_mms.py` now connects
the reference pieces:

```text
source = div(u x B)
d phi / dn = (u x B) . n
J = sigma (-grad(phi) + u x B)
```

The manufactured construction sets

```text
u x B = grad(phi_exact) + J_target / sigma
```

where `J_target` is generated from a streamfunction, so it is divergence-free
and has zero wall-normal current. The potential solve should recover
`phi_exact`, and the reconstructed current should recover `J_target`.

Default metrics:

| quantity | value |
|---|---:|
| max `|phi - phi_exact|` | `1.102e-16` |
| max `|J - J_target|` | `3.553e-10 A/m^2` |
| max interior `|div J|` | `2.235e-8 A/m^3` |
| max wall-normal current | `3.553e-10 A/m^2` |

This is the strongest Python-side reference target so far for a future
constant-conductivity inductionless implementation. It still does not include
variable conductivity, material interfaces, conducting-wall Robin closure, or a
finite-element weak form.

## Smooth variable-conductivity reference

The diagnostic `cases/lm_mhd_coupling/inductionless_variable_sigma_mms.py`
verifies coefficient handling for

```text
div(sigma grad(phi)) = source
```

with

```text
sigma = 1 + 0.3 x + 0.2 z
phi = sin(pi x) sin(pi z)
```

and analytic `source = div(sigma grad(phi))`. It uses Dirichlet walls to isolate
the coefficient operator before combining variable conductivity with Neumann
current closure.

Observed max-error convergence:

| n | h | max error | observed order |
|---:|---:|---:|---:|
| 17 | `6.250e-2` | `3.216e-3` | |
| 25 | `4.167e-2` | `1.428e-3` | `2.003` |
| 33 | `3.125e-2` | `8.029e-4` | `2.001` |
| 49 | `2.083e-2` | `3.568e-4` | `2.001` |
| 65 | `1.562e-2` | `2.007e-4` | `2.000` |

## Two-region conductivity-jump reference

The diagnostic `cases/lm_mhd_coupling/inductionless_two_region_sigma_mms.py`
checks a face-aligned material interface. The source is zero, the exact
potential is piecewise linear, and the normal current is continuous across the
interface:

```text
sigma_left = 1
sigma_right = 4
target interface flux = 1.6
```

The variable-conductivity reference uses harmonic face averages, which recover
the interface flux for this face-aligned two-region problem.

Default result:

| quantity | value |
|---|---:|
| max potential error | `6.217e-14` |
| max interface-flux error | `1.061e-13` |

## Smooth variable-conductivity Neumann-flux reference

The diagnostic
`cases/lm_mhd_coupling/inductionless_variable_sigma_neumann_mms.py` verifies
the variable-coefficient Neumann flux solve:

```text
div(sigma grad(phi)) = source
normal flux = sigma d(phi)/dn
```

with smooth `sigma = 1 + 0.3 x + 0.2 z` and a quadratic manufactured
potential. The wall flux is nonzero on every side and the mean value fixes the
Neumann gauge.

Observed max-error convergence:

| n | h | max error | observed order |
|---:|---:|---:|---:|
| 17 | `6.250e-2` | `4.470e-4` | |
| 25 | `4.167e-2` | `1.989e-4` | `1.997` |
| 33 | `3.125e-2` | `1.120e-4` | `1.998` |
| 49 | `2.083e-2` | `4.980e-5` | `1.998` |
| 65 | `1.562e-2` | `2.802e-5` | `1.999` |

## Full variable-conductivity inductionless MMS

The diagnostic
`cases/lm_mhd_coupling/inductionless_variable_sigma_full_mms.py` connects the
variable-coefficient pieces into the full smooth scalar-conductivity
inductionless reference path:

```text
source = div(sigma u x B)
wall flux = sigma (u x B) . n
J = sigma (-grad(phi) + u x B)
```

The manufactured construction sets
`u x B = grad(phi_exact) + J_target/sigma`, with `J_target` generated from a
streamfunction so that it is divergence-free and has zero wall-normal current.
This checks source generation, insulating wall-flux generation, the
variable-coefficient Neumann solve, current reconstruction, `div J`, and wall
current diagnostics in one loop.

Default result:

| quantity | value |
|---|---:|
| max `|phi - phi_exact|` | `3.587e-7` |
| max `|J - J_target|` | `5.390e-7` |
| max interior `|div J|` | `2.260e-7` |
| max wall-normal current | `5.390e-7` |

This is the strongest Python-side smooth-conductivity reference target so far
for a future inductionless electric-potential implementation. It is still a
manufactured finite-difference diagnostic, not a MUG solve and not a coupled
plasma/liquid-metal simulation.

## Two-region full inductionless MMS

The diagnostic `cases/lm_mhd_coupling/inductionless_two_region_full_mms.py`
combines the two-region conductivity jump with the full inductionless
manufactured path:

```text
source = div(sigma u x B)
wall flux = sigma (u x B) . n
J = sigma (-grad(phi) + u x B)
```

The exact potential is piecewise linear with continuous
`sigma d(phi)/dx` across a face-aligned vertical material interface, and a
streamfunction supplies a nonzero divergence-free target current that crosses
the interface. The diagnostic checks bulk pointwise current, wall current,
bulk `div J`, and a conservative face-current reconstruction at the material
interface.

Default result:

| quantity | value |
|---|---:|
| max `|phi - phi_exact|` | `2.759e-14` |
| max bulk `|J - J_target|` | `1.255e-13` |
| max bulk `|div J|` | `1.858e-12` |
| max wall-normal current | `3.264e-14` |
| max conservative interface-face current error | `1.228e-13` |
| max all-point pointwise current error | `5.250e-1` |

The last row is intentionally reported. It is the expected artifact from using
nodal central differences at a discontinuous conductivity interface. The
interface validation claim is therefore the conservative face-current row, not
the all-point pointwise current row.

## Offline CSV postprocessor

The utility `cases/lm_mhd_coupling/inductionless_csv_postprocess.py` provides a
structured-grid bridge for future MUG/TokaMaker field extraction. It reads
`x,z,u_x,u_z,sigma,phi` plus optional magnetic-field and target-current columns,
then computes reconstructed current, `div J`, wall-normal current, and target
current error.

The default demo writes a manufactured input CSV and recovers the target
current to roundoff:

| quantity | value |
|---|---:|
| max interior `|div J|` | `5.551e-16` |
| max wall-normal current | `8.706e-18` |
| max `|J - J_target|` | `6.939e-18` |

This is not a finite-element projection and should not be used for
conservative interface-current claims. It is an offline extraction and
sanity-check bridge.

The design note `docs/dev/lm_mhd_inductionless_weak_form_design.md` records the
weak-form signs, gauge constraint, axisymmetric metric, material-interface
caveats, and a MUG field-extraction checklist.

## Q1 weak-form prototype

The standalone script
`cases/lm_mhd_coupling/inductionless_fe_potential_mms.py` assembles the scalar
finite-element weak form

```text
int sigma grad(w).grad(phi) = int sigma grad(w).(u x B)
```

with all-insulating natural current closure and a nodal mean-value gauge. This
is a Python sign/gauge prototype only; it does not add a MUG field or Fortran
assembly path.

Observed max-error convergence:

| n | max error | observed order |
|---:|---:|---:|
| 9 | `9.333e-6` | |
| 17 | `2.338e-6` | `1.997` |
| 25 | `1.040e-6` | `1.999` |
| 33 | `5.849e-7` | `1.999` |
| 49 | `2.600e-7` | `1.999` |

The finest-grid nodal postprocessing gives max interior `|div J| = 7.546e-13`
and max wall-normal current `3.617e-7`. The recovered boundary-face
`|J.n| = 5.827e-4` on the finest grid and decreases approximately first order
over the sweep, as expected for simple Q1 element-gradient flux recovery.

## Q1 two-region interface-flux prototype

The script `cases/lm_mhd_coupling/inductionless_fe_two_region_flux_mms.py`
aligns a scalar conductivity jump with Q1 element faces, solves the same
insulating weak form, and recovers normal current on both sides of the material
interface. The exact potential is piecewise linear, so it is in the conforming
Q1 space.

Default result:

| n | max potential error | max interface target error | max interface jump |
|---:|---:|---:|---:|
| 17 | `5.551e-16` | `4.270e-15` | `3.109e-15` |
| 33 | `2.554e-15` | `1.271e-14` | `7.661e-15` |
| 49 | `1.238e-14` | `6.240e-14` | `8.882e-15` |

This is still a standalone Python prototype, but it is the strongest
material-interface sign check in the reference stack because it avoids nodal
central differences at the conductivity jump.

## Axisymmetric Q1 weak-form prototype

The script
`cases/lm_mhd_coupling/inductionless_axisymmetric_fe_potential_mms.py`
assembles the scalar weak form with the `R` metric:

```text
int R sigma grad(w).grad(phi) dR dZ
  = int R sigma grad(w).(u x B) dR dZ
```

The manufactured current is axisymmetrically divergence-free:

```text
(1/R) d(R J_R)/dR + dJ_Z/dZ = 0
```

Observed max-error convergence:

| n | max error | observed order |
|---:|---:|---:|
| 9 | `1.943e-5` | |
| 17 | `4.886e-6` | `1.992` |
| 25 | `2.175e-6` | `1.995` |
| 33 | `1.225e-6` | `1.997` |
| 49 | `5.448e-7` | `1.998` |

Finest-grid nodal postprocessing gives max axisymmetric `|div J| = 3.123e-7`
and max wall-normal current `1.356e-6`. These are structured-grid
postprocessing diagnostics, not conservative FE flux recovery.

## MUG restart export-readiness manifest

The utility `cases/lm_mhd_coupling/mug_restart_manifest.py` inspects the
available resolved-Hartmann MUG HDF5 restart before attempting any offline
inductionless exporter. On the default Ha=2 restart it finds:

| item | measured value |
|---|---:|
| restart time | `0` |
| restart state length | `644` |
| profile sidecar rows | `725` |
| restart/profile direct index match | `False` |
| coordinates/connectivity in restart | `False` |
| electric potential in restart | `False` |

This means a real MUG-to-CSV exporter cannot honestly be built from this
restart alone. It needs mesh/plot coordinates plus state fields, and current
`xmhd_2d` does not contain the inductionless electric potential field required
by `inductionless_csv_postprocess.py`.

## Axisymmetric toroidal-feedback diagnostic

The script
`cases/lm_mhd_coupling/axisymmetric_inductionless_toroidal_feedback.py`
connects the axisymmetric inductionless current reconstruction to the existing
`CouplingSnapshot` magnetic-feedback machinery. It reconstructs

```text
J = sigma(-grad(phi) + u x B + E_phi e_phi)
```

and then estimates poloidal magnetic-field perturbations from `J_phi` using
circular-loop superposition. The selected default case has `E_phi = 0.02 V/m`,
representative `sigma0 = 8e5 S/m`, and `B_phi = 5 T`.

| quantity | value |
|---|---:|
| total toroidal current | `1.599e4 A` |
| max `|delta B|/|B_ref|` | `1.656e-3` |
| feedback class | `two-way-candidate` |
| max `|J_phi|` | `1.696e4 A/m^2` |
| max `|J x B|` | `3.917e3 N/m^3` |
| max Joule heating | `3.524e2 W/m^3` |
| max interior axisymmetric `|div J|` | `1.858e-10` |
| max wall-normal current | `0.000e0` |

This is the most directly two-way-coupling-relevant diagnostic in the Python
stack so far, but it is still manufactured. The potential and toroidal electric
field are prescribed; no Grad-Shafranov response is solved.

The companion script `cases/lm_mhd_coupling/axisymmetric_loop_voltage_feedback.py`
sets the toroidal electric field from loop voltage:

```text
E_phi(R) = V_loop / (2*pi*R)
```

For `V_loop = 0.2 V`, the selected default result is:

| quantity | value |
|---|---:|
| `E_phi` range | `1.592e-2` to `3.183e-2 V/m` |
| total toroidal current | `1.756e4 A` |
| max `|delta B|/|B_ref|` | `1.909e-3` |
| feedback class | `two-way-candidate` |
| max `|J_phi|` | `2.469e4 A/m^2` |
| max `|J x B|` | `4.977e3 N/m^3` |
| max Joule heating | `7.972e2 W/m^3` |
| max interior axisymmetric `|div J|` | `1.858e-10` |
| max wall-normal current | `0.000e0` |

The threshold scan
`cases/lm_mhd_coupling/axisymmetric_loop_voltage_threshold_scan.py` estimates
the loop voltage needed to cross the planning thresholds as conductivity
varies. For the PbLi-like row nearest `sigma0 = 8e5 S/m`, it finds:

| threshold | estimated loop voltage |
|---|---:|
| `|delta B|/B_ref = 1e-3` | `1.046e-1 V` |
| `|delta B|/B_ref = 1e-2` | not reached up to `1 V` |
| max feedback at `1 V` | `9.534e-3` |

## Why this matters

The roadmap correctly identifies inductionless electric-potential MHD as the
main computational unlock for blanket flows. It is also invasive: it needs a
new potential field, a charge-conservation solve, boundary-current closure, and
careful coupling of `J x B` into momentum. This slice does not implement that
solver. It creates a low-risk reference layer that future solver work can use
to check signs, shape conventions, units, divergence diagnostics, wall-normal
current diagnostics, Lorentz force, and Joule dissipation.

## Red-team limits

- A scalar electric-potential Poisson equation is solved only in uniform-grid
  finite-difference reference problems: Dirichlet, homogeneous Neumann,
  prescribed Neumann wall data, constant-conductivity `u x B` MMS, and smooth
  scalar-conductivity `sigma u x B` MMS.
- No finite-element weak form is implemented.
- The full inductionless MMS cases have no conductivity tensors, no
  conducting-wall Robin condition, and no resolved Hartmann/Shercliff/Hunt
  velocity field.
- The two-region conductivity-jump MMS is face-aligned and one-dimensional in
  x; the full two-region reconstruction handles this aligned case only.
- Pointwise central-difference current reconstruction is not treated as valid
  at discontinuous material interfaces; interface-current claims use a
  conservative face-current calculation.
- The CSV postprocessor assumes a structured grid and is not a replacement for
  a finite-element weak-form residual or conservative interface flux.
- The Q1 finite-element prototype is standalone Python. It checks weak-form
  sign, gauge, and convergence. Boundary-face flux is recovered from element
  gradients; the two-region Q1 prototype adds conforming material-interface
  flux recovery for a face-aligned scalar conductivity jump.
- The axisymmetric Q1 prototype includes the scalar `R` metric but omits
  toroidal electric field, tensor conductivity, and live TokaMaker feedback.
- The MUG restart manifest is an export-readiness check only. The inspected
  restart is not the final developed-flow profile and lacks mesh coordinates.
- The axisymmetric toroidal-feedback diagnostic estimates magnetic feedback
  from manufactured `J_phi`; it is not a TokaMaker equilibrium update.
- The loop-voltage continuation derives `E_phi(R)` from `V_loop`, but still
  prescribes velocity and potential fields.
- The threshold scan is a coarse manufactured sensitivity map, not a design
  envelope.
- The MMS current is manufactured and is not a resolved Hartmann, Shercliff, or
  Hunt flow.
- The structured-grid divergence diagnostic is for reference/postprocessing; it
  is not the operator that MUG or TokaMaker would use.

## Next step

The next non-invasive reference step should put the loop-voltage-driven
`E_phi(R)` into the standalone axisymmetric weak-form prototype, or
locate/generate MUG plot output with coordinates plus final-time state fields
for a conservative CSV exporter.
