# LM-MHD autonomous progress ledger, 2026-07-01

This ledger summarizes the autonomous LM-MHD work after the effective-drag
closure request. All new results here are analytic, reduced-model, or Python
reference diagnostics. They are not MUG validation, not a coupled TokaMaker/MUG
solve, and not a claim of new liquid-metal solver capability.

## Local commits from this run

| commit | subject |
|---|---|
| `3e79f32` | `lm-mhd: add effective drag closure diagnostics` |
| `00cbb5b` | `lm-mhd: add axisymmetric feedback diagnostics` |
| `2c82a1f` | `lm-mhd: add inductionless reference diagnostics` |
| `fabe450` | `lm-mhd: add inductionless poisson reference mms` |
| `98ea81d` | `lm-mhd: add neumann potential reference mms` |
| `27594a4` | `lm-mhd: add nonzero neumann potential reference` |
| `9d056db` | `lm-mhd: document autonomous diagnostic progress` |
| `4fa5ad5` | `lm-mhd: add full inductionless potential mms` |
| `ad2f2c8` | `lm-mhd: add variable conductivity potential mms` |
| `bf3c722` | `lm-mhd: add two-region conductivity reference mms` |
| `357fcd7` | `lm-mhd: add variable conductivity neumann reference` |

No push was performed. The continuation after `357fcd7` is intentionally left
uncommitted because the latest overnight instruction says to leave changes
uncommitted for review.

## What was added

### Effective drag closure infrastructure

Added `OpenFUSIONToolkit.LM_MHD.closures` with Hartmann-layer resolution,
mean-response effective drag, projection, implicit drag-update, and drag-power
diagnostics. The saved scan artifacts are under
`cases/lm_mhd_coupling/results/`, including:

- `hartmann_drag_closure_scan.csv`
- `hartmann_layer_resolution_map.png`
- `hartmann_effective_drag_alpha.png`
- `drag_closure_validity_matrix.png`
- `hartmann_drag_closure_summary.md`

Interpretation: reduced drag is for underresolved full-device meshes, not for
resolved standalone Hartmann verification.

### Two-way feedback triage

Added axisymmetric circular-loop Biot-Savart feedback diagnostics using the
existing `CouplingSnapshot` HDF5 schema. The baseline manufactured current patch
with peak `J_phi = 1e5 A/m^2` gives max probe
`|delta B|/|B_ref| = 2.576e-3`, labeled `two-way-candidate`.

Key artifacts:

- `axisymmetric_feedback_scan.csv`
- `axisymmetric_feedback_snapshot.h5`
- `axisymmetric_feedback_field_map.png`
- `axisymmetric_feedback_scaling.png`
- `axisymmetric_feedback_summary.md`
- `docs/dev/lm_mhd_two_way_feedback_diagnostics.md`

Interpretation: two-way coupling is strategically important, but not always
dominant. A feedback ratio below `1e-3` supports one-way imposed-field scans as
a first pass; values near or above `1e-3` deserve coupled sensitivity studies;
values near `1e-2` are strong prompts for reverse coupling.

### Inductionless reference infrastructure

Added `OpenFUSIONToolkit.LM_MHD.inductionless` with tested reference utilities
for:

- `J = sigma (-grad(phi) + u x B)`
- `f = J x B`
- `q = |J|^2 / sigma`
- structured `div J` diagnostics
- streamfunction-generated divergence-free currents
- uniform-grid Dirichlet Poisson reference solves
- uniform-grid Neumann Poisson reference solves, including nonzero wall-normal
  derivative data
- variable-conductivity Dirichlet, Neumann, interface, and full
  `sigma (u x B)` manufactured references

Key artifacts:

- `inductionless_mms_metrics.csv`
- `inductionless_mms_fields.png`
- `inductionless_mms_midplane.png`
- `inductionless_poisson_mms_metrics.csv`
- `inductionless_poisson_mms_solution.png`
- `inductionless_poisson_mms_convergence.png`
- `inductionless_neumann_mms_metrics.csv`
- `inductionless_neumann_mms_solution.png`
- `inductionless_neumann_mms_convergence.png`
- `inductionless_neumann_wall_mms_metrics.csv`
- `inductionless_neumann_wall_mms_solution.png`
- `inductionless_full_mms_metrics.csv`
- `inductionless_full_mms_potential.png`
- `inductionless_full_mms_current.png`
- `inductionless_variable_sigma_mms_metrics.csv`
- `inductionless_variable_sigma_mms_solution.png`
- `inductionless_variable_sigma_mms_convergence.png`
- `inductionless_two_region_sigma_mms_metrics.csv`
- `inductionless_two_region_sigma_mms_solution.png`
- `inductionless_two_region_sigma_mms_flux.png`
- `inductionless_variable_sigma_neumann_mms_metrics.csv`
- `inductionless_variable_sigma_neumann_mms_solution.png`
- `inductionless_variable_sigma_neumann_mms_convergence.png`
- `inductionless_variable_sigma_full_mms_metrics.csv`
- `inductionless_variable_sigma_full_mms_potential.png`
- `inductionless_variable_sigma_full_mms_current.png`
- `inductionless_two_region_full_mms_metrics.csv`
- `inductionless_two_region_full_mms_potential.png`
- `inductionless_two_region_full_mms_current.png`
- `inductionless_two_region_full_mms_interface.png`
- `inductionless_postprocess_demo_input.csv`
- `inductionless_postprocess_demo_metrics.csv`
- `inductionless_postprocess_demo_current.png`
- `inductionless_fe_potential_mms_metrics.csv`
- `inductionless_fe_potential_mms_solution.png`
- `inductionless_fe_potential_mms_convergence.png`
- `inductionless_fe_two_region_flux_mms_metrics.csv`
- `inductionless_fe_two_region_flux_mms_solution.png`
- `inductionless_fe_two_region_flux_mms_interface.png`
- `inductionless_axisymmetric_fe_potential_mms_metrics.csv`
- `inductionless_axisymmetric_fe_potential_mms_solution.png`
- `inductionless_axisymmetric_fe_potential_mms_convergence.png`
- `mug_restart_manifest.csv`
- `mug_restart_field_ranges.png`
- `axisymmetric_inductionless_toroidal_feedback_scan.csv`
- `axisymmetric_inductionless_toroidal_feedback_maps.png`
- `axisymmetric_inductionless_toroidal_feedback_scan.png`
- `axisymmetric_inductionless_toroidal_feedback_snapshot.h5`
- `axisymmetric_loop_voltage_feedback_scan.csv`
- `axisymmetric_loop_voltage_feedback_maps.png`
- `axisymmetric_loop_voltage_feedback_scan.png`
- `axisymmetric_loop_voltage_feedback_snapshot.h5`
- `axisymmetric_loop_voltage_threshold_scan.csv`
- `axisymmetric_loop_voltage_thresholds.csv`
- `axisymmetric_loop_voltage_threshold_heatmap.png`
- `axisymmetric_loop_voltage_thresholds.png`
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
- `docs/dev/lm_mhd_inductionless_reference_plan.md`
- `docs/dev/lm_mhd_inductionless_weak_form_design.md`

Interpretation: this is the safest on-ramp toward the roadmap's highest-payoff
feature, an inductionless electric-potential mode. It establishes signs,
operators, gauges, compatibility checks, and MMS artifacts before any Fortran
field-structure work.

Latest uncommitted continuation: added global axisymmetric integral utilities
for `f.u` mechanical power density, `R f_phi` torque density, volume
integration, integrated power, and integrated torque. The selected
constant-`E_phi` manufactured case reports integrated Joule power
`3.022e3 W`, integrated mechanical power `2.969e0 W`, and net torque
`9.311e0 N m`; the selected loop-voltage case reports integrated Joule power
`3.513e3 W`, integrated mechanical power `-4.551e-1 W`, and net torque
`9.311e0 N m`. These numbers are reduced diagnostic audits of the manufactured
fields only.

The same continuation added signed outward axisymmetric boundary-current
integrals. The regenerated selected constant-`E_phi` and loop-voltage cases
both report max integrated boundary current `0.000e0 A` and net integrated
boundary current `0.000e0 A`, complementing the differential `div J` residual
and pointwise wall-normal current extrema.

The latest handoff-facing continuation added
`cases/lm_mhd_coupling/coupling_snapshot_audit.py` and the tested helper
`structured_snapshot_current_grid(snapshot)`. The default audit reads the
loop-voltage `CouplingSnapshot`, detects an `81 x 65` structured grid, recovers
max `|delta B|/|B_ref| = 1.909e-3`, classifies the snapshot as
`two-way-candidate`, and reports zero integrated boundary leakage plus zero
cell-volume relative error against `2*pi*R*dR*dZ`.

The latest red-team continuation added
`cases/lm_mhd_coupling/axisymmetric_feedback_probe_sensitivity.py`. For the
same loop-voltage snapshot, 50% of probe/regularization rows cross the `1e-3`
two-way-candidate threshold, no rows cross `1e-2`, the strongest row is
`2.047e-3`, and the weakest row is `2.361e-4`. This makes the priority for
reverse-coupling diagnostics plausible but not universal.

The latest regime-priority continuation added
`cases/lm_mhd_coupling/lm_mhd_regime_priority_map.py`. At the representative
point nearest `U = 0.2 m/s`, `B = 5 T`, PbLi is `drag-important` with feedback
`4.966e-4`, while Li is `two-way-candidate` with feedback `1.986e-3`. Across
the scanned B-U grid, PbLi is `drag-important` in 90.2% of rows and
`two-way-candidate` in 9.8%; Li is `drag-important` in 65.5%,
`two-way-candidate` in 19.5%, and `full-induction-candidate` in 15.0%.

The latest artifact-audit continuation added
`cases/lm_mhd_coupling/lm_mhd_artifact_manifest.py`. The generated manifest
contains 32 entries: 9 CSV files, 14 PNG files, 2 HDF5 files, and 7 Markdown
summaries; all 14 PNG files are nonblank. It records SHA-256 hashes and row
counts for the generated diagnostic artifacts.

The latest model-choice continuation added
`cases/lm_mhd_coupling/full_induction_low_rm_bridge.py` and focused tests in
`src/tests/physics/test_lm_mhd_full_induction_bridge.py`. In the selected
PbLi-like 1D bridge, `Rm = 1.005e-2`, `max |b_x|/B0 = 1.290e-3`, and the
Ohm-vs-Ampere current relative error is `1.875e-5`. This does not replace full
MUG induction, but it gives a parent-model check that the inductionless current
closure is the correct low-`Rm` limit for this controlled case.

The same continuation added the pointwise inductionless power identity
`q_J + f.u - J.E = 0`. The selected constant-`E_phi` and loop-voltage cases
close the integrated residual to `~1e-14 W`.

## Verification

Latest focused test command:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src/python pytest -q \
  src/tests/physics/test_lm_mhd_coupling.py \
  src/tests/physics/test_lm_mhd_closures.py \
  src/tests/physics/test_lm_mhd_feedback.py \
  src/tests/physics/test_lm_mhd_inductionless.py \
  src/tests/physics/test_lm_mhd_full_induction_bridge.py
```

Latest result before the global-integral continuation: `38 passed`.
Inductionless-only result after adding global-integral utilities: `26 passed`.
Inductionless-only result after adding boundary-current integrals: `27 passed`.
Full focused LM-MHD Python result after the global-integral continuation:
`45 passed`.
Full focused LM-MHD Python result after the boundary-current continuation:
`46 passed`.
Full focused LM-MHD Python result after the `CouplingSnapshot` audit
continuation: `48 passed`.
Full focused LM-MHD Python result after the probe-sensitivity continuation:
`48 passed`.
Full focused LM-MHD Python result after the regime-priority continuation:
`48 passed`.
Full focused LM-MHD Python result after audit cell-volume hardening:
`48 passed`.
Full focused LM-MHD Python result after artifact-manifest continuation:
`48 passed`.
Expanded focused LM-MHD Python result after the power-identity and
full-induction bridge continuation: `53 passed`.

Additional checks performed repeatedly:

- `git diff --check`
- Re-ran each diagnostic script after edits.
- PIL image-size/nonblank checks for new PNGs.
- Visual inspection of the principal saved plots.

## Dirty worktree boundary

The following preexisting dirty items were intentionally left out of these
commits:

- modified `cases/hartmann/results/*` and `cases/hartmann/results/report/*`
- modified `src/physics/xmhd_2d.F90`
- untracked Hartmann packaging/deliverable files and local notes

Those are separate from the committed LM-MHD Python diagnostics above.

## Best next technical step

The next non-invasive reference step originally identified here was to connect
the pieces:

1. Choose an MMS velocity field and uniform applied field.
2. Generate the scalar source `div(sigma u x B)`.
3. Generate insulating wall data `d phi / dn = (u x B) . n`.
4. Solve `div(sigma grad phi) = div(sigma u x B)` with the nonzero Neumann
   reference.
5. Reconstruct `J = sigma(-grad phi + u x B)` and verify `div J = 0` plus
   zero wall-normal current.

That would be the most rigorous Python-side target before touching MUG's
Fortran field structure for an inductionless mode.

## Continuation update

This step was implemented in the later autonomous continuation as the
constant-conductivity diagnostic
`cases/lm_mhd_coupling/inductionless_full_mms.py`. It generates
`u x B = grad(phi_exact) + J_target/sigma`, solves the Neumann potential
problem with wall data `(u x B).n`, reconstructs `J`, and verifies charge
conservation plus zero wall-normal current.

Default result:

| quantity | value |
|---|---:|
| max `|phi - phi_exact|` | `1.102e-16` |
| max `|J - J_target|` | `3.553e-10 A/m^2` |
| max interior `|div J|` | `2.235e-8 A/m^3` |
| max wall-normal current | `3.553e-10 A/m^2` |

The next best non-invasive target is now variable-conductivity inductionless
MMS: `div(sigma grad phi) = div(sigma u x B)` with smooth or two-region
`sigma(x,z)`, current-continuity checks, and no Fortran edits until that
reference is solid.

## Additional Continuation Update

The smooth variable-conductivity reference was implemented as
`cases/lm_mhd_coupling/inductionless_variable_sigma_mms.py`. It solves
`div(sigma grad(phi)) = source` with
`sigma = 1 + 0.3 x + 0.2 z`, Dirichlet walls, and an analytic MMS source.

Default convergence:

| n | max error | observed order |
|---:|---:|---:|
| 17 | `3.216e-3` | |
| 25 | `1.428e-3` | `2.003` |
| 33 | `8.029e-4` | `2.001` |
| 49 | `3.568e-4` | `2.001` |
| 65 | `2.007e-4` | `2.000` |

The two-region conductivity-jump reference was then implemented as
`cases/lm_mhd_coupling/inductionless_two_region_sigma_mms.py`. It uses
harmonic face averaging across a face-aligned interface with
`sigma_left = 1`, `sigma_right = 4`, zero source, and piecewise-linear exact
potential.

Default result:

| quantity | value |
|---|---:|
| target interface flux | `1.600e0` |
| max potential error | `6.217e-14` |
| max interface-flux error | `1.061e-13` |

The next best target is now variable-coefficient Neumann current closure:
`div(sigma grad phi) = div(sigma u x B)` with insulating `J.n = 0` wall checks
for spatially varying or regional sigma.

## Variable-Coefficient Neumann Update

The smooth variable-conductivity Neumann-flux reference was implemented as
`cases/lm_mhd_coupling/inductionless_variable_sigma_neumann_mms.py`. It solves
`div(sigma grad(phi)) = source` with nonzero prescribed outward wall flux
`sigma dphi/dn` and a mean-value gauge.

Default convergence:

| n | max error | observed order |
|---:|---:|---:|
| 17 | `4.470e-4` | |
| 25 | `1.989e-4` | `1.997` |
| 33 | `1.120e-4` | `1.998` |
| 49 | `4.980e-5` | `1.998` |
| 65 | `2.802e-5` | `1.999` |

This target was completed in the latest uncommitted continuation as
`cases/lm_mhd_coupling/inductionless_variable_sigma_full_mms.py`. It generates
source and wall flux from `sigma (u x B)`, solves the variable-coefficient
Neumann potential problem, reconstructs `J = sigma(-grad phi + u x B)`, and
checks `div J` plus insulating wall current.

Default result:

| quantity | value |
|---|---:|
| max `|phi - phi_exact|` | `3.587e-7` |
| max `|J - J_target|` | `5.390e-7` |
| max interior `|div J|` | `2.260e-7` |
| max wall-normal current | `5.390e-7` |

Artifacts:

- `inductionless_variable_sigma_full_mms_metrics.csv`
- `inductionless_variable_sigma_full_mms_potential.png`
- `inductionless_variable_sigma_full_mms_current.png`
- `inductionless_variable_sigma_full_mms_summary.md`

That next target was completed as
`cases/lm_mhd_coupling/inductionless_two_region_full_mms.py`. It combines a
face-aligned conductivity jump, `u x B` source generation, insulating wall
flux, current reconstruction, conservative interface normal-current
continuity, and `div J` checks.

Default result:

| quantity | value |
|---|---:|
| max `|phi - phi_exact|` | `2.759e-14` |
| max bulk `|J - J_target|` | `1.255e-13` |
| max bulk `|div J|` | `1.858e-12` |
| max wall-normal current | `3.264e-14` |
| max conservative interface-face current error | `1.228e-13` |
| max all-point pointwise current error | `5.250e-1` |

The final row is deliberately retained as a warning: nodal central differences
are not a conservative current operator at a discontinuous material interface.
The interface claim is based on the conservative face current.

That target was advanced with
`docs/dev/lm_mhd_inductionless_weak_form_design.md` and
`cases/lm_mhd_coupling/inductionless_csv_postprocess.py`. The note records the
electric-potential weak form, gauge constraint, axisymmetric metric,
material-interface caveats, and MUG extraction conventions. The postprocessor
reads structured `x,z,u_x,u_z,sigma,phi` CSV data and computes reconstructed
current, `div J`, wall-normal current, and optional target-current errors.

Default postprocessor demo:

| quantity | value |
|---|---:|
| max interior `|div J|` | `5.551e-16` |
| max wall-normal current | `8.706e-18` |
| max `|J - J_target|` | `6.939e-18` |

The standalone finite-element scalar-potential matrix prototype was then added
as `cases/lm_mhd_coupling/inductionless_fe_potential_mms.py`. It assembles the
Q1 weak form with all-insulating natural current closure and a nodal
mean-value gauge, and now also recovers boundary-face `J.n` from Q1 element
gradients.

Default convergence:

| n | max error | observed order |
|---:|---:|---:|
| 9 | `9.333e-6` | |
| 17 | `2.338e-6` | `1.997` |
| 25 | `1.040e-6` | `1.999` |
| 33 | `5.849e-7` | `1.999` |
| 49 | `2.600e-7` | `1.999` |

Finest-grid recovered boundary-face `|J.n| = 5.827e-4`, decreasing
approximately first order over the sweep. This is a flux-recovery diagnostic,
not a production FE residual.

Material-interface FE flux recovery was then added as
`cases/lm_mhd_coupling/inductionless_fe_two_region_flux_mms.py`. It aligns the
conductivity jump to Q1 element faces and recovers normal current from both
sides of the interface.

Default result:

| n | max potential error | max interface target error | max interface jump |
|---:|---:|---:|---:|
| 17 | `5.551e-16` | `4.270e-15` | `3.109e-15` |
| 33 | `2.554e-15` | `1.271e-14` | `7.661e-15` |
| 49 | `1.238e-14` | `6.240e-14` | `8.882e-15` |

The axisymmetric `R`-metric prototype was then added as
`cases/lm_mhd_coupling/inductionless_axisymmetric_fe_potential_mms.py`. It
uses an axisymmetrically divergence-free manufactured current and the weighted
weak form `int R sigma grad(w).grad(phi) = int R sigma grad(w).(u x B)`.

Default convergence:

| n | max error | observed order |
|---:|---:|---:|
| 9 | `1.943e-5` | |
| 17 | `4.886e-6` | `1.992` |
| 25 | `2.175e-6` | `1.995` |
| 33 | `1.225e-6` | `1.997` |
| 49 | `5.448e-7` | `1.998` |

A measured MUG export-readiness check was added as
`cases/lm_mhd_coupling/mug_restart_manifest.py`. On the default resolved
Hartmann Ha=2 restart it reports:

| item | measured value |
|---|---:|
| restart time | `0` |
| restart state length | `644` |
| profile sidecar rows | `725` |
| direct restart/profile index match | `False` |
| coordinates/connectivity in restart | `False` |
| electric potential in restart | `False` |

The next best target is now locating or generating MUG plot output with
final-time fields plus coordinates/connectivity, then exporting only measured
fields into the structured CSV schema used by the postprocessor. Toroidal
electric-field extensions to the axisymmetric weak-form prototype are another
separate path.

## Axisymmetric Toroidal-Feedback Update

The next auto-mode slice added axisymmetric inductionless current utilities and
`cases/lm_mhd_coupling/axisymmetric_inductionless_toroidal_feedback.py`.
The reconstruction includes an explicit toroidal electric field:

```text
J = sigma(-grad(phi) + u x B + E_phi e_phi)
```

Default selected result:

| quantity | value |
|---|---:|
| `E_phi` | `2.000e-2 V/m` |
| total toroidal current | `1.599e4 A` |
| max `|delta B|/|B_ref|` | `1.656e-3` |
| feedback class | `two-way-candidate` |
| max `|J_phi|` | `1.696e4 A/m^2` |
| max `|J x B|` | `3.917e3 N/m^3` |
| max Joule heating | `3.524e2 W/m^3` |
| max interior axisymmetric `|div J|` | `1.858e-10` |
| max wall-normal current | `0.000e0` |

The loop-voltage/conductivity threshold scan
`cases/lm_mhd_coupling/axisymmetric_loop_voltage_threshold_scan.py` then
estimated the feedback-transition voltage. For the PbLi-like row nearest
`sigma0 = 8e5 S/m`, it finds:

| threshold | estimated loop voltage |
|---|---:|
| `|delta B|/B_ref = 1e-3` | `1.046e-1 V` |
| `|delta B|/B_ref = 1e-2` | not reached up to `1 V` |
| max feedback at `1 V` | `9.534e-3` |

This strengthens the earlier two-way-coupling answer: two-way feedback is not
automatically dominant, but plausible toroidal-current drives can cross the
`1e-3` planning threshold, so they deserve explicit diagnostics in any serious
coupled LM/plasma workflow.

The loop-voltage companion diagnostic
`cases/lm_mhd_coupling/axisymmetric_loop_voltage_feedback.py` then replaced the
constant toroidal electric field with `E_phi(R) = V_loop/(2*pi*R)`.

Default selected result:

| quantity | value |
|---|---:|
| `V_loop` | `2.000e-1 V` |
| `E_phi` range | `1.592e-2` to `3.183e-2 V/m` |
| total toroidal current | `1.756e4 A` |
| max `|delta B|/|B_ref|` | `1.909e-3` |
| feedback class | `two-way-candidate` |
| max `|J_phi|` | `2.469e4 A/m^2` |
| max `|J x B|` | `4.977e3 N/m^3` |
| max Joule heating | `7.972e2 W/m^3` |
| max interior axisymmetric `|div J|` | `1.858e-10` |
| max wall-normal current | `0.000e0` |
