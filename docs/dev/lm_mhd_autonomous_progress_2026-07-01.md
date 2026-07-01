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

No push was performed.

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
- `docs/dev/lm_mhd_inductionless_reference_plan.md`

Interpretation: this is the safest on-ramp toward the roadmap's highest-payoff
feature, an inductionless electric-potential mode. It establishes signs,
operators, gauges, compatibility checks, and MMS artifacts before any Fortran
field-structure work.

## Verification

Latest focused test command:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src/python pytest -q \
  src/tests/physics/test_lm_mhd_coupling.py \
  src/tests/physics/test_lm_mhd_closures.py \
  src/tests/physics/test_lm_mhd_feedback.py \
  src/tests/physics/test_lm_mhd_inductionless.py
```

Latest result: `30 passed`.

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
