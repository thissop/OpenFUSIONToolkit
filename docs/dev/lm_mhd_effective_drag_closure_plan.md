# LM-MHD effective drag closure diagnostics

Prepared on branch `lm-mhd-upgrades`, 2026-07-01.

## Status and scope

This is an analytic/reduced-model infrastructure pass, not a MUG validation run.
No Fortran solver physics was changed. The goal was to make the unresolved
Hartmann-layer question executable and auditable for future full-device
plasma/blanket work.

The motivation follows Sophia's clarification: a standalone channel can resolve
Hartmann layers and does not need a reduced drag term for the resolved
benchmark. A full-device coupled plasma/blanket mesh generally cannot afford to
resolve every Hartmann layer, so it needs a calibrated effective drag or
pressure-drop closure in underresolved LM regions.

## Added capability

The new `OpenFUSIONToolkit.LM_MHD.closures` module provides:

- `hartmann_layer_thickness(a, Ha) = a / Ha`.
- `cells_per_hartmann_layer(dx, a, Ha)`.
- `classify_hartmann_resolution(...)` with labels `resolved`, `marginal`, and
  `underresolved`.
- `hartmann_channel_effective_drag(a, Ha, nu)`, the mean-flow equivalent
  damping rate for an insulating Hartmann channel driven by a uniform body force.
- Drag projection helpers: `perpendicular_velocity`, `parallel_velocity`,
  `implicit_drag_update`, and `drag_power_density`.

The effective drag is defined by matching the analytic mean velocity, not the
pointwise velocity profile:

```text
u_mean = f a^2 / nu * [1 - tanh(Ha)/Ha] / Ha^2
alpha_eff = f / u_mean
          = (nu/a^2) Ha^2 / [1 - tanh(Ha)/Ha]
```

Its zero-field limit is the Poiseuille mean-flow damping rate `3 nu/a^2`.
At large Ha it approaches `(nu/a^2) Ha^2`. This is a useful calibration target
for a reduced balance, but it must not be presented as a pointwise
Hartmann-layer model.

## Generated artifacts

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 cases/lm_mhd_coupling/hartmann_drag_closure_scan.py
```

Outputs under `cases/lm_mhd_coupling/results/`:

- `hartmann_drag_closure_scan.csv`
- `hartmann_layer_resolution_map.png`
- `hartmann_effective_drag_alpha.png`
- `drag_closure_validity_matrix.png`
- `manufactured_current_feedback.h5`
- `manufactured_current_feedback.png`
- `hartmann_drag_closure_summary.md`

The default scan uses `a = 0.05 m`, `nu = 1.8e-7 m^2/s`, and a resolved
threshold of eight cells per Hartmann layer. In that grid, the generated summary
reports:

| label | count |
|---|---:|
| resolved | 5202 |
| marginal | 1270 |
| underresolved | 5288 |

The manufactured current-feedback diagnostic writes and reads a
`CouplingSnapshot` HDF5 file, then plots synthetic `J_phi`, a sheet-current
field estimate, and a feedback flag. The default manufactured case reports
`max |delta B|/|B_ref| = 1.567e-02` and is classified as
`two-way-candidate`. This is a schema and triage check only, not a coupled
TokaMaker result.

## Verification performed

Focused tests:

```text
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src/python pytest -q \
  src/tests/physics/test_lm_mhd_closures.py \
  src/tests/physics/test_lm_mhd_coupling.py
```

Result:

```text
12 passed
```

Additional checks performed before commit:

- Regenerated all new closure artifacts with the script above.
- Used PIL to confirm each generated PNG has nonzero dimensions and nonblank
  bounding boxes.
- Ran `git diff --check`.

## Red-team constraints for future use

- Do not use `alpha_eff` to claim pointwise Hartmann profile agreement.
- Do not use this closure in resolved standalone channel verification; resolved
  layers are the reference.
- Do not claim the closure supplies wall-current closure. In the coupled
  plasma/blanket formulation, current closure should come through the
  TokaMaker/finite-element field continuity path.
- Do not use rounded default material properties for design claims.
- Before wiring this into Fortran solver physics, calibrate against resolved
  Hartmann/Shercliff references and verify energy dissipation.

## Recommended next implementation step

Add a default-off calibrated drag model interface on top of the existing
`use_wall_drag` machinery only after the calibration artifacts are reviewed.
The minimum safe Fortran-facing extension would be a typed model selector plus
metadata for `(a, Ha, alpha_eff)` so runs can record whether drag is constant,
calibrated, or disabled. Existing tests must continue to pass with the drag path
off.
