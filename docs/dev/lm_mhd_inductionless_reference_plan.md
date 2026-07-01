# LM-MHD inductionless reference infrastructure

This note records the third autonomous LM-MHD upgrade slice from 2026-07-01.
It addresses the roadmap's highest-impact long-term upgrade, the low magnetic
Reynolds number electric-potential mode, but only at the Python reference and
diagnostic level. No Fortran solver physics was changed.

## What changed

Added `OpenFUSIONToolkit.LM_MHD.inductionless` with small tested utilities:

- `ohms_law_current_density(sigma, velocity, magnetic_field, electric_potential_gradient)`
- `lorentz_force_density(current_density, magnetic_field)`
- `joule_heating_density(current_density, sigma)`
- `structured_gradient_2d(scalar, x, z)`
- `motional_electric_field(velocity, magnetic_field)`
- `potential_source_from_motional_emf(motional_emf, x, z)`
- `insulating_wall_normal_gradient(motional_emf, x, z)`
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
  prescribed Neumann wall data, and constant-conductivity `u x B` MMS.
- No finite-element weak form is implemented.
- The constant-conductivity `u x B` MMS has no material jumps, conductivity
  tensors, conducting-wall Robin condition, or current continuity across
  regions.
- The MMS current is manufactured and is not a resolved Hartmann, Shercliff, or
  Hunt flow.
- The structured-grid divergence diagnostic is for reference/postprocessing; it
  is not the operator that MUG or TokaMaker would use.

## Next step

The next non-invasive reference step should be a variable-conductivity MMS for
`div(sigma grad phi) = div(sigma u x B)`, including either smooth
`sigma(x,z)` or a two-region jump with current-continuity checks. That would
exercise the part of the eventual finite-element mode most likely to matter for
multi-material blankets while still avoiding risky Fortran edits.
