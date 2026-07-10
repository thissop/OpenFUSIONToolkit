# LM-MHD inductionless weak-form design notes

Prepared on branch `lm-mhd-upgrades`, 2026-07-01.

## Scope

This is a design and postprocessing note for a future low magnetic Reynolds
number electric-potential mode. It is not a MUG implementation, not a
TokaMaker/MUG coupled solve, and not a validation claim. The current code added
in this pass is limited to Python reference diagnostics and an offline CSV
postprocessor.

## Strong form and signs

The Python reference utilities use:

```text
J = sigma (-grad(phi) + u x B)
f = J x B
q = |J|^2 / sigma
```

Charge conservation gives:

```text
div J = 0
div(sigma grad(phi)) = div(sigma u x B)
```

For an outward normal `n`, the normal current is:

```text
j_n = J . n = sigma (-(d phi/dn) + (u x B).n)
```

Insulating walls impose `j_n = 0`, equivalently
`sigma d(phi)/dn = sigma (u x B).n`.

## Weak form

Let `e = u x B` and let `w` be a scalar test function. Starting from
`div(sigma(-grad(phi) + e)) = 0`, the weak form with prescribed outward normal
current `j_n` is:

```text
Integral_Omega sigma grad(w).grad(phi) dV
  = Integral_Omega sigma grad(w).e dV
    - Integral_boundary w j_n dS
```

For insulating walls, `j_n = 0`, so the boundary term vanishes. This is the
cleanest finite-element target because insulating current closure is natural in
the weak form. The nullspace must still be removed with either a mean-value
constraint or a single potential gauge point.

For conducting walls, the boundary model must be specified before
implementation. A simple Robin-style closure might relate `j_n` to wall
potential or wall conductance, but the sign, units, wall-coordinate metric, and
connection to TokaMaker field continuity need review before it is wired into
MUG.

## Axisymmetric metric

For axisymmetric geometry, the same sign convention should be retained but the
weak form must carry the volume and surface metric:

```text
Integral_Omega 2*pi*R sigma grad(w).grad(phi) dR dZ
  = Integral_Omega 2*pi*R sigma grad(w).(u x B) dR dZ
    - Integral_boundary 2*pi*R w j_n dS
```

The postprocessing reference currently uses a Cartesian x-z grid. Any
axisymmetric implementation should add a separate manufactured solution with
the `R` metric rather than silently reusing the Cartesian operator.

## Material interfaces

For scalar conductivity jumps on a conforming mesh, the weak form should impose
continuity of potential and normal current weakly:

```text
phi_left = phi_right
J_left.n = J_right.n
```

The Python two-region diagnostics show why this distinction matters. Nodal
central differences produce a large pointwise current artifact at a
discontinuous conductivity interface, even when the conservative face current
is correct to roundoff. A future finite-element implementation should report
interface current using a conservative face/weak-flux diagnostic, not a naive
nodal gradient at the jump.

## Offline extraction path

The added utility
`cases/lm_mhd_coupling/inductionless_csv_postprocess.py` reads a structured
x-z CSV and computes:

- reconstructed current `J = sigma(-grad(phi) + u x B)`,
- structured `div J`,
- wall-normal current extrema,
- optional error against `target_jx`, `target_jy`, `target_jz`.

Demo command:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 cases/lm_mhd_coupling/inductionless_csv_postprocess.py --demo
```

Demo outputs under `cases/lm_mhd_coupling/results/`:

- `inductionless_postprocess_demo_input.csv`
- `inductionless_postprocess_demo_metrics.csv`
- `inductionless_postprocess_demo_current.png`
- `inductionless_postprocess_demo_summary.md`

The default demo reports roundoff-level charge/current recovery:

| metric | value |
|---|---:|
| max interior `|div J|` | `5.551e-16` |
| max wall-normal current | `8.706e-18` |
| max target-current error | `6.939e-18` |

## Standalone Q1 weak-form prototype

The script `cases/lm_mhd_coupling/inductionless_fe_potential_mms.py` assembles
the Cartesian Q1 finite-element weak form directly:

```text
Integral_Omega sigma grad(w).grad(phi) dV
  = Integral_Omega sigma grad(w).(u x B) dV
```

It uses all-insulating natural current closure and an augmented nodal
mean-value gauge. The default manufactured solution gives second-order
potential convergence:

| n | max error | observed order |
|---:|---:|---:|
| 9 | `9.333e-6` | |
| 17 | `2.338e-6` | `1.997` |
| 25 | `1.040e-6` | `1.999` |
| 33 | `5.849e-7` | `1.999` |
| 49 | `2.600e-7` | `1.999` |

This confirms the weak-form sign and gauge for a scalar smooth-conductivity
case before any Fortran field-structure change. The script also recovers
boundary-face `J.n` from Q1 element gradients; the finest-grid max is
`5.827e-4` and decreases approximately first order. Material-interface flux
recovery remains a next step.

The companion script
`cases/lm_mhd_coupling/inductionless_fe_two_region_flux_mms.py` covers that
next interface check for a face-aligned scalar conductivity jump. Because the
piecewise-linear exact potential is in the conforming Q1 space, the default run
recovers the interface current to roundoff:

| n | max potential error | max interface target error | max interface jump |
|---:|---:|---:|---:|
| 17 | `5.551e-16` | `4.270e-15` | `3.109e-15` |
| 33 | `2.554e-15` | `1.271e-14` | `7.661e-15` |
| 49 | `1.238e-14` | `6.240e-14` | `8.882e-15` |

The script
`cases/lm_mhd_coupling/inductionless_axisymmetric_fe_potential_mms.py` adds
the scalar axisymmetric `R` metric to the weak form:

```text
Integral_Omega R sigma grad(w).grad(phi) dR dZ
  = Integral_Omega R sigma grad(w).(u x B) dR dZ
```

It uses an axisymmetrically divergence-free manufactured current. The default
run gives second-order potential convergence:

| n | max error | observed order |
|---:|---:|---:|
| 9 | `1.943e-5` | |
| 17 | `4.886e-6` | `1.992` |
| 25 | `2.175e-6` | `1.995` |
| 33 | `1.225e-6` | `1.997` |
| 49 | `5.448e-7` | `1.998` |

The axisymmetric utility layer now also includes an explicit toroidal electric
field term in the reconstructed current:

```text
J = sigma(-grad(phi) + u x B + E_phi e_phi)
```

The diagnostic
`cases/lm_mhd_coupling/axisymmetric_inductionless_toroidal_feedback.py`
uses this path to create a `CouplingSnapshot` and estimate magnetic feedback
from `J_phi`. In the default selected case, `E_phi = 0.02 V/m` gives
max `|delta B|/|B_ref| = 1.656e-3`, classified as
`two-way-candidate`. This is a manufactured current diagnostic, not a coupled
TokaMaker solve, but it is the clearest current path to preserve in future
weak-form work.

The companion diagnostic
`cases/lm_mhd_coupling/axisymmetric_loop_voltage_feedback.py` derives the
toroidal electric field from `E_phi(R) = V_loop/(2*pi*R)`. The default selected
case, `V_loop = 0.2 V`, gives max `|delta B|/|B_ref| = 1.909e-3`, again
classified as `two-way-candidate`.

The same axisymmetric reference layer now includes global integral checks using
`2*pi*R*dR*dZ` control-volume weights:

```text
P_mech = Integral f.u dV
P_Joule = Integral |J|^2/sigma dV
T_axis = Integral R f_phi dV
```

For the manufactured default cases, the constant-`E_phi` drive gives
`P_Joule = 3.022e3 W`, `P_mech = 2.969e0 W`, and
`T_axis = 9.311e0 N m`; the `V_loop = 0.2 V` drive gives
`P_Joule = 3.513e3 W`, `P_mech = -4.551e-1 W`, and
`T_axis = 9.311e0 N m`. These are diagnostic integrals of prescribed fields.
A future weak-form implementation should reproduce these integral conventions
before using force density outputs for two-way momentum or thermal coupling.

The boundary-current closure audit uses signed outward surface integrals:

```text
I_rmin = Integral -J_R 2*pi*R_min dZ
I_rmax = Integral  J_R 2*pi*R_max dZ
I_zmin = Integral -J_Z 2*pi*R dR
I_zmax = Integral  J_Z 2*pi*R dR
```

Both selected manufactured axisymmetric cases currently give zero max
integrated boundary current and zero net integrated boundary current. A future
finite-element implementation should preserve this integrated current balance
even when pointwise nodal-gradient diagnostics are noisy near material
interfaces.

For future MUG output, export or interpolate fields to a documented structured
CSV with:

```text
x,z,u_x,u_z,sigma,phi
```

Optional columns are:

```text
u_y,b_x,b_y,b_z,target_jx,target_jy,target_jz
```

For HDF5 handoff files, use
`cases/lm_mhd_coupling/coupling_snapshot_audit.py` on the exported
`CouplingSnapshot` before feeding it into a coupled workflow. The audit checks
the `J_phi dA` feedback metric, detects complete structured `R-Z` grids, and
applies the signed integrated boundary-current closure test. The default
loop-voltage snapshot audit currently reports max `|delta B|/|B_ref| =
1.909e-3`, `two-way-candidate`, and zero integrated boundary leakage.

If extracting from `xmhd_2d.F90`, keep these conventions visible:

```text
B = grad(psi) x yhat + by*yhat + B_0
B_x = -d(psi)/dz
B_y = by + B_0_y
B_z =  d(psi)/dx
```

Do not assume the `eta` input is an electrical resistivity without checking the
normalization. The resolved Hartmann work used the MUG mapping
`Ha = B0*a/sqrt(mu0*eta*rho*nu)`, which should be written beside any extracted
`sigma` value.

The utility `cases/lm_mhd_coupling/mug_restart_manifest.py` was added to check
export readiness before guessing. On the default resolved-Hartmann Ha=2
restart it reports:

| item | measured value |
|---|---:|
| restart time | `0` |
| restart state length | `644` |
| profile sidecar rows | `725` |
| restart/profile direct index match | `False` |
| coordinates/connectivity in restart | `False` |
| electric potential in restart | `False` |

So the current restart alone is not sufficient for the inductionless CSV
postprocessor. A real exporter needs final-time state fields together with
mesh/plot coordinates, and it cannot fabricate `phi` because `xmhd_2d` is not
currently an inductionless electric-potential solve.

## Implementation ladder

1. Keep using the Python MMS suite as the sign and operator referee.
2. Add an offline MUG field extraction script that produces the CSV schema
   above from existing plot/HDF5 output, without changing solver physics.
3. Put loop-voltage-driven `E_phi(R)` into the axisymmetric weak-form prototype,
   and locate/generate MUG plot output containing final-time state fields plus
   coordinates/connectivity. Do not fabricate missing electric potential data.
4. Add tests matching the existing Python references: smooth scalar
   conductivity, nonzero Neumann flux, full `sigma u x B` MMS, and two-region
   interface current.
5. Only after those tests are stable, decide whether the potential field belongs
   in standalone MUG, the coupled TokaMaker/blanket layer, or both.

## Red-team checklist

- Verify the sign of `u x B`, `grad(phi)`, `J x B`, and boundary normals with a
  manufactured solution before any solver result is trusted.
- Check the Neumann compatibility condition and gauge in every all-insulating
  case.
- Confirm Joule heating is non-negative and drag/Lorentz work has the expected
  sign in energy diagnostics.
- Do not report pointwise interface current from nodal gradients at
  discontinuous conductivity jumps.
- Keep reduced Hartmann drag off in resolved Hartmann/Shercliff/Hunt
  verification cases.
- Treat two-way plasma/liquid-metal feedback as a regime question: the
  axisymmetric feedback diagnostics suggest one-way studies are reasonable
  below about `1e-3` field perturbation, while `1e-3` to `1e-2` deserves
  coupled sensitivity studies and `~1e-2` is a strong two-way-coupling prompt.
