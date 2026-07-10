# Axisymmetric LM-MHD global integral diagnostics

Prepared on branch `lm-mhd-upgrades`, 2026-07-02.

## Scope

This note records a reduced Python diagnostic layer for integrating LM-MHD
force, power, and torque quantities on structured axisymmetric `R-Z` grids.
It is not a MUG result, not a TokaMaker result, and not a blanket loading
prediction. The intent is to make future two-way coupling audits harder to
fool: local `|J x B|` maxima and magnetic-feedback ratios must be accompanied
by global work, Joule power, and torque checks.

## Definitions

The added utilities in `OpenFUSIONToolkit.LM_MHD.inductionless` use the
existing sign convention:

```text
J = sigma(-grad(phi) + u x B + E_phi e_phi)
f = J x B
q_J = |J|^2 / sigma
```

The new integral diagnostics are:

```text
p_mech = f.u
p_elec = J.E
tau_axis = R f_phi
dV = 2*pi*R*dR*dZ
P_mech = Integral p_mech dV
P_elec = Integral p_elec dV
P_J = Integral q_J dV
T_axis = Integral tau_axis dV
```

For inductionless Ohm's law, the pointwise power identity is:

```text
q_J + f.u - J.E = 0
```

This is a sign-convention red-team check tying Joule heating, Lorentz
mechanical work, and electrical/loop-voltage input together.

The functions are:

- `mechanical_power_density(force_density, velocity)`
- `axisymmetric_torque_density(force_density, r, z)`
- `axisymmetric_volume_integral(density, r, z)`
- `axisymmetric_integrated_power(power_density, r, z)`
- `axisymmetric_integrated_torque(force_density, r, z)`

The same continuation adds conservative current-closure audits:

```text
I_rmin = Integral -J_R 2*pi*R_min dZ
I_rmax = Integral  J_R 2*pi*R_max dZ
I_zmin = Integral -J_Z 2*pi*R dR
I_zmax = Integral  J_Z 2*pi*R dR
I_leak = I_rmin + I_rmax + I_zmin + I_zmax
```

The functions are `axisymmetric_boundary_current_integrals(J, r, z)` and
`axisymmetric_net_boundary_current(J, r, z)`.

## Manufactured Cases

The existing axisymmetric feedback scripts now emit the integral columns in
their CSV and Markdown summaries:

- `cases/lm_mhd_coupling/axisymmetric_inductionless_toroidal_feedback.py`
- `cases/lm_mhd_coupling/axisymmetric_loop_voltage_feedback.py`
- `cases/lm_mhd_coupling/axisymmetric_loop_voltage_threshold_scan.py`

Selected constant-`E_phi` case (`E_phi = 0.02 V/m`):

| quantity | value |
|---|---:|
| integrated Joule power | `3.022e3 W` |
| integrated mechanical power `f.u` | `2.969e0 W` |
| integrated electric power `J.E` | `3.025e3 W` |
| integrated power-balance residual | `-2.659e-15 W` |
| max pointwise power-balance residual | `1.705e-13 W/m^3` |
| net axisymmetric torque | `9.311e0 N m` |
| max integrated boundary current | `0.000e0 A` |
| net integrated boundary current | `0.000e0 A` |
| max `|delta B|/|B_ref|` | `1.656e-3` |
| feedback class | `two-way-candidate` |

Selected loop-voltage case (`V_loop = 0.2 V`):

| quantity | value |
|---|---:|
| integrated Joule power | `3.513e3 W` |
| integrated mechanical power `f.u` | `-4.551e-1 W` |
| integrated electric power `J.E` | `3.513e3 W` |
| integrated power-balance residual | `-1.153e-14 W` |
| max pointwise power-balance residual | `3.411e-13 W/m^3` |
| net axisymmetric torque | `9.311e0 N m` |
| max integrated boundary current | `0.000e0 A` |
| net integrated boundary current | `0.000e0 A` |
| max `|delta B|/|B_ref|` | `1.909e-3` |
| feedback class | `two-way-candidate` |

## Red-Team Notes

- These integrals are over prescribed manufactured fields, not solved MUG
  fields.
- The selected constant-`E_phi` and loop-voltage cases have the same net
  torque because the manufactured poloidal-current contribution dominates the
  azimuthal Lorentz component in this setup. That is a property of this
  manufactured case, not a general LM blanket conclusion.
- The sign of `P_mech = Integral f.u dV` is convention-sensitive and useful:
  the selected loop-voltage case gives a small negative value even while local
  `|J x B|` and Joule heating are positive.
- The power identity closes to roundoff only because the same electric-field
  convention is used in current reconstruction and power accounting. Future
  MUG handoff fields should be checked against the same identity when `E` is
  available.
- Joule power is strictly non-negative pointwise, but numerical integration
  does not by itself validate the conductivity, wall closure, or current path.
- These diagnostics should be rerun on any future MUG-exported current and
  velocity fields before claiming two-way force, power, or thermal coupling.
- A small pointwise `div J` residual does not by itself prove conservative
  closure. The boundary-current integral is the explicit signed leakage audit
  for exported axisymmetric current fields.

## Verification

Focused inductionless tests after adding the integral utilities:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src/python pytest -q \
  src/tests/physics/test_lm_mhd_inductionless.py
```

Result: `26 passed`.

Full focused LM-MHD Python suite:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src/python pytest -q \
  src/tests/physics/test_lm_mhd_coupling.py \
  src/tests/physics/test_lm_mhd_closures.py \
  src/tests/physics/test_lm_mhd_feedback.py \
  src/tests/physics/test_lm_mhd_inductionless.py \
  src/tests/physics/test_lm_mhd_full_induction_bridge.py
```

Result: `45 passed`.
Current result after adding signed boundary-current integrals: `46 passed`.
Current expanded focused result after adding the power-identity audit and
full-induction bridge checks: `53 passed`.

The tests check exact annular volume integration for constant density, vector
dot-product power, axisymmetric torque-density convention, integrated torque,
signed outward boundary-current metrics, net boundary-current cancellation,
and shape rejection. The three axisymmetric diagnostic scripts were rerun, and
the affected saved PNGs were checked with PIL for nonblank image content.
