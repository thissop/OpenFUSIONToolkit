# Full-induction / inductionless low-Rm bridge

Prepared on branch `lm-mhd-upgrades`, 2026-07-02.

## Scope

This note records a reduced full-induction reference bridge for a 1D channel.
It is not a MUG run and not a blanket prediction. Its purpose is to answer the
modeling question: when the magnetic Reynolds number is small, does the
full-induction parent model recover the same current as the inductionless
Ohm's-law current?

## Model

The bridge uses a prescribed parabolic channel flow:

```text
u = u_x(z) e_x
B = B0 e_z + b_x(z) e_x
```

The inductionless current with net streamwise current constrained to zero is:

```text
J_y = sigma (E_y - u_x B0)
E_y = B0 mean(u_x)
```

The steady linear full-induction equation gives:

```text
0 = B0 du_x/dz + eta_m d^2 b_x/dz^2
J_y = (1/mu0) db_x/dz
eta_m = 1/(mu0 sigma)
```

Thus, in this controlled geometry, Ampere current from the induced magnetic
field must match the inductionless current. The induced-field amplitude scales
like `Rm B0`.

## Diagnostic

Script:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
  cases/lm_mhd_coupling/full_induction_low_rm_bridge.py
```

Generated artifacts:

- `cases/lm_mhd_coupling/results/full_induction_low_rm_bridge_scan.csv`
- `cases/lm_mhd_coupling/results/full_induction_low_rm_bridge_profiles.png`
- `cases/lm_mhd_coupling/results/full_induction_low_rm_bridge_scaling.png`
- `cases/lm_mhd_coupling/results/full_induction_low_rm_bridge_transient.png`
- `cases/lm_mhd_coupling/results/full_induction_low_rm_bridge_summary.md`

## Selected Result

Default selected case:

| quantity | value |
|---|---:|
| half-width `a` | `5.000e-2 m` |
| `B0` | `5.000e0 T` |
| velocity scale `U` | `2.000e-1 m/s` |
| conductivity `sigma` | `8.000e5 S/m` |
| `Rm = mu0 sigma U a` | `1.005e-2` |
| max `|b_x|/B0` | `1.290e-3` |
| Ohm-vs-Ampere current relative error | `1.875e-5` |
| net current integral | `-7.276e-12 A/m` |
| wall `b_x` mismatch | `8.186e-18 T` |

The scan reaches max `Rm = 2.011e-1` for the highest-conductivity/highest-velocity
case, and the induced-field ratio scales linearly with `Rm` in this bridge.

## Interpretation

This is the rigorous answer to the “why inductionless?” question:

- full induction is the parent model,
- inductionless is the low-`Rm` asymptotic model,
- for PbLi-like default values in this bridge, `Rm ≈ 1e-2` and the induced
  field is about `1e-3 B0`,
- the current from full-induction Ampere's law matches the inductionless
  Ohm's-law current to discretization error.

That does **not** mean full induction is unimportant. It means the model choice
must be regime-driven. Higher conductivity, larger velocity, larger length
scale, transients, or nonlinear induced-field dynamics push the problem toward
full induction.

## Red-Team Limits

- This is a 1D prescribed-flow bridge, not a solved MUG flow.
- The nonlinear `u x b` feedback is not active in this geometry; the bridge is
  a low-`Rm` asymptotic/sign check.
- Boundary and wall-current physics are simplified to a net-current-zero
  channel constraint.
- The result justifies the inductionless limit only for regimes with small
  `Rm` and small induced-field ratio.

## Verification

Added focused tests in
`src/tests/physics/test_lm_mhd_full_induction_bridge.py`:

- Ampere current recovers inductionless current.
- induced-field ratio scales linearly with `Rm`,
- transient magnetic diffusion relaxes to the steady induced field.

Focused result:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src/python pytest -q \
  src/tests/physics/test_lm_mhd_full_induction_bridge.py
```

Result: `3 passed`.
