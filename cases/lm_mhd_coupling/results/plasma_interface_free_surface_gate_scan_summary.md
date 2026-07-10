# Plasma/LM Free-Surface Gate Scan

This is a linear capillary-gravity-magnetic tractability diagnostic, not a nonlinear free-surface simulation.

| quantity | selected value |
|---|---:|
| material | `Li` |
| electron density | `9.713e+18 m^-3` |
| electron temperature | `2.000e+01 eV` |
| wavelength | `2.045e-02 m` |
| magnetic field | `5.000e+00 T` |
| layer thickness | `2.000e-03 m` |
| normal stress | `3.112e+01 Pa` |
| surface stiffness | `6.112e+09 Pa/m` |
| magnetic stiffness | `6.112e+09 Pa/m` |
| static displacement | `5.092e-09 m` |
| displacement / layer | `2.546e-06` |
| wave period | `1.025e-04 s` |
| explicit dt gate | `3.125e-07 s` |
| magnetic damping ratio | `1.305e+00` |
| thermal penetration / layer | `4.880e-01` |
| cells per wavelength | `2.045e+02` |
| cells across layer | `1.000e+02` |
| full gate pass fraction | `0.775` |

Gate interpretation:
- `linear_ok` requires `|eta|/h` below the configured small-displacement limit.
- `surface_resolved` requires enough cells per surface wavelength.
- `normal_resolved` requires enough cells across the liquid layer.
- `timestep_ok` requires the explicit wave/damping step to exceed the configured practical floor.

The magnetic stiffness scale uses `coefficient * k B^2 / mu`; the coefficient is a model input because exact prefactors depend on geometry and electromagnetic boundary conditions.

Verification backing: `src/tests/physics/test_lm_mhd_plasma_interface.py`.
