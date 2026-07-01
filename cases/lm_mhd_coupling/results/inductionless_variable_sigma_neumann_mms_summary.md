# Variable-conductivity Neumann-flux potential MMS

This is a Python reference diagnostic for smooth variable-conductivity
Neumann flux closure. It does not run MUG.

Manufactured setup:
- `sigma = 1 + 0.3 x + 0.2 z`.
- `phi = (x - 1/2)^2 + (z - 1/2)^2 - mean(phi)`.
- `source = div(sigma grad(phi))` evaluated analytically.
- outward wall flux is `sigma dphi/dn`.
- mean(phi) is fixed to zero.

Finest-grid result:
- n = 65.
- max error = 2.802e-05.
- rms error = 1.046e-05.
- solved mean(phi) = 3.364e-18.

Generated files:
- `inductionless_variable_sigma_neumann_mms_metrics.csv`
- `inductionless_variable_sigma_neumann_mms_solution.png`
- `inductionless_variable_sigma_neumann_mms_convergence.png`

Red-team interpretation:
- Smooth scalar conductivity only.
- Prescribed flux data only; not generated from `u x B` yet.
- No discontinuous material interface in this Neumann diagnostic.
- Uniform-grid finite differences only, not finite elements.
