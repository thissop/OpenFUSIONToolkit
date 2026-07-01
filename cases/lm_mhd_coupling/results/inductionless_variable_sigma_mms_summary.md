# Variable-conductivity Dirichlet potential MMS

This is a Python reference diagnostic for coefficient handling in the future
inductionless electric-potential mode. It does not run MUG.

Manufactured setup:
- `sigma = 1 + 0.3 x + 0.2 z`.
- `phi = sin(pi x) sin(pi z)`.
- `source = div(sigma grad(phi))` evaluated analytically.
- Dirichlet boundary values are zero.

Finest-grid result:
- n = 65.
- max error = 2.032e-04.
- rms error = 1.000e-04.

Generated files:
- `inductionless_variable_sigma_mms_metrics.csv`
- `inductionless_variable_sigma_mms_solution.png`
- `inductionless_variable_sigma_mms_convergence.png`

Red-team interpretation:
- Smooth scalar conductivity only.
- Dirichlet boundary values only.
- No material jump/interface condition or Robin wall conductance.
- Uniform-grid finite differences only, not finite elements.
