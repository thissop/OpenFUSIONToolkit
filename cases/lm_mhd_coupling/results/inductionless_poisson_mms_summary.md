# Dirichlet electric-potential Poisson MMS

This is a scalar finite-difference reference solve for future inductionless
LM-MHD electric-potential work. It is not a MUG result.

Equation and manufactured solution:
- Solve `laplacian(phi) = source` on `[0,1] x [0,1]`.
- `phi = sin(pi x) sin(pi z)`.
- `source = -2 pi^2 phi`.
- Dirichlet boundary values are zero on all sides.

Finest-grid result:
- n = 65.
- h = 1.562e-02.
- max error = 2.008e-04.
- rms error = 9.887e-05.

Generated files:
- `inductionless_poisson_mms_metrics.csv`
- `inductionless_poisson_mms_solution.png`
- `inductionless_poisson_mms_convergence.png`

Red-team interpretation:
- This verifies a uniform-grid Dirichlet scalar reference solve only.
- It does not implement insulating Neumann or conducting-wall Robin current closure.
- It does not use finite elements and should be treated as a manufactured reference target.
