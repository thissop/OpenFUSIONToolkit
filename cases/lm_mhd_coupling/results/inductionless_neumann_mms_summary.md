# Homogeneous Neumann electric-potential Poisson MMS

This verifies a mean-gauged scalar finite-difference reference solve with
zero normal derivative on all rectangle walls. It is not a MUG result.

Equation and manufactured solution:
- Solve `laplacian(phi) = source` on `[0,1] x [0,1]`.
- `phi = cos(pi x) cos(pi z)`, which has homogeneous Neumann walls.
- `source = -2 pi^2 phi`, whose grid mean is zero to roundoff.
- The nullspace is fixed by enforcing mean(phi) = 0.

Finest-grid result:
- n = 65.
- h = 1.562e-02.
- max mean-free error = 2.008e-04.
- rms mean-free error = 1.020e-04.
- solved mean(phi) = 1.337e-16.

Generated files:
- `inductionless_neumann_mms_metrics.csv`
- `inductionless_neumann_mms_solution.png`
- `inductionless_neumann_mms_convergence.png`

Red-team interpretation:
- This is homogeneous Neumann only: no nonzero motional-wall Neumann data.
- It has no Robin conducting-wall closure and no material interfaces.
- It is a finite-difference reference, not the finite-element equation needed in MUG.
