# Inductionless manufactured-solution diagnostic

This is a Python reference artifact for the future electric-potential LM-MHD
mode. It does not run MUG and does not validate an inductionless solver.

Manufactured setup:
- A streamfunction creates a rectangular-domain current with zero discrete divergence.
- The wall-normal current is zero on all four rectangular boundaries.
- With uniform `B_y` and zero electric potential, a velocity field is chosen so
  `sigma (u x B)` reconstructs the manufactured current.

Headline metrics:
- sigma = 8.000e+05 S/m.
- B_y = 5.000e+00 T.
- max |J| = 1.000e+05 A/m^2.
- max |u| = 2.501e-02 m/s.
- max interior |div J| = 2.328e-08 A/m^3.
- max wall-normal current = 1.225e-11 A/m^2.

Generated files:
- `inductionless_mms_metrics.csv`
- `inductionless_mms_fields.png`
- `inductionless_mms_midplane.png`

Red-team interpretation:
- This checks signs, algebra, and diagnostics only.
- It does not solve the electric-potential Poisson equation.
- It does not include material jumps, conducting walls, or finite-element operators.
- The next implementation step is a real Poisson reference solve with Neumann/Robin
  wall-current closure and manufactured finite-element residual checks.
