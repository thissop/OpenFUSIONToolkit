# Full-induction / inductionless low-Rm bridge

This reduced 1D channel diagnostic solves the linear full-induction
steady relation for induced `b_x` and compares Ampere current
`(1/mu0) db_x/dz` to the inductionless Ohm's-law current with net
streamwise current constrained to zero. It is not a MUG run.

Selected case:
- half-width `a = 5.000e-02 m`.
- B0 = 5.000e+00 T.
- U = 2.000e-01 m/s.
- sigma = 8.000e+05 S/m.
- Rm = 1.005e-02.
- max |b_x|/B0 = 1.290e-03.
- current bridge Linf relative error = 1.875e-05.
- net current integral = -7.276e-12 A/m.
- wall b_x mismatch = 8.186e-18 T.
- final transient relative error at 3 tau_m = 6.523e-13.

Scan extrema:
- max Rm in scan = 2.011e-01.
- max induced-field ratio in scan = 2.580e-02.

Generated files:
- `full_induction_low_rm_bridge_scan.csv`
- `full_induction_low_rm_bridge_profiles.png`
- `full_induction_low_rm_bridge_scaling.png`
- `full_induction_low_rm_bridge_transient.png`

Red-team interpretation:
- Full induction is the parent model; inductionless is justified only
  when induced-field ratios and Rm are small for the target regime.
- This bridge uses a prescribed 1D flow and linear full-induction
  equation, so it does not validate MUG or nonlinear 3D induction.
- The result is a sign/asymptotic consistency check: in this controlled
  low-Rm setting, Ampere current from induced `b` and Ohm's-law current
  must agree.
