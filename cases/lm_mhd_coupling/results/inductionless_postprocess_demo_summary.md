# Offline inductionless CSV diagnostics

Input CSV: `/Users/tkiker/Documents/GitHub/OpenFUSIONToolkit/cases/lm_mhd_coupling/results/inductionless_postprocess_demo_input.csv`

This postprocesses structured x-z fields into inductionless current
diagnostics. It does not solve for the potential and does not run MUG.

Metrics:
- max_abs_current = 6.288e-02.
- max_abs_divJ_interior = 5.551e-16.
- rms_divJ_interior = 1.240e-16.
- max_wall_normal_current = 8.706e-18.
- max_abs_current_error = 6.939e-18.
- rms_current_error = 1.566e-18.

Expected CSV columns:
- Required: `x`, `z`, `u_x`, `u_z`, `sigma`, `phi`.
- Optional: `u_y`, `b_x`, `b_y`, `b_z`, `target_jx`, `target_jy`, `target_jz`.

Red-team interpretation:
- Structured uniform-grid diagnostic only.
- Not a finite-element projection and not a conservative interface-current
  operator at material jumps.
- For MUG use, first export/interpolate fields onto a documented structured
  grid and keep the interpolation script with the result.
