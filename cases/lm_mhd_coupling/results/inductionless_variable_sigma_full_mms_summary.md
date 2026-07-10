# Full variable-conductivity inductionless MMS

This is a Python reference diagnostic connecting variable conductivity,
motional EMF source generation, insulating wall flux, potential solve, and
current reconstruction. It does not run MUG.

Headline metrics:
- max |phi error| = 3.587e-07.
- max |J - J_target| = 5.390e-07.
- max interior |div J| = 2.260e-07.
- max wall-normal current = 5.390e-07.

Generated files:
- `inductionless_variable_sigma_full_mms_metrics.csv`
- `inductionless_variable_sigma_full_mms_potential.png`
- `inductionless_variable_sigma_full_mms_current.png`

Red-team interpretation:
- Smooth scalar conductivity only.
- Uniform-grid finite differences only.
- No material jump in the full inductionless reconstruction.
- No conducting-wall Robin closure.
