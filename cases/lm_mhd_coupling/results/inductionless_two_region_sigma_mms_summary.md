# Two-region conductivity-jump potential MMS

This is a finite-difference reference diagnostic for material-interface
coefficient handling. It does not run MUG.

Manufactured setup:
- Two conductivities: sigma_left = 1, sigma_right = 4.
- The interface is aligned with an x-grid face at x = 0.5.
- The exact potential is piecewise linear.
- The source is zero.
- Potential and normal current are continuous across the interface.

Result:
- target interface flux = 1.600000e+00.
- max potential error = 6.217e-14.
- max interface-flux error = 1.061e-13.

Generated files:
- `inductionless_two_region_sigma_mms_metrics.csv`
- `inductionless_two_region_sigma_mms_solution.png`
- `inductionless_two_region_sigma_mms_flux.png`

Red-team interpretation:
- The interface is face-aligned and one-dimensional in x.
- This checks harmonic face averaging and flux continuity only.
- It is not a finite-element multi-region implementation.
