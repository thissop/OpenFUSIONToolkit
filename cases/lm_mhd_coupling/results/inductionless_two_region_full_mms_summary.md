# Two-region full inductionless MMS

This is a Python finite-difference reference diagnostic connecting a
face-aligned conductivity jump, motional-EMF source generation,
insulating wall flux, Neumann potential recovery, reconstructed current,
and conservative interface-current checks. It does not run MUG.

Headline metrics:
- max |phi error| = 2.759e-14.
- max bulk |J - J_target| = 1.255e-13.
- max bulk |div J| = 1.858e-12.
- max wall-normal current = 3.264e-14.
- max conservative interface-face current error = 1.228e-13.
- max all-point pointwise current error = 5.250e-01.

Generated files:
- `inductionless_two_region_full_mms_metrics.csv`
- `inductionless_two_region_full_mms_potential.png`
- `inductionless_two_region_full_mms_current.png`
- `inductionless_two_region_full_mms_interface.png`

Red-team interpretation:
- The interface is face-aligned and exactly vertical.
- Bulk pointwise current errors exclude a narrow interface band because
  nodal central differences are not a conservative current operator at
  discontinuous conductivity jumps.
- The interface claim is therefore based on a conservative face-current
  reconstruction using the same harmonic-face conductivity as the
  reference potential operator.
- No tensor conductivity, curved material boundary, conducting-wall
  Robin closure, finite-element weak form, or MUG solve is included.
