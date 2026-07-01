# Full inductionless potential MMS

This is a constant-conductivity Python reference diagnostic for the future
electric-potential LM-MHD mode. It does not run MUG.

The manufactured construction sets `u x B = grad(phi) + J_target/sigma`,
where `J_target` is streamfunction-generated, divergence-free, and has zero
normal current at the walls. The scalar solve uses
`laplacian(phi) = div(u x B)` with insulating wall data
`dphi/dn = (u x B).n`.

Headline metrics:
- max |phi error| = 1.102e-16.
- max |J - J_target| = 3.553e-10 A/m^2.
- max interior |div J| = 2.235e-08 A/m^3.
- max wall-normal current = 3.553e-10 A/m^2.

Generated files:
- `inductionless_full_mms_metrics.csv`
- `inductionless_full_mms_potential.png`
- `inductionless_full_mms_current.png`
- `inductionless_full_mms_midplane.png`

Red-team interpretation:
- Constant sigma only; no variable-coefficient or material-interface solve.
- Uniform-grid finite differences only; no finite-element weak form.
- No conducting-wall Robin closure.
- This is a reference target for future solver work, not a solver validation.
