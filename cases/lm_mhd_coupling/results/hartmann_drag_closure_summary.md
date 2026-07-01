# Hartmann effective-drag closure diagnostics

This is an analytic/reduced-model infrastructure artifact. It is not a MUG run,
not a TokaMaker coupling result, and not a claim of standalone solver validation.

Purpose:
- Identify when a full-device mesh under-resolves delta_H = a/Ha.
- Provide the mean-flow equivalent alpha_eff for an insulating Hartmann channel.
- Exercise the CouplingSnapshot HDF5 path with a manufactured current distribution.

Assumptions:
- Channel half-width a = 0.05 m.
- Kinematic viscosity nu = 1.8e-07 m^2/s.
- Resolved threshold = 8 cells per Hartmann layer.
- Marginal means half the target threshold to the target threshold.
- Reduced drag is intended for full-device underresolved-layer meshes, not for
  replacing resolved standalone Hartmann verification.

Resolution matrix counts:

| label | count |
|---|---:|
| resolved | 5202 |
| marginal | 1270 |
| underresolved | 5288 |

Manufactured current-feedback diagnostic:

- HDF5 file: `manufactured_current_feedback.h5`
- max |J_phi| = 2.493e+06 A/m^2
- max sheet-estimate |delta B|/|B_ref| = 1.567e-02
- triage classification = `two-way-candidate`

Generated plots:
- `hartmann_layer_resolution_map.png`
- `hartmann_effective_drag_alpha.png`
- `drag_closure_validity_matrix.png`
- `manufactured_current_feedback.png`

Red-team notes:
- alpha_eff is calibrated to the mean velocity/forcing balance only.
- It must not be used to claim pointwise boundary-layer profile agreement.
- A feedback flag from the manufactured current is only a schema/triage check.
