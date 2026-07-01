# LM-MHD two-way feedback diagnostics

This note records the second autonomous LM-MHD upgrade slice from 2026-07-01.
It adds analytic/current-geometry infrastructure for estimating whether liquid
metal blanket currents could produce a magnetic perturbation large enough to
justify a reverse coupling path into TokaMaker. It is not a MUG result, not a
TokaMaker result, and not a plasma-equilibrium response calculation.

## What changed

- Added axisymmetric circular-loop Biot-Savart utilities to
  `OpenFUSIONToolkit.LM_MHD.coupling`.
- Added helpers for converting axisymmetric cell volumes to poloidal area
  weights, integrating `J_phi dA`, computing current centroids, and returning
  reduced feedback metrics.
- Added unit tests that check the on-axis loop-field closed form, singular
  probe rejection, linearity in current, signed versus absolute current
  integration, and the axisymmetric volume-to-area convention.
- Added `cases/lm_mhd_coupling/axisymmetric_feedback_scan.py`, which creates a
  manufactured outboard blanket current patch, writes/reads it through the
  existing HDF5 `CouplingSnapshot` schema, and saves smplotlib plots.

## Results generated

Artifacts are in `cases/lm_mhd_coupling/results/`:

- `axisymmetric_feedback_snapshot.h5`
- `axisymmetric_feedback_scan.csv`
- `axisymmetric_feedback_field_map.png`
- `axisymmetric_feedback_scaling.png`
- `axisymmetric_feedback_summary.md`

For the baseline manufactured current patch, with peak `J_phi = 1.0e5 A/m^2`
and `B_ref = 5 T`, the scan reports:

| quantity | value |
|---|---:|
| integrated `|J_phi| dA` | `2.163e4 A` |
| max probe `|delta B|` | `1.288e-2 T` |
| max probe `|delta B|/|B_ref|` | `2.576e-3` |
| classification | `two-way-candidate` |

The amplitude scan intentionally crosses the thresholds: `1e4` and `3e4 A/m^2`
peak current density remain `weak-feedback`, `1e5` and `3e5 A/m^2` are
`two-way-candidate`, and `1e6` and `3e6 A/m^2` are
`strong-two-way-candidate`.

## Interpretation

Two-way plasma/liquid-metal coupling is not automatically important for every
blanket MHD calculation. If induced or motional currents generate only
`|delta B|/|B_ref| << 1e-3` at plasma-relevant locations, one-way imposed-field
LM-MHD is a defensible first calculation. If the ratio is around `1e-3` or
larger, the reverse magnetic perturbation is large enough to deserve a coupled
sensitivity study. Ratios around `1e-2` should be treated as strong prompts for
two-way coupling, subject to a real current-closure and equilibrium solve.

This is why two-way coupling is strategically significant but not the very first
implementation target. Effective drag and inductionless current closure make
the blanket calculation tractable; the feedback metric says when the resulting
currents need to be passed back to the plasma equilibrium.

## Red-team limits

- The field estimate uses circular current-loop geometry, not the TokaMaker
  Green's function implementation.
- The manufactured currents do not solve Ohm's law, electric potential,
  wall-current closure, or the MUG induction equations.
- The diagnostic ignores vessel currents, conducting structures, iron,
  plasma-response screening/amplification, and equilibrium constraints.
- Probe points are outside the current patch, so this is not a finite-element
  self-field or near-singular integration test.
- Thresholds are engineering triage conventions, not validation criteria.

## Verification

Focused test command:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src/python pytest -q \
  src/tests/physics/test_lm_mhd_coupling.py \
  src/tests/physics/test_lm_mhd_closures.py \
  src/tests/physics/test_lm_mhd_feedback.py
```

Current result: `19 passed`.

The scan command:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 cases/lm_mhd_coupling/axisymmetric_feedback_scan.py
```

The generated PNG files were checked with PIL for nonblank image content and
visually inspected.
