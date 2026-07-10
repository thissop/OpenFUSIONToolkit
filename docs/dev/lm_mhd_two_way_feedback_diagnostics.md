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

## Inductionless continuation

A later continuation added
`cases/lm_mhd_coupling/axisymmetric_inductionless_toroidal_feedback.py`, which
keeps the same feedback metric but generates the current through an
axisymmetric inductionless Ohm's-law reconstruction:

```text
J = sigma(-grad(phi) + u x B + E_phi e_phi)
```

This matters because an explicit toroidal electric field or loop-voltage drive
can create `J_phi`, and `J_phi` is the component that produces poloidal
magnetic perturbations in the circular-loop feedback diagnostic.

For the selected default case, `E_phi = 0.02 V/m`, representative
`sigma0 = 8e5 S/m`, and `B_phi = 5 T`, the diagnostic reports:

| quantity | value |
|---|---:|
| total toroidal current | `1.599e4 A` |
| max `|delta B|/|B_ref|` | `1.656e-3` |
| classification | `two-way-candidate` |
| max interior axisymmetric `|div J|` | `1.858e-10` |
| max wall-normal current | `0.000e0` |

This strengthens the triage message: two-way coupling is still regime
dependent, but plausible inductionless toroidal-current drives can cross the
`1e-3` planning threshold and should be preserved in the coupling workflow.

The companion script `cases/lm_mhd_coupling/axisymmetric_loop_voltage_feedback.py`
derives `E_phi(R) = V_loop/(2*pi*R)`. For `V_loop = 0.2 V`, it reports
max `|delta B|/|B_ref| = 1.909e-3`, also classified as
`two-way-candidate`, with max interior axisymmetric `|div J| = 1.858e-10`.

The threshold scan
`cases/lm_mhd_coupling/axisymmetric_loop_voltage_threshold_scan.py` estimates
that, for the PbLi-like row nearest `sigma0 = 8e5 S/m`, the manufactured
two-way threshold `|delta B|/B_ref = 1e-3` occurs near
`V_loop = 1.046e-1 V`. The strong `1e-2` threshold is not reached up to
`V_loop = 1 V` in that row.

The latest continuation adds global integral audits to these manufactured
axisymmetric cases. For the selected constant-`E_phi` case, the diagnostic
reports integrated Joule power `3.022e3 W`, integrated mechanical power
`2.969e0 W`, and net axisymmetric torque `9.311e0 N m`. For the selected
`V_loop = 0.2 V` case, it reports integrated Joule power `3.513e3 W`,
integrated mechanical power `-4.551e-1 W`, and net torque `9.311e0 N m`.
The sign change in integrated mechanical power is a useful reminder that
local force magnitudes and magnetic-feedback ratios are not substitutes for
power/force-balance checks.

The same continuation adds signed outward boundary-current integrals with the
axisymmetric surface metrics. The selected constant-`E_phi` and loop-voltage
cases both report max integrated boundary current `0.000e0 A` and net
integrated boundary current `0.000e0 A`, matching the manufactured insulating
closure.

The handoff-facing audit script
`cases/lm_mhd_coupling/coupling_snapshot_audit.py` reads a `CouplingSnapshot`
HDF5 file and applies the same reduced feedback plus current-closure checks.
On the default loop-voltage snapshot it detects an `81 x 65` structured grid,
recovers max `|delta B|/|B_ref| = 1.909e-3`, classifies the file as
`two-way-candidate`, and reports zero integrated boundary leakage.

The red-team script
`cases/lm_mhd_coupling/axisymmetric_feedback_probe_sensitivity.py` varies
probe side, distance from the current annulus, and circular-loop
regularization. For the same loop-voltage snapshot, the default-like inboard
row gives max `|delta B|/|B_ref| = 1.968e-3`; the strongest scanned row is
`2.047e-3`, the weakest is `2.361e-4`, 50% of rows cross `1e-3`, and no rows
cross `1e-2`. This keeps the two-way-candidate priority plausible but
explicitly not universal.

The broader B-U material map
`cases/lm_mhd_coupling/lm_mhd_regime_priority_map.py` reinforces that nuance.
At the representative point nearest `U = 0.2 m/s`, `B = 5 T`, default PbLi is
classified `drag-important` with feedback `4.966e-4`, while default Li is
classified `two-way-candidate` with feedback `1.986e-3`. Over the scanned grid,
PbLi is `drag-important` in 90.2% of rows and `two-way-candidate` in 9.8%;
Li is `drag-important` in 65.5%, `two-way-candidate` in 19.5%, and
`full-induction-candidate` in 15.0%.

## Red-team limits

- The field estimate uses circular current-loop geometry, not the TokaMaker
  Green's function implementation.
- The manufactured currents do not solve Ohm's law, electric potential,
  wall-current closure, or the MUG induction equations.
- The later inductionless continuation reconstructs current through Ohm's law,
  but still prescribes the potential and `E_phi`; it is not a MUG or TokaMaker
  solve.
- The global power and torque integrals are computed on the same manufactured
  fields. They are consistency diagnostics, not validated predictions of
  blanket loading or plasma torque.
- The signed integrated boundary-current checks are conservative leakage
  audits; they still depend on the prescribed manufactured current path.
- The `CouplingSnapshot` audit is a handoff-file consistency check. It is not a
  replacement for a TokaMaker response calculation.
- Probe-sensitivity results depend on chosen probe families. They are a
  robustness check on the reduced metric, not a plasma-response bound.
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
