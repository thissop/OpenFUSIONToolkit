# Axisymmetric feedback probe-sensitivity red-team scan

Prepared on branch `lm-mhd-upgrades`, 2026-07-02.

## Scope

This note records a reduced probe-location sensitivity scan for the default
loop-voltage `CouplingSnapshot`. It is a red-team check on the circular-loop
feedback metric, not a MUG run, not a TokaMaker solve, and not a plasma
equilibrium sensitivity calculation.

## Diagnostic

Script:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
  cases/lm_mhd_coupling/axisymmetric_feedback_probe_sensitivity.py
```

Input snapshot:

```text
cases/lm_mhd_coupling/results/axisymmetric_loop_voltage_feedback_snapshot.h5
```

The scan varies:

- probe side: `inboard` and `outboard`,
- offset from the current annulus: `0.05, 0.1, 0.2, 0.35, 0.5 m`,
- circular-loop regularization: `0.01, 0.02, 0.03, 0.05, 0.1 m`,
- five vertical probe locations from `Z = -0.4 m` to `Z = 0.4 m`.

Generated artifacts:

- `cases/lm_mhd_coupling/results/axisymmetric_feedback_probe_sensitivity.csv`
- `cases/lm_mhd_coupling/results/axisymmetric_feedback_probe_sensitivity_heatmap.png`
- `cases/lm_mhd_coupling/results/axisymmetric_feedback_probe_sensitivity_curves.png`
- `cases/lm_mhd_coupling/results/axisymmetric_feedback_probe_sensitivity_summary.md`

## Results

The default-like row nearest offset `0.1 m` and regularization `0.03 m`
reports:

| quantity | value |
|---|---:|
| side | `inboard` |
| max `|delta B|/|B_ref|` | `1.968e-3` |
| classification | `two-way-candidate` |

Extremes over all 50 rows:

| case | value |
|---|---:|
| strongest feedback | `2.047e-3` |
| strongest location | `inboard`, offset `0.05 m`, eps `0.01 m` |
| weakest feedback | `2.361e-4` |
| weakest location | `outboard`, offset `0.5 m`, eps `0.1 m` |
| fraction at or above `1e-3` | `0.500` |
| fraction at or above `1e-2` | `0.000` |

## Interpretation

The two-way-candidate conclusion is not just one exact probe point, but it is
also not universal. Inboard and closer probe sets cross the `1e-3` planning
threshold, while distant outboard probes can fall back to `weak-feedback`.
None of the rows reach the `1e-2` strong-candidate threshold.

This supports a conservative prioritization: preserve and test the reverse
current-to-field coupling path, but do not claim strong two-way plasma response
from this reduced model alone. A real TokaMaker response calculation is still
needed for equilibrium impact.

## Red-Team Limits

- The scan varies probe geometry and regularization only; it does not vary
  current closure, vessel currents, plasma response, or blanket geometry.
- Circular-loop regularization is a numerical diagnostic knob, not a physical
  wall or finite-element model.
- The current field comes from a manufactured loop-voltage diagnostic, not MUG.
- Thresholds `1e-3` and `1e-2` are planning conventions, not validation
  criteria.

## Verification

The script was rerun after implementation. The CSV contains `50` rows. Both
saved PNGs were checked with PIL for nonblank image content.
