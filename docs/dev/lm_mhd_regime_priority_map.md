# LM-MHD B-U regime priority map

Prepared on branch `lm-mhd-upgrades`, 2026-07-02.

## Scope

This note records a reduced material-property scan over velocity `U` and
magnetic field `B`. It is a planning diagnostic for prioritizing drag,
inductionless, and two-way feedback work. It is not a MUG result, not a
TokaMaker response calculation, and not a blanket design prediction.

## Diagnostic

Script:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
  cases/lm_mhd_coupling/lm_mhd_regime_priority_map.py
```

Default assumptions:

- representative length `L = 0.05 m`,
- conducting thickness `t = 0.02 m`,
- feedback reference field `B_ref = 5 T`,
- sheet-current closure factor `0.25`,
- `U` scanned from `0.01` to `1.0 m/s`,
- `B` scanned from `0.5` to `10 T`.

Generated artifacts:

- `cases/lm_mhd_coupling/results/lm_mhd_regime_priority_map.csv`
- `cases/lm_mhd_coupling/results/lm_mhd_regime_priority_summary.csv`
- `cases/lm_mhd_coupling/results/lm_mhd_regime_priority_metrics.png`
- `cases/lm_mhd_coupling/results/lm_mhd_regime_priority_classification.png`
- `cases/lm_mhd_coupling/results/lm_mhd_regime_priority_summary.md`

## Representative Result

Nearest grid point to `U = 0.2 m/s`, `B = 5 T`:

| material | Ha | Re | N | Rm | feedback | class |
|---|---:|---:|---:|---:|---:|---|
| PbLi | `5.465e3` | `5.430e4` | `5.499e2` | `9.827e-3` | `4.966e-4` | `drag-important` |
| Li | `2.131e4` | `1.086e4` | `4.180e4` | `3.931e-2` | `1.986e-3` | `two-way-candidate` |

Classification fractions over the scanned B-U grid:

| material | low | drag-important | two-way-candidate | full-induction-candidate |
|---|---:|---:|---:|---:|
| PbLi | `0.000` | `0.902` | `0.098` | `0.000` |
| Li | `0.000` | `0.655` | `0.195` | `0.150` |

## Interpretation

For these default assumptions, Lorentz-force/drag physics dominates a broad
region for both representative liquid metals. Two-way magnetic feedback is not
automatically the leading correction for PbLi at the representative point:
the sheet estimate is below `1e-3`, so one-way imposed-field studies can be a
defensible first pass if the current closure is similar. Lithium is more
aggressive: higher conductivity pushes the representative point into
`two-way-candidate`, and part of the scanned grid crosses the `Rm = 0.1`
full-induction scrutiny threshold.

This supports the current priority order:

1. Keep resolved Hartmann/Shercliff verification and drag/underresolution
   infrastructure clean.
2. Build inductionless current closure and conservative current audits.
3. Preserve a two-way current-to-field path, especially for higher-conductivity
   liquids, larger velocities, larger thicknesses, or stronger current closure.
4. Treat full induction as a regime-specific concern rather than the default
   for PbLi-like low-`Rm` cases.

## Red-Team Limits

- Feedback is a sheet-current scale with a closure factor, not a solved current
  path.
- The material properties are rounded representative values.
- The classification uses planning thresholds, not validation criteria.
- Geometry, wall conductance, vessel currents, plasma response, and blanket
  segmentation are absent.
- The map is useful for choosing what to simulate, not for claiming a device
  response.

## Verification

The script was rerun after implementation. Both saved PNGs were checked with
PIL for nonblank image content. The formulas used by the script are covered by
the focused LM-MHD coupling unit tests.
