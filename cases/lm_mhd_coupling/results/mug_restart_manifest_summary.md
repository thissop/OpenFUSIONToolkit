# MUG restart export-readiness manifest

Restart inspected: `/Users/tkiker/Documents/GitHub/OpenFUSIONToolkit/cases/hartmann/mug/Ha_2/xmhd2d_00000.rst`

Datasets:

| name | shape | dtype | max_abs |
|---|---:|---|---:|
| `OFT_idx_Version` | `1` | `int32` | `1.000e+00` |
| `U_T` | `644` | `float64` | `1.000e+00` |
| `U_by` | `644` | `float64` | `0.000e+00` |
| `U_n` | `644` | `float64` | `1.000e+00` |
| `U_psi` | `644` | `float64` | `0.000e+00` |
| `U_velx` | `644` | `float64` | `0.000e+00` |
| `U_vely` | `644` | `float64` | `0.000e+00` |
| `U_velz` | `644` | `float64` | `0.000e+00` |
| `dt` | `1` | `float64` | `2.000e-02` |
| `t` | `1` | `float64` | `0.000e+00` |

Profile sidecar:
- path: `/Users/tkiker/Documents/GitHub/OpenFUSIONToolkit/cases/hartmann/mug/Ha_2/hartmann_mug.profile`
- exists: `true`
- rows: `725`
- coordinate span: x=[-0.25, 0.25], z=[-1, 1]

Export-readiness conclusion:
- coordinates/connectivity present in restart: `False`.
- electric potential field present in restart: `False`.
- restart time `t = 0`.
- restart state size = `644`, profile rows = `725`, direct index match = `False`.
- The restart contains xmhd_2d state fields (`U_n`, velocity, `U_T`,
  `U_psi`, `U_by`) but not the electric potential required by the
  inductionless CSV postprocessor.
- The inspected Hartmann restart appears to be an initial or checkpoint
  state, not the final developed-flow profile used for comparison.
- Coordinates for this Hartmann case are available only through the
  profile sidecar, not the restart itself, and the sidecar row count
  does not match the restart state vector length.
- A real MUG exporter should therefore read mesh/plot output together
  with restart/state fields and should not fabricate `phi` or `sigma`
  without a documented normalization.
