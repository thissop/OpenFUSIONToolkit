# LM-MHD reference matrix run

Suite: `plasma-interface`

This driver reruns reduced diagnostics and records command status. It is a reproducibility artifact, not a physics validation result.

| command | return code | elapsed [s] | log |
|---|---:|---:|---|
| `plasma_interface_sheath_scan` | `0` | `0.915` | `cases/lm_mhd_coupling/results/reference_matrix_plasma_interface_sheath_scan.log` |
| `plasma_interface_thermal_response_scan` | `0` | `0.872` | `cases/lm_mhd_coupling/results/reference_matrix_plasma_interface_thermal_response_scan.log` |
| `plasma_interface_current_feedback_scan` | `0` | `1.148` | `cases/lm_mhd_coupling/results/reference_matrix_plasma_interface_current_feedback_scan.log` |
| `plasma_interface_marangoni_scan` | `0` | `0.970` | `cases/lm_mhd_coupling/results/reference_matrix_plasma_interface_marangoni_scan.log` |
| `plasma_interface_free_surface_gate_scan` | `0` | `1.229` | `cases/lm_mhd_coupling/results/reference_matrix_plasma_interface_free_surface_gate_scan.log` |

Result: PASS (all commands returned zero).

Run `cases/lm_mhd_coupling/lm_mhd_artifact_manifest.py` after this driver to refresh artifact hashes and PNG sanity metadata.
