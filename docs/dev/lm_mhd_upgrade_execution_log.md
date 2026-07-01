# LM-MHD upgrade execution log

This log records the 2026-07-01 continuation of the LM-MHD upgrade work on
branch `lm-mhd-upgrades`.

## 2026-07-01

- Read `deliverables/lm_mhd_for_yuchen/OFT_LM_MHD_state_and_roadmap.md` and the
  attached "Upgrade roadmap, ranked by effort to impact" excerpt.
- Checked branch and status first. Branch was `lm-mhd-upgrades`. The tree was
  already dirty with prior Claude/user Hartmann artifacts and an uncommitted
  `src/physics/xmhd_2d.F90` conducting-wall stub.
- Read `AGENTS.md`, `docs/dev/lm_mhd_upgrade_plan.md`,
  `docs/dev/lm_mhd_validation_plan.md`, TokaMaker Python interfaces,
  `gs_save_mug`, and the existing Hartmann scripts.
- Incorporated Sophia Guizzo's clarification from the user-supplied screenshots:
  online standalone 2D MUG is compressible; incompressibility is in Sophia's
  pending merge; wall-current closure in coupled simulations comes through
  TokaMaker/FE continuity; reduced drag is mainly relevant for full-device
  coupled meshes where Hartmann layers cannot all be resolved.
- Decided not to edit Fortran solver physics in this pass. Rationale: high risk
  of colliding with pending incompressible merge and existing dirty Fortran
  artifacts; better immediate value from executable triage and exchange
  contracts.
- Added pure-Python `OpenFUSIONToolkit.LM_MHD` package with:
  - dimensionless groups: `Re`, `Rm`, `Ha`, interaction parameter `N`;
  - time scales: magnetic diffusion and flow time;
  - upper-bound sheet-current feedback estimate;
  - `CouplingSnapshot` HDF5 schema for future MUG-to-TokaMaker current feedback.
- Added `src/tests/physics/test_lm_mhd_coupling.py`.
- Added `cases/lm_mhd_coupling/two_way_coupling_scan.py`.
- First test attempt failed because importing the source-tree top-level
  `OpenFUSIONToolkit` package requires a compiled `liboftpy` shared library in
  this Mac-side shell. This was a packaging/import issue, not a formula failure.
- Added a source-tree fallback loader in the test and scan script; fixed a
  Python 3.13 dataclass loader issue by registering the fallback module in
  `sys.modules`.
- Verified:

```text
PYTHONPATH=src/python pytest -q src/tests/physics/test_lm_mhd_coupling.py
...... [100%]
6 passed in 0.38s
```

- Regenerated smplotlib artifacts:

```text
python3 cases/lm_mhd_coupling/two_way_coupling_scan.py
Wrote LM-MHD coupling scan artifacts to .../cases/lm_mhd_coupling/results
```

- Generated files:
  - `cases/lm_mhd_coupling/results/coupling_regime_scan.csv`
  - `cases/lm_mhd_coupling/results/pbli_regime_heatmap.png`
  - `cases/lm_mhd_coupling/results/li_regime_heatmap.png`
  - `cases/lm_mhd_coupling/results/two_way_threshold_lines.png`
  - `cases/lm_mhd_coupling/results/summary.md`

## Honest limitations

- No MUG Fortran physics was changed in this pass.
- No new MUG or TokaMaker coupled simulation was run.
- The coupling scan is an order-of-magnitude triage artifact using rounded
  representative material properties. It is not a solver result and not a design
  claim.
- Existing dirty Hartmann result files and the existing `xmhd_2d.F90` stub were
  left untouched except for being observed in status.
