# LM-MHD Two-Paper Execution Plan

Started: 2026-07-09
Branch: `lm-mhd-upgrades`

This is a working execution plan, not a capability claim. Results are promoted
only after they pass the verification gates below. The reduced wall-drag closure
remains a separate, default-off model; resolved-layer benchmarks must keep it
off unless the specific study is closure calibration.

## Current Measured Baseline

- Local focused LM-MHD Python suite after the L0/L1/L2/L3 plasma-interface
  additions: `65 passed in 0.66s`.
- Local focused LM-MHD Python suite after the L4 linear free-surface gate
  addition: `67 passed in 0.70s`.
- Local plasma-interface reference matrix driver: all five reduced diagnostics
  returned zero; focused suite recheck `67 passed in 0.86s`; local manifest now
  reports 49 files.
- Ginsburg Python/reference baseline after the L0 plasma-interface addition:
  job `8902072`, `59 passed in 2.08s`, Hartmann analytic self-test passed.
- Ginsburg Python/reference baseline after the L0-L3 plasma-interface
  additions: job `8910539`, `65 passed in 38.88s`, Hartmann analytic self-test
  passed, and the remote manifest refreshed to 44 files.
- Ginsburg Python/reference baseline after the L4 linear free-surface gate:
  job `8910555`, `67 passed in 2.08s`, Hartmann analytic self-test passed, all
  five plasma-interface scans ran, and the remote manifest refreshed to 47
  files.
- Ginsburg reference-matrix wrapper baseline: job `8910567`, `67 passed in
  2.40s`, Hartmann analytic self-test passed, all five reference-matrix
  commands returned zero, and the remote manifest refreshed to 49 files.
- Ginsburg resolved-MUG Hartmann probe: Ha=10 packed-grid job `8910589`
  completed successfully (`0:0`) with drag off, `NY=160`, `packing=8.0`,
  `dt=5e-3`, and `NSTEPS=400`. `analyze_mug.py` reported developed-flow
  nonuniformity `1.347e-7`, `Ha_fit/Ha_pred=1.0009`, normalized-profile L2
  `3.406e-4` vs fitted-Ha analytic, and `u_mean/u_c=0.9004` vs analytic
  `0.9002`. The larger `dt=2e-2` probe job `8910576` was cancelled after
  stalling before the first completed timestep.
- Ginsburg resolved-MUG Hartmann probe: Ha=20 packed-grid job `8910691`
  completed successfully (`0:0`) with drag off, `NY=320`, `packing=8.0`,
  `dt=1.25e-3`, and `NSTEPS=1600`. `analyze_mug.py` reported developed-flow
  nonuniformity `1.312e-3`, `Ha_fit/Ha_pred=1.0004`, normalized-profile L2
  `5.548e-5` vs fitted-Ha analytic, and `u_mean/u_c=0.9500` vs analytic
  `0.9500`. Because the x-uniformity is looser than the strict plan gate and
  no mesh/time-step convergence has been done, this is a successful probe, not
  a final convergence-grade high-Ha result.
- Ginsburg resolved-MUG Hartmann Ha=50 bounded feasibility attempts did not
  produce a profile. Job `8910721` (`NY=500`, `packing=12.0`, `dt=2e-4`,
  `NSTEPS=5000`) failed before timestepping with OFT run status 99:
  `Ground point repeated on same domain` in `mesh_global_igrnd` during square
  mesh generation. Retry job `8910734` (`NY=500`, `packing=8.0`, same `dt` and
  `NSTEPS`) was cancelled after 00:08:42 with metadata only, no mesh statistics
  or timestep output, and low CPU accumulation. Ha=50 remains blocked/open for
  this probe family.
- Ginsburg resolved-MUG Hartmann Ha=20 `NY=400`, `packing=8.0` mesh-refinement
  repeat did not produce a profile. Job `8910763` on `g063` was cancelled after
  00:02:23 with metadata only, no mesh/timestep output, and low CPU
  accumulation. Retry job `8910766` with `--exclude=g063` landed on `g090` and
  was cancelled after 00:11:56; it was CPU-active but still emitted only the
  initial branch/revision lines, with no OFT initialization, mesh statistics,
  timestep output, profile, or status file. Ha=20 convergence remains open.
- Ginsburg resolved-MUG Hartmann Ha=20 timestep-refinement repeat also did not
  produce a profile. Job `8911505` (`NY=320`, `packing=8.0`, `dt=6.25e-4`,
  `NSTEPS=3200`, `--exclude=g063,g090`) was cancelled on `g198` after
  00:05:44 with metadata only, initial branch/revision lines only, no OFT
  initialization, no mesh/timestep output, and low CPU accumulation.
- Ginsburg artifact refresh job `8902095` generated the sheath-load scan and
  refreshed the remote manifest to 35 files.
- Ginsburg OFT build/regression baseline: replacement job `8902073` completed
  successfully. OFT built via source-dependency fallback and targeted physics
  regressions passed: `21 passed, 36 skipped` in `2079.55s`.

## Paper 1: Rigorous LM-MHD Verification Capability

Working title:
`Verification of Liquid-Metal MHD Building Blocks in the Open FUSION Toolkit`

The paper should be framed as a verified workflow and reference-infrastructure
paper. It should not claim that the whole MUG liquid-metal program is validated
until full resolved duct, wall-current, and conservation tests pass.

### Minimum Publishable Claim

OFT can reproduce canonical low-Rm/resolved LM-MHD verification problems and
provides tested reference utilities for inductionless current reconstruction,
power accounting, wall-current closure diagnostics, and plasma-equilibrium
feedback handoff metrics.

### Required Verification Matrix

| rung | artifact | pass criterion |
|---|---|---|
| Analytic Hartmann | `cases/hartmann/analytic_hartmann.py` | Poiseuille limit `u_mean/u_c -> 2/3`; high-Ha flow-rate scaling `Q* ~ 1/Ha`; numeric quadrature agreement below `1e-8` on the sweep. |
| OpenFOAM Hartmann cross-check | `cases/hartmann/openfoam/` results | Insulating-wall profile L2 below `2e-3`; `u_mean/u_c` relative error below `2e-3` for `Ha = 1, 5, 10, 20, 50`. |
| MUG resolved Hartmann | `xmhd_2d.F90` driven channel, drag off | Profile L2 below `5e-3`; `Ha_fit / Ha_input - 1` below `2e-3`; x-uniformity below `1e-6`; Hartmann layer resolved by at least 8 wall-normal cells. |
| MUG high-Ha extension | graded/refined wall-normal mesh | Ha=10 and Ha=20 packed-grid probes now pass profile/Ha-mapping checks; repeat/refine for convergence and extend to `Ha = 50` rather than extrapolating from single probes. |
| Hunt/Shercliff duct reference | `cases/hunt/` and/or `cases/shercliff/` | Conducting/insulating wall-current topology distinguishes Hunt side jets from insulating channel; profile and flow-rate errors reported against a trusted analytic or independently solved reference. |
| Inductionless potential MMS | `cases/lm_mhd_coupling/inductionless*_mms.py` | Second-order potential convergence for smooth Q1 cases; conservative interface flux recovery to roundoff for face-aligned conductivity jumps. |
| Axisymmetric metric MMS | `inductionless_axisymmetric_fe_potential_mms.py` | Second-order convergence with the `R` metric; zero integrated boundary leakage for manufactured divergence-free currents. |
| Power identity | `inductionless_power_balance_residual` tests | Pointwise and integrated `q_J + f.u - J.E = 0` to roundoff for manufactured fields. |
| Current handoff audit | `CouplingSnapshot` HDF5 and audit tools | Schema validation, cell-volume consistency, toroidal-current integration, loop-field feedback scaling, and boundary-current leakage checks all pass. |
| Full-induction low-Rm bridge | `full_induction_low_rm_bridge.py` | Ampere current from induced field matches inductionless Ohm current; induced-field amplitude scales linearly with `Rm` in the low-Rm regime. |

### Paper 1 Implementation Tasks

1. Finish the Ginsburg build/regression baseline.
   - Required: `test_alfven_2d.py`, `test_sound_2d.py`,
     `test_drag_decay.py` pass unchanged with all new options default off.
   - Record exact Slurm job IDs, modules, compiler versions, and failure modes.

2. Replace uniform high-Ha MUG Hartmann meshes with wall-normal refinement.
   - Current status: packed-grid Ginsburg probes at `Ha=10` and `Ha=20` meet
     the profile and Ha-mapping checks. `Ha=20` has a measured x-nonuniformity
     of `1.312e-3`, so it should remain labeled as a probe until repeated with
     convergence controls. Two bounded `Ha=50` feasibility attempts failed
     before producing a profile: one mesh-packing failure at `packing=12.0`,
     and one `packing=8.0` setup/first-solve stall cancelled after 00:08:42.
     A `Ha=20`, `NY=400`, `packing=8.0` refinement repeat also stalled before
     OFT initialization/mesh output in this wrapper and was cancelled. A
     `Ha=20`, `NY=320`, `dt=6.25e-4` timestep-refinement repeat showed the
     same pre-initialization stall pattern and was cancelled.
   - Acceptance: `Ha = 50` resolved MUG profile meets the same error gates as
     the lower-Ha points, and `Ha=10`/`Ha=20` are repeated under at least one
     refinement or time-step check.
   - If MUG mesh generation cannot support this cleanly, document it as a
     meshing limitation and keep the paper claim at low/moderate Ha.

3. Implement a MUG duct geometry for out-of-plane axial flow.
   - Required for Hunt/Shercliff MUG comparisons.
   - Separate velocity component mapping from the Task-A in-plane channel.
   - Add conducting Hartmann-wall EM boundary handling only behind explicit
     default-off case/test configuration.

4. Promote existing Python references into reproducible commands.
   - Add a single `cases/lm_mhd_coupling/run_reference_matrix.py` driver or a
     documented command list.
   - Regenerate CSV/PNG/Markdown artifacts and update the manifest hashes.

5. Add a final Paper 1 report table.
   - One row per rung: analytic, OpenFOAM, MUG, MMS, conservation, feedback.
   - Every row must include status, metric, pass/fail, and artifact path.

### Paper 1 Stop Conditions

- Publishable if Hartmann analytic/OpenFOAM/MUG low-to-moderate Ha, MMS,
  conservation, and handoff diagnostics pass with reproducible artifacts.
- Stronger if high-Ha MUG and Hunt/Shercliff MUG pass.
- Do not claim wall-conductance or full free-surface capability from the current
  reduced-diagnostic layer alone.

## Paper 2: Plasma / Liquid-Metal Coupling

Working title:
`A Verified Coupling Ladder for Plasma-Facing Liquid-Metal MHD in OFT`

The paper should be staged around a model ladder. The first honest contribution
is a verified coupling contract and regime map. Full nonlinear free-surface
plasma/LM interaction is tractable on Ginsburg only if the lower rungs pass and
the geometry remains 2D/axisymmetric or modest 3D.

### Coupling Ladder

| level | model | publishable claim if passed |
|---|---|---|
| L0 | Interface load contract | Plasma sheath heat/current/mass/momentum loads are converted into LM boundary data with tested SI units and conservation identities. |
| L1 | One-way plasma-to-LM response | Prescribed plasma heat/current loads drive inductionless LM fields, Joule heating, Lorentz force, and linear surface response with verified energy/current accounting. |
| L2 | LM-to-plasma magnetic feedback | Axisymmetric LM toroidal current snapshots perturb a TokaMaker/GS equilibrium or a reduced Green's-function proxy; feedback is classified by measured `max |delta B|/|B_ref|`. |
| L3 | Linear free-surface fallback | Capillary-gravity/thermocapillary response is solved and verified against analytic static modes and capillary-wave frequencies. |
| L4 | Nonlinear free surface | ALE/VOF or equivalent nonlinear surface evolution with surface tension, gravity, Lorentz forcing, sheath heat flux, and evaporation terms. |
| L5 | Full induction | Only activated when `Rm`, transient time scales, or induced-field amplitude break the inductionless gate. |

### Paper 2 Gates

- Use inductionless LM if `Rm < 0.1` and the full-induction bridge gives
  `max |b|/B0 < 1e-2` for the case family.
- Preserve two-way magnetic feedback if `max |delta B|/|B_ref| >= 1e-3`;
  promote to strong feedback if it exceeds `1e-2`.
- Treat Lorentz coupling as dynamically important if the interaction parameter
  `N >= 1`, even when reverse magnetic feedback is weak.
- Do not run nonlinear free-surface cases before L0-L3 tests pass.
- If L4 is unstable or too expensive, publish L0-L3 as the rigorous paper and
  label L4 as future work or a single exploratory demonstration.

### Paper 2 Verification Tests

| rung | tests |
|---|---|
| L0 interface contract | Bohm speed formula; `Gamma_i = n_e c_s`; `j_i = Z e Gamma_i`; `q = gamma_sh Gamma_i e T_e`; incident mass flux; evaporative heat/mass signs; vectorized inputs; invalid states rejected. |
| L1 energy/current closure | Coupled surface load -> inductionless current -> Joule/mechanical/electric power identity; integrated boundary current leakage zero for insulating cases. |
| L1 thermal response | Heat equation MMS with sheath/evaporation Robin or flux boundary; integrated energy balance closes to tolerance. |
| L2 magnetic feedback | Circular-loop field analytic checks; TokaMaker Green's-function or GS perturbation comparison; current centroid and total-current convergence under mesh refinement. |
| L3 surface response | Static Young-Laplace/capillary-gravity mode amplitude; capillary wave frequency; thermocapillary shear-driven flow benchmark if included. |
| L4 nonlinear surface | Volume conservation, bounded curvature force, capillary-wave convergence, stress-balance residual, and energy budget with evaporation. |
| L5 full induction | Low-Rm bridge recovers inductionless current; finite-Rm transient benchmark against analytic diffusion or a trusted solver. |

### Work Executed In This Pass

Added `src/python/OpenFUSIONToolkit/LM_MHD/plasma_interface.py` and
`src/tests/physics/test_lm_mhd_plasma_interface.py`.

The new tested L0/L1 utilities include:

- Bohm sound speed,
- sheath particle flux,
- incident ion-current density,
- sheath heat-transmission flux,
- ion momentum flux,
- floating/biasable plasma-sheath surface loads,
- evaporative heat and mass accounting,
- linear capillary-gravity stiffness, displacement, and wave frequency,
- representative liquid-metal thermal properties,
- thermal diffusivity and penetration-depth diagnostics,
- lumped surface-layer energy balance,
- semi-infinite constant-heat-flux analytic surface temperature,
- thermal-regime classification by penetration depth,
- surface-current ohmic heat-flux accounting for biased/electrode-style
  interface currents,
- Marangoni shear stress,
- thermocapillary Couette profile and mean velocity for a thin-layer linear
  fallback,
- order-of-magnitude magnetic free-surface stiffness `k B^2 / mu`,
- low-Rm magnetic damping rate `sigma B^2 / rho`,
- linear capillary-gravity-magnetic wave frequency/period,
- conservative explicit time-step gate `safety / max(omega, gamma)`,
- linear free-surface gate metrics for displacement fraction, thermal
  penetration fraction, magnetic damping ratio, and stiffness.

Generated plasma-interface artifacts:

- `cases/lm_mhd_coupling/plasma_interface_sheath_scan.py`
- `cases/lm_mhd_coupling/plasma_interface_thermal_response_scan.py`
- `cases/lm_mhd_coupling/plasma_interface_current_feedback_scan.py`
- `cases/lm_mhd_coupling/plasma_interface_marangoni_scan.py`
- `cases/lm_mhd_coupling/plasma_interface_free_surface_gate_scan.py`
- `cases/lm_mhd_coupling/run_reference_matrix.py`
- `cases/lm_mhd_coupling/results/plasma_interface_sheath_scan.csv`
- `cases/lm_mhd_coupling/results/plasma_interface_sheath_scan.png`
- `cases/lm_mhd_coupling/results/plasma_interface_sheath_scan_summary.md`
- `cases/lm_mhd_coupling/results/plasma_interface_thermal_response_scan.csv`
- `cases/lm_mhd_coupling/results/plasma_interface_thermal_response_scan.png`
- `cases/lm_mhd_coupling/results/plasma_interface_thermal_response_scan_summary.md`
- `cases/lm_mhd_coupling/results/plasma_interface_current_feedback_scan.csv`
- `cases/lm_mhd_coupling/results/plasma_interface_current_feedback_scan.png`
- `cases/lm_mhd_coupling/results/plasma_interface_current_feedback_scan_summary.md`
- `cases/lm_mhd_coupling/results/plasma_interface_marangoni_scan.csv`
- `cases/lm_mhd_coupling/results/plasma_interface_marangoni_scan.png`
- `cases/lm_mhd_coupling/results/plasma_interface_marangoni_scan_summary.md`
- `cases/lm_mhd_coupling/results/plasma_interface_free_surface_gate_scan.csv`
- `cases/lm_mhd_coupling/results/plasma_interface_free_surface_gate_scan.png`
- `cases/lm_mhd_coupling/results/plasma_interface_free_surface_gate_scan_summary.md`
- `cases/lm_mhd_coupling/results/reference_matrix_runs.csv`
- `cases/lm_mhd_coupling/results/reference_matrix_summary.md`

Measured local result after adding the L3 linear thermocapillary fallback rung:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src/python pytest -q \
  src/tests/physics/test_lm_mhd_coupling.py \
  src/tests/physics/test_lm_mhd_closures.py \
  src/tests/physics/test_lm_mhd_feedback.py \
  src/tests/physics/test_lm_mhd_inductionless.py \
  src/tests/physics/test_lm_mhd_full_induction_bridge.py \
  src/tests/physics/test_lm_mhd_plasma_interface.py
```

Result: `65 passed in 0.66s`.

The refreshed local artifact manifest reports 44 files: 13 CSV files and 18
PNG files with 18/18 nonblank.

Measured local result after adding the L4 linear free-surface tractability gate:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src/python pytest -q \
  src/tests/physics/test_lm_mhd_coupling.py \
  src/tests/physics/test_lm_mhd_closures.py \
  src/tests/physics/test_lm_mhd_feedback.py \
  src/tests/physics/test_lm_mhd_inductionless.py \
  src/tests/physics/test_lm_mhd_full_induction_bridge.py \
  src/tests/physics/test_lm_mhd_plasma_interface.py
```

Result: `67 passed in 0.70s`.

The refreshed local artifact manifest reports 47 files: 14 CSV files and 19
PNG files with 19/19 nonblank. The L4 gate is still a reduced diagnostic: it
quantifies stiffness, displacement scale, thermal penetration, damping, and
mesh/time-step pressure before attempting any nonlinear surface evolution.

After adding the reproducibility driver
`cases/lm_mhd_coupling/run_reference_matrix.py`, all five plasma-interface
artifact commands returned zero locally. The focused suite still passed
(`67 passed in 0.86s`) and the manifest reports 49 files: 15 CSV files and 19
PNG files with 19/19 nonblank.

### Ginsburg Tractability

Tractable now:

- Python reference/MMS/diagnostic sweeps.
- OFT build/regression once dependency setup is stable.
- 2D/axisymmetric inductionless current and feedback scans.
- Linear surface response and reduced thermal/sheath scans.

Tractable with careful scope:

- Resolved 2D MUG Hartmann/Hunt/Shercliff cases with refined meshes.
- Axisymmetric or quasi-2D plasma/LM feedback scans with restartable Slurm
  arrays.
- A modest nonlinear free-surface demonstration if L0-L3 pass first.

Risky for this hardware/time budget:

- Broad 3D nonlinear free-surface MHD parameter studies.
- Fully coupled plasma equilibrium + nonlinear LM surface + evaporation +
  full-induction dynamics in one monolithic solve.

The recommended Paper 2 execution path is therefore L0-L3 first, with L4/L5 as
strictly gated stretch goals.
