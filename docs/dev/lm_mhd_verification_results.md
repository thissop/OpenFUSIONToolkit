# LM-MHD Verification Results — Hartmann Channel (and Hunt)

**Date:** 2026-06-03 (overnight unattended run)
**Branch:** `lm-mhd-upgrades`
**Environment:** Lima Ubuntu 25.10 aarch64, gfortran 15.2, OpenFOAM v1912 (`mhdFoam`), MPICH 4.2.3.

> **Scope / honesty statement.** This is a **personal verification / learning exercise**.
> It independently reproduces the classic Hartmann (and Hunt) duct-flow benchmarks that the
> team has already done, to learn the workflow and cross-check tools. It does **not** validate
> MUG, claim new MUG capability, or advance its physics. The **reduced wall-drag closure was
> OFF (`use_wall_drag=.FALSE.`) throughout** — the entire point is that a *resolved* channel
> does not need it.

---

## 2026-07-02 continuation — reduced coupling diagnostics

The later autonomous continuation did **not** add new resolved MUG/Hartmann
claims and did **not** touch the reduced-drag setting. It added Python-side
reference infrastructure and red-team diagnostics for inductionless
current-feedback workflows:

- axisymmetric inductionless current reconstruction with explicit
  `E_phi`/loop-voltage drive,
- global Joule-power, mechanical-power, and axisymmetric torque integrals,
- signed outward boundary-current integrals and net leakage checks,
- a reusable `CouplingSnapshot` HDF5 audit tool,
- feedback probe/regularization sensitivity scans,
- B-U material regime maps for drag/Rm/two-way prioritization,
- a reduced full-induction/inductionless low-`Rm` bridge,
- a generated artifact manifest with hashes, CSV row counts, and PNG sanity
  metadata.

Key measured reduced-model results:

| diagnostic | measured result |
|---|---|
| constant `E_phi = 0.02 V/m` axisymmetric current | max `|delta B|/|B_ref| = 1.656e-3`, `two-way-candidate` |
| loop-voltage `V_loop = 0.2 V` axisymmetric current | max `|delta B|/|B_ref| = 1.909e-3`, `two-way-candidate` |
| loop-voltage threshold, `sigma0 = 8e5 S/m` | `1e-3` feedback near `V_loop = 1.046e-1 V`; `1e-2` not reached up to `1 V` |
| selected constant-`E_phi` integrals | `P_J = 3.022e3 W`, `P_mech = 2.969e0 W`, `T_axis = 9.311e0 N m` |
| selected loop-voltage integrals | `P_J = 3.513e3 W`, `P_mech = -4.551e-1 W`, `T_axis = 9.311e0 N m` |
| inductionless power identity | selected cases close `q_J + f.u - J.E` to `~1e-14 W` integrated residual |
| selected boundary-current audit | max integrated boundary current `0.000e0 A`, net `0.000e0 A` |
| `CouplingSnapshot` audit | `81 x 65` grid, zero cell-volume relative error, zero boundary leakage |
| probe-sensitivity scan | 50% of rows cross `1e-3`; 0% cross `1e-2` |
| B-U regime map at `U ≈ 0.2 m/s`, `B ≈ 5 T` | PbLi `drag-important` (`feedback = 4.966e-4`); Li `two-way-candidate` (`feedback = 1.986e-3`) |
| full-induction bridge, PbLi-like selected case | `Rm = 1.005e-2`, `max |b_x|/B0 = 1.290e-3`, Ohm-vs-Ampere current error `1.875e-5` |

Interpretation: two-way plasma/liquid-metal magnetic feedback is
regime-dependent. For PbLi-like default assumptions, Lorentz-force/drag and
inductionless current closure remain the first priorities; reverse coupling
becomes plausible near specific loop-voltage/current-drive, probe, material,
or conductivity regimes and should be preserved as a tested handoff path.

Latest focused validation for the Python LM-MHD utilities:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src/python pytest -q \
  src/tests/physics/test_lm_mhd_coupling.py \
  src/tests/physics/test_lm_mhd_closures.py \
  src/tests/physics/test_lm_mhd_feedback.py \
  src/tests/physics/test_lm_mhd_inductionless.py \
  src/tests/physics/test_lm_mhd_full_induction_bridge.py
```

Result: `53 passed`. `git diff --check` was clean after the latest
continuation. New reduced-diagnostic artifacts are indexed by
`cases/lm_mhd_coupling/results/lm_mhd_artifact_manifest.csv`; the manifest
contains 32 entries, including 14 nonblank PNG files.

Detailed continuation notes:

- `docs/dev/lm_mhd_axisymmetric_inductionless_feedback.md`
- `docs/dev/lm_mhd_axisymmetric_global_integrals.md`
- `docs/dev/lm_mhd_coupling_snapshot_audit.md`
- `docs/dev/lm_mhd_feedback_probe_sensitivity.md`
- `docs/dev/lm_mhd_regime_priority_map.md`
- `docs/dev/lm_mhd_full_induction_bridge.md`
- `docs/dev/lm_mhd_inductionless_reference_plan.md`

---

## 2026-07-09 continuation — two-paper execution and plasma-interface L0/L1

This continuation locked the two-paper plan into an execution-gated form and
started the first conservative Paper 2 implementation rung. It does **not** add
new resolved MUG/Hartmann claims and does **not** change the reduced-drag
scope.

New planning/report artifact:

- `docs/dev/lm_mhd_two_paper_execution_plan.md`

New plasma/liquid-metal interface reference layer:

- `src/python/OpenFUSIONToolkit/LM_MHD/plasma_interface.py`
- `src/tests/physics/test_lm_mhd_plasma_interface.py`
- `cases/lm_mhd_coupling/plasma_interface_sheath_scan.py`
- `cases/lm_mhd_coupling/plasma_interface_thermal_response_scan.py`
- `cases/lm_mhd_coupling/plasma_interface_current_feedback_scan.py`
- `cases/lm_mhd_coupling/plasma_interface_marangoni_scan.py`

The new L0/L1 interface utilities provide tested formula-level accounting for:

- Bohm sound speed,
- sheath particle flux,
- incident ion-current density,
- sheath heat-transmission flux,
- ion momentum flux,
- floating/biasable plasma-sheath surface loads,
- evaporative heat and mass signs,
- linear capillary-gravity stiffness, displacement, and wave frequency,
- representative LM thermal properties,
- thermal diffusivity and penetration depth,
- lumped surface-layer energy closure,
- semi-infinite constant-heat-flux analytic surface temperature,
- thermal-regime classification,
- surface-current ohmic heating for biased/electrode-style interface currents,
- Marangoni shear stress and thermocapillary Couette velocity for a linear
  free-surface fallback.

Measured results from this pass:

| diagnostic | measured result |
|---|---|
| Ginsburg Python/reference baseline before L0 addition | Slurm job `8901355`: `53 passed in 1.86s`; Hartmann analytic self-test passed; artifact manifest refreshed |
| Ginsburg expanded Python/reference baseline after L0 addition | Slurm job `8902072`: `59 passed in 2.08s`; Hartmann analytic self-test passed |
| local expanded LM-MHD Python suite after L0 addition | `59 passed in 0.76s` |
| local expanded LM-MHD Python suite after L1 thermal addition | `62 passed in 1.00s` |
| local expanded LM-MHD Python suite after L2 current-feedback addition | `63 passed in 0.68s` |
| local expanded LM-MHD Python suite after L3 thermocapillary fallback addition | `65 passed in 0.66s` |
| plasma-interface sheath scan | generated CSV/PNG/Markdown artifacts under `cases/lm_mhd_coupling/results/` |
| plasma-interface thermal-response scan | generated CSV/PNG/Markdown artifacts under `cases/lm_mhd_coupling/results/` |
| plasma-interface biased-current feedback scan | generated CSV/PNG/Markdown artifacts under `cases/lm_mhd_coupling/results/`; floating sheath is zero-current by default, biased scans use explicit nonzero current fraction |
| plasma-interface thermocapillary fallback scan | generated CSV/PNG/Markdown artifacts under `cases/lm_mhd_coupling/results/`; linear Marangoni/Couette diagnostic only |
| artifact manifest after L3 thermocapillary scan | local manifest refreshed to 44 entries, 13 CSV, 18 PNG with 18/18 nonblank; remote manifest reached 35 entries before SSH access was lost |
| Ginsburg OFT build/regression baseline | Direct configure reached but missing module METIS; source-dependency fallback built dependencies after a SuperLU `lib64` path workaround. Replacement job `8902073` completed successfully: OFT build PASS, targeted physics regressions `21 passed, 36 skipped` in `2079.55s` |

The plasma-interface additions are deliberately a **boundary-load contract**,
analytic heat-response references, biased-current feedback triage, and linear
surface-response references, not a sheath solver and not a nonlinear
free-surface solver. They are the first tested rungs needed before a
plasma/liquid-metal coupling paper can honestly advance to a real TokaMaker
feedback solve, solved thermal transport, or nonlinear free-surface motion.

---

## 1. One-paragraph summary

Three independent solutions of the **insulating-wall Hartmann channel** were computed and
cross-checked: (1) the analytic profile, (2) the resolved 2D MUG solver (`xmhd_2d.F90`, drag
off, driven by a new default-off body-force term), and (3) OpenFOAM `mhdFoam`. **All three
agree to ~0.1%** on the normalized velocity profile and on the discriminating scalar
`u_mean/u_c` across the Hartmann-number range. For MUG, the Hartmann number *fitted from the
computed profile shape* matches the Hartmann number *predicted from the input parameters*
(`Ha = B0·a/√(μ0·η·ρ·ν)`) to **0.1%**, independently confirming both the geometry mapping and
the parameter scaling. The analytic + mhdFoam legs cover `Ha ∈ {1,5,10,20,50}`; the resolved-MUG
leg is confirmed at `Ha ∈ {2,3,5}` plus packed-grid Ginsburg probes at `Ha=10` and `Ha=20`.
`Ha=20` is still labeled a probe because the x-uniformity is looser than the strict plan gate and
no mesh/time-step convergence has been done; MUG at `Ha=50` remains open after two bounded
feasibility attempts failed before producing a profile; see §6. Task B (Hunt conducting-wall duct)
is addressed with a finite-difference reference; see §7.

**Presentation-quality figures** (smplotlib-styled) are in `cases/hartmann/results/report/`:
`threeway_overlay_Ha5.png` (analytic/mhdFoam/MUG with a residual inset showing ~1e-4 agreement),
`umean_uc_vs_Ha.png` (cross-code scaling of the discriminating scalar), `error_vs_Ha.png`
(both codes' L2 vs analytic, ≲1e-3), and `flowrate_vs_Ha.png` (analytic Q*~1/Ha suppression).

---

## 2. Task A — three-way error table (THE HEADLINE)

Insulating Hartmann channel, profiles normalized by centerline velocity. `u_mean/u_c` is the
discriminating scalar (→2/3 Poiseuille at Ha→0, →1 plug core at Ha→∞). Full machine-readable
table: `cases/hartmann/results/error_table.md`. Overlay plots: `cases/hartmann/results/overlay_Ha*.png`.

| Ha | u_mean/u_c (analytic) | mhdFoam L2 | mhdFoam Linf | mhdFoam u_mean/u_c (relerr) | MUG L2 | MUG u_mean/u_c (relerr) |
|----|----|----|----|----|----|----|
| 1  | 0.6774 | 1.04e-3 | 1.08e-3 | 0.6783 (1.3e-3) | — | — |
| 2  | 0.7055 | — (no OF run) | — | — | 1.33e-4 | 0.7056 (7.2e-5) |
| 5  | 0.8109 | 1.95e-4 | 2.94e-4 | 0.8118 (1.1e-3) | 1.25e-4 | 0.8108 (1.2e-4) |
| 10 | 0.9001 | 6.90e-4 | 1.19e-3 | 0.9012 (1.3e-3) | 4.31e-4 | 0.9004 (3.7e-4) |
| 20 | 0.9500 | 5.79e-4 | 1.46e-3 | 0.9511 (1.2e-3) | 7.29e-5 | 0.9500 (2.0e-5) |
| 50 | 0.9800 | 5.97e-4 | 2.48e-3 | 0.9810 (1.1e-3) | — | — |

(Ha=2 has a MUG point but no mhdFoam run; Ha=5, Ha=10, and Ha=20 now have both mhdFoam and MUG
points. Ha=10 and Ha=20 are single packed-grid Ginsburg probes, not mesh-convergence studies.
The MUG `u_mean/u_c`
relerr is vs analytic at the integer Ha in the regenerated table; see §3 for the MUG-internal
Ha fit/predict check.)

**mhdFoam vs analytic:** L2 ≤ 1.0e-3, Linf ≤ 2.5e-3, `u_mean/u_c` relerr ~1.1e-3 across the
whole sweep. The discriminating scalar tracks analytic from near-Poiseuille (Ha=1) to plug-core
(Ha=50). This is the solid, guaranteed result.

---

## 3. Task A — resolved MUG leg (Ha mapping confirmed from first principles)

Driver: `src/tests/physics/test_hartmann_mug.F90`. Setup: x-periodic channel, no-slip walls,
applied field `B_0=(0,0,B0)` (wall-normal, in-plane), `n,T` frozen (incompressible), insulating
`ψ` BC, driven by `body_force=(fx,0,0)`; marched to steady state; `velx(y)` extracted.

| Ha (pred) | B0 | mesh (y cells) | Ha (fit) | fit/pred | profile L2 vs analytic | u_mean/u_c (MUG / analytic) |
|----|----|----|----|----|----|----|
| 2.002 | 2.244e-3 | 80  | 2.002 | 1.000 | 4.67e-5 | 0.7056 / 0.7056 |
| 3.002 | 3.366e-3 | 80  | 3.003 | 1.000 | 6.61e-5 | 0.7421 / 0.7421 |
| 4.995 | 5.60e-3  | 80  | 4.998 | 1.001 | 1.25e-4 | 0.8108 / 0.8109 |
| 10.000 | 1.121e-2 | 160 packed (`packing=8`) | 10.009 | 1.001 | 3.41e-4 | 0.9004 / 0.9002 |
| 20.000 | 2.242e-2 | 320 packed (`packing=8`) | 20.007 | 1.000 | 5.55e-5 | 0.9500 / 0.9500 |
| 50 | — | — | _not run/finalized_ | — | — | — |

The Hartmann number **fitted from the computed profile shape** matches the **predicted Ha from
inputs** to ≤0.1% at Ha=2, 3, and 5, the profiles match analytic to L2 ~1e-4, and `u_mean/u_c`
matches to ≤1e-4 for Ha≤5. The new Ginsburg Ha=10 packed-grid probe also matches analytic well:
developed-flow nonuniformity `max std(velx)/max|u| = 1.347e-7`, `Ha_fit/Ha_pred = 1.0009`,
profile L2 `3.406e-4` vs analytic at fitted Ha, and `u_mean/u_c = 0.9004` vs analytic `0.9002`.
The Ginsburg Ha=20 packed-grid probe completed successfully with developed-flow nonuniformity
`1.312e-3`, `Ha_fit/Ha_pred = 1.0004`, profile L2 `5.548e-5` vs analytic at fitted Ha, and
`u_mean/u_c = 0.9500` vs analytic `0.9500`. This confirms both the geometry mapping (x–z mesh
plane, y out-of-plane; induced `b_x = -dψ/dz`) and the scaling
`Ha = B0·a/√(μ0·η·ρ·ν)` through Ha=20. It does **not** yet prove high-Ha mesh convergence;
the attempted Ha=20 `NY=400` refinement did not reach solver initialization/output in this wrapper,
the attempted Ha=20 smaller-time-step repeat showed the same pre-initialization stall pattern, and
Ha=50 remains a meshing/solver-cost extension (§6).

---

## 4. Rungs: completed / blocked

| Rung | Status | Notes |
|------|--------|-------|
| A1 analytic Hartmann | ✅ complete | `analytic_hartmann.py` self-test passes (Poiseuille 2/3, plug→1, Q*~1/Ha, exact vs numeric <1e-8). |
| A3 mhdFoam Hartmann | ✅ complete | Ha∈{1,5,10,20,50}, agreement ~1e-3. |
| A2 resolved MUG Hartmann | ✅ complete through Ha=20 probe | Ha∈{2,3,5} plus packed-grid Ginsburg Ha=10 and Ha=20; Ha mapping confirmed to 0.1%. Ha=20 still needs convergence controls; Ha=50 not finalized (mesh/time-step cost, §6). |
| B analytic Hunt | see §7 | conducting-wall duct, side jets. |
| B MUG Hunt | not attempted | Hunt needs a *different* geometry (out-of-plane axial flow `vely(x,z)`, conducting Hartmann-wall EM BC). The Task-A MUG geometry (in-plane channel) does not transfer directly; flagged as next step rather than risk a long debug overnight. |

**Known issue worked around (logged):** the OpenFOAM v1912 *sampling functionObject* fails with
an `IOstream "sha1"` error (reproduced on local disk too, not a mount issue), so velocity profiles
were extracted via `foamToVTK` + `meshio` cell-centroid selection (`cases/hartmann/extract_profile.py`).
**Self-inflicted issue (logged):** editing a running `runTimeModifiable` case's `controlDict`
(live `endTime` change) raced with the solver's re-read and crashed the Ha=1 mhdFoam run at t=10;
its t=10 data was already written and converged (L2=1.0e-3), so it was used.

---

## 5. Code added to MUG (default OFF) and regression status

`src/physics/xmhd_2d.F90` — new **uniform body-force momentum source**, default off:
```fortran
LOGICAL  :: use_body_force = .FALSE.
REAL(r8) :: body_force(3)  = [0,0,0]   ! per-unit-mass acceleration
```
Added to the nonlinear residual only (constant ⇒ no Jacobian contribution), both Cartesian and
cylindrical paths, sign-consistent with the pressure-gradient term. Used to drive fully-developed
channel flow for this verification.

**Regression (drag and body-force both default off):**
`test_alfven_2d` (9), `test_sound_2d` (9), `test_drag_decay` (3) → **21 passed, 36 skipped** —
unchanged from before the body-force change. The new code paths are inert when disabled.

---

## 6. What remains for the high-Ha resolved-MUG leg

The Hartmann layer thickness is ~`a/Ha`. Resolving it with ≥8 cells needs wall-normal spacing
`≲ a/(8·Ha)`. The 2D MUG cube mesh is **uniform** (no wall-normal grading), so resolution must be
added everywhere unless packing/refinement is used. Ha=5 is comfortable on 80 uniform wall-normal
cells. A Ginsburg probe reached Ha=10 using `NY=160`, `packing=8.0`, and a smaller stable
time step `dt=5e-3`; the earlier `dt=2e-2` trial stalled before the first completed timestep and
was cancelled as a timestep-stiffness probe failure. A second Ginsburg probe reached Ha=20 using
`NY=320`, `packing=8.0`, `dt=1.25e-3`, and `NSTEPS=1600`, with excellent profile/Ha agreement but
measured x-nonuniformity `1.312e-3`. Ha=50 (layer 0.02) would need either about 400 uniform cells
across the channel or a stronger wall-normal packing/refinement strategy. Two bounded Ha=50
attempts were made after the Ha=20 success:

- Job `8910721`: `NY=500`, `packing=12.0`, `dt=2e-4`, `NSTEPS=5000`; failed before timestepping
  with OFT run status 99, `Ground point repeated on same domain` in `mesh_global_igrnd` during
  square mesh generation.
- Job `8910734`: `NY=500`, `packing=8.0`, same `dt` and `NSTEPS`; cancelled after 00:08:42 with
  metadata only, no mesh statistics or timestep output, and low CPU accumulation.

An Ha=20 refinement repeat was also attempted:

- Job `8910763`: `Ha=20`, `NY=400`, `packing=8.0`, `dt=1.25e-3`, `NSTEPS=1600`; cancelled on
  `g063` after 00:02:23 with metadata only, no mesh/timestep output, and low CPU accumulation.
- Job `8910766`: same case with `--exclude=g063`; cancelled on `g090` after 00:11:56 while
  CPU-active, but still with only the initial branch/revision lines and no OFT initialization,
  mesh statistics, timestep output, profile, or status file.
- Job `8911505`: `Ha=20`, `NY=320`, `packing=8.0`, `dt=6.25e-4`, `NSTEPS=3200`,
  `--exclude=g063,g090`; cancelled on `g198` after 00:05:44 with metadata only, initial
  branch/revision lines only, no OFT initialization, mesh statistics, timestep output, profile,
  or status file, and low CPU accumulation.

mhdFoam, by contrast, used a graded mesh (wall cell ~0.002) to reach Ha=50 cheaply. **Next step
for MUG:** add a true mesh/time-step convergence check for Ha=10/20 and extend the packed-grid
probe to Ha=50 with a mesh/setup strategy that avoids the failures above before claiming a
high-Ha resolved-MUG result. This is a meshing/solver-cost limitation of *this exercise's setup*,
not a statement about the solver's correctness.

---

## 7. Task B — Hunt conducting-wall duct (analytic/numerical reference)

Hunt's flow: rectangular duct, transverse field, **perfectly conducting Hartmann walls**
(⊥ B) + **insulating side walls** (∥ B). Signature physics: high-velocity **side jets** near
the insulating side walls and a near-uniform core — qualitatively different from the jet-free
insulating Hartmann channel of Task A.

**Method (honest):** rather than the closed-form double series (error-prone to transcribe
correctly from memory), I used a **finite-difference solve of Hunt's governing inductionless
MHD duct equations** as the referee (`cases/hunt/hunt_duct.py`):
```
lap(u) + Ha·∂b/∂y = -1 ,   lap(b) + Ha·∂u/∂y = 0
u=0 all walls ;  b=0 on insulating side walls ;  ∂b/∂y=0 on conducting Hartmann walls
```
This is the exact problem Hunt solved; the FD solution unambiguously reproduces the signature
physics. The overnight-run instruction explicitly permits verifying the well-known qualitative
features over a shaky series — logged.

**Result — side jets confirmed and growing with Ha** (`cases/hunt/results/hunt_sidejet_table.md`,
mid-plane y=0 peak/centerline ratio):

| Ha | peak/centerline (y=0) | x of peak | side jets? |
|----|----|----|----|
| 5  | 1.008 | −0.300 | incipient |
| 10 | 1.793 | −0.657 | yes |
| 20 | 6.225 | +0.771 | yes |
| 50 | 9.556 | +0.857 | strong |

The side-jet peak grows from ~1 (no jets, near-Hartmann) at Ha=5 to ~10× the centerline at Ha=50,
with the peak migrating toward the side walls (x→±1) as Ha rises — the textbook Hunt behavior, and
clearly distinct from the insulating Hartmann channel (no jets, flat core). Contour and mid-plane
plots: `cases/hunt/results/hunt_contour_Ha*.png`, `hunt_midplane_Ha*.png`.

**MUG Hunt run: not attempted.** Hunt's duct requires a *different* MUG geometry from the Task-A
channel — axial flow out-of-plane (`vely(x,z)` over the 2-D cross-section), walls on all four
sides, and a *conducting* Hartmann-wall EM BC (`ψ` Dirichlet vs the insulating natural BC used in
Task A). The in-plane channel mapping does not transfer directly. To avoid a long unattended debug
spiral, this is flagged as the top next step rather than attempted tonight (§8).

---

## 8. Next steps

1. Wall-normal-graded 2D mesh to add convergence controls at Ha=20 and extend the resolved-MUG Hartmann leg to Ha=50.
2. MUG Hunt: implement the out-of-plane axial-flow duct geometry (`vely(x,z)`) with a conducting
   Hartmann-wall EM BC (`ψ` Dirichlet) and compare to the analytic side-jet profile.
3. Fix or report the OpenFOAM v1912 sampling `sha1` bug upstream (workaround in place).
4. (Separate track) the reduced wall-drag closure remains validated only by `test_drag_decay`;
   calibrating its `drag_coeff` against these resolved Hartmann results is a natural follow-up.

---

_Artifacts: `cases/hartmann/` (analytic module, OpenFOAM template, extraction + compare + analyze
scripts, `results/`), `cases/hunt/`, running log `docs/dev/lm_mhd_verification_overnight_log.md`._

---

## 9. 2026-07-10 Paper 2 continuation: plasma/LM reduced coupling ladder

This continuation added a verified reduced plasma/liquid-metal interface ladder
for paper planning. It is **not** a nonlinear plasma/free-surface solver and it
does **not** validate two-way plasma/LM dynamics. It provides formula-level
boundary-load contracts, thermal/current/surface-response gates, and
reproducible artifacts for deciding which full simulations are worth attempting.

Completed rungs:

- L0 sheath/load contract: Bohm speed, particle/current/heat/mass/momentum
  loads, floating net-current default, and evaporative heat/mass accounting.
- L1 thermal response: representative Li/PbLi thermal scales, penetration
  depth, lumped-layer energy balance, and semi-infinite constant-flux response.
- L2 biased-current feedback: explicit non-floating current path and ohmic
  surface-current heating scale.
- L3 linear thermocapillary fallback: Marangoni shear and thin-layer Couette
  response.
- L4 gate only: linear capillary-gravity-magnetic stiffness, low-Rm magnetic
  damping, explicit timestep pressure, displacement fraction, and thermal
  penetration diagnostics before any nonlinear free-surface calculation.

Measured verification:

- Local focused LM-MHD suite after L4 gate: `67 passed in 0.70s`.
- Local artifact manifest: 47 files, 14 CSV files, 19 PNG files with 19/19
  nonblank.
- Local reference-matrix driver: all five plasma-interface diagnostics returned
  zero; focused suite recheck `67 passed in 0.86s`; manifest refreshed to 49
  files, 15 CSV files, 19 PNG files with 19/19 nonblank.
- Ginsburg L0-L3 Python/reference/artifact baseline job `8910539`: completed
  `0:0` on `g256`; `65 passed in 38.88s`; Hartmann analytic self-test passed;
  remote manifest refreshed to 44 files.
- Ginsburg L4 Python/reference/artifact baseline job `8910555`: completed
  `0:0` on `g178`; `67 passed in 2.08s`; Hartmann analytic self-test passed;
  all five plasma-interface scans ran; remote manifest refreshed to 47 files.
- Ginsburg reference-matrix wrapper baseline job `8910567`: completed `0:0`
  on `g174`; `67 passed in 2.40s`; Hartmann analytic self-test passed; all
  five reference-matrix commands returned zero; remote manifest refreshed to 49
  files.
- Ginsburg OFT build/regression baseline job `8902073`: completed `0:0` on
  `g099`; OFT build passed via source-dependency fallback; targeted physics
  regressions passed (`21 passed, 36 skipped`) in `2079.55s`.

New primary artifacts:

- `src/python/OpenFUSIONToolkit/LM_MHD/plasma_interface.py`
- `src/tests/physics/test_lm_mhd_plasma_interface.py`
- `cases/lm_mhd_coupling/plasma_interface_sheath_scan.py`
- `cases/lm_mhd_coupling/plasma_interface_thermal_response_scan.py`
- `cases/lm_mhd_coupling/plasma_interface_current_feedback_scan.py`
- `cases/lm_mhd_coupling/plasma_interface_marangoni_scan.py`
- `cases/lm_mhd_coupling/plasma_interface_free_surface_gate_scan.py`
- `cases/lm_mhd_coupling/run_reference_matrix.py`
- `docs/dev/lm_mhd_two_paper_execution_plan.md`
- `docs/dev/lm_mhd_paper_execution_log.md`

Next honest step: for Paper 1, extend the packed-grid resolved-MUG Hartmann
probe from Ha=20 to Ha=50 and add mesh/time-step convergence; for Paper 2,
keep the current result at L0-L3 plus an L4 gate unless a tightly scoped
nonlinear surface demonstration passes its conservation/stress-balance tests.
