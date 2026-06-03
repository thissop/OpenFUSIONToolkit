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

## 1. One-paragraph summary

Three independent solutions of the **insulating-wall Hartmann channel** were computed and
cross-checked: (1) the analytic profile, (2) the resolved 2D MUG solver (`xmhd_2d.F90`, drag
off, driven by a new default-off body-force term), and (3) OpenFOAM `mhdFoam`. **All three
agree to ~0.1%** on the normalized velocity profile and on the discriminating scalar
`u_mean/u_c` across the Hartmann-number range. For MUG, the Hartmann number *fitted from the
computed profile shape* matches the Hartmann number *predicted from the input parameters*
(`Ha = B0·a/√(μ0·η·ρ·ν)`) to **0.1%**, independently confirming both the geometry mapping and
the parameter scaling. The analytic + mhdFoam legs cover `Ha ∈ {1,5,10,20,50}`; the resolved-MUG
leg is confirmed at `Ha ∈ {2,5}` (an attempted `Ha=10` run on a finer 160-cell uniform mesh did
not reach steady state within the run window and was stopped — higher Ha needs wall-normal mesh
grading the uniform 2D mesh does not provide; see §6). Task B (Hunt conducting-wall duct) is
addressed with a finite-difference reference; see §7.

---

## 2. Task A — three-way error table (THE HEADLINE)

Insulating Hartmann channel, profiles normalized by centerline velocity. `u_mean/u_c` is the
discriminating scalar (→2/3 Poiseuille at Ha→0, →1 plug core at Ha→∞). Full machine-readable
table: `cases/hartmann/results/error_table.md`. Overlay plots: `cases/hartmann/results/overlay_Ha*.png`.

| Ha | u_mean/u_c (analytic) | mhdFoam L2 | mhdFoam Linf | mhdFoam u_mean/u_c (relerr) | MUG L2 | MUG u_mean/u_c (relerr) |
|----|----|----|----|----|----|----|
| 1  | 0.6774 | 1.04e-3 | 1.08e-3 | 0.6783 (1.3e-3) | — | — |
| 2  | 0.7355 | — (no OF run) | — | — | 4.67e-5 | 0.7056 (≈0) |
| 5  | 0.8109 | 1.95e-4 | 2.94e-4 | 0.8118 (1.1e-3) | 1.25e-4 | 0.8108 (1.2e-4) |
| 10 | 0.9001 | 6.90e-4 | 1.19e-3 | 0.9012 (1.3e-3) | — (not finalized) | — |
| 20 | 0.9500 | 5.79e-4 | 1.46e-3 | 0.9511 (1.2e-3) | — | — |
| 50 | 0.9800 | 5.97e-4 | 2.48e-3 | 0.9810 (1.1e-3) | — | — |

(Ha=2 has a MUG point but no mhdFoam run; Ha=5 is the full three-way overlap — see
`cases/hartmann/results/overlay_Ha5.png`. The MUG `u_mean/u_c` relerr is vs analytic at the
fitted Ha; see §3 for the MUG-internal Ha fit/predict check.)

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
| 4.995 | 5.60e-3  | 80  | 4.998 | 1.001 | 1.25e-4 | 0.8108 / 0.8109 |
| 9.99  | 1.122e-2 | 160 | _not finalized_ | — | — | — |

The Hartmann number **fitted from the computed profile shape** matches the **predicted Ha from
inputs** to ≤0.1% at both Ha=2 and Ha=5, the profiles match analytic to L2 ~1e-4, and `u_mean/u_c`
matches to ≤1e-4. The fully-developed assumption is verified: `velx` is uniform in x to <1e-6.
This confirms both the geometry mapping (x–z mesh plane, y out-of-plane; induced `b_x = -dψ/dz`)
and the scaling `Ha = B0·a/√(μ0·η·ρ·ν)`. **The Ha=10 run** on a finer 160-cell uniform mesh did
not reach steady state within the run window (the implicit solve stiffens markedly on the finer
uniform mesh) and was stopped — see §6. The two confirmed points already demonstrate the profile
match and the linear `Ha`–`B0` scaling; extending to higher Ha is a meshing task (§6), not a
solver question.

---

## 4. Rungs: completed / blocked

| Rung | Status | Notes |
|------|--------|-------|
| A1 analytic Hartmann | ✅ complete | `analytic_hartmann.py` self-test passes (Poiseuille 2/3, plug→1, Q*~1/Ha, exact vs numeric <1e-8). |
| A3 mhdFoam Hartmann | ✅ complete | Ha∈{1,5,10,20,50}, agreement ~1e-3. |
| A2 resolved MUG Hartmann | ✅ complete (Ha≤5) | Ha∈{2,5}; Ha mapping confirmed to 0.1%. Ha=10 attempted but not finalized (slow implicit solve on finer uniform mesh, §6). |
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

## 6. Why the resolved-MUG leg stops at Ha≈5 here

The Hartmann layer thickness is ~`a/Ha`. Resolving it with ≥8 cells needs wall-normal spacing
`≲ a/(8·Ha)`. The 2D MUG cube mesh is **uniform** (no wall-normal grading), so resolution must be
added everywhere: Ha=5 is comfortable on 80 cells, but Ha=10 already pushed to 160 uniform cells,
and the **implicit Newton–Krylov solve stiffens sharply** on the finer uniform mesh — the Ha=10
run did not reach steady state within the run window and was stopped (no faked result). Ha=50
(layer 0.02) would need ~400 uniform cells across the channel. mhdFoam, by contrast, used a graded
mesh (wall cell ~0.002) to reach Ha=50 cheaply. **Next step for MUG:** a wall-normal-graded or
boundary-layer-refined 2D mesh (or anisotropic refinement) to push the resolved-MUG leg to
Ha=10–50 without a uniform-mesh blow-up in DOF count and solver cost. This is a meshing/solver-cost
limitation of *this exercise's setup*, not a statement about the solver's correctness — the two
confirmed points (Ha=2, 5) already match analytic to ~1e-4 with the Ha mapping verified to 0.1%.

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

1. Wall-normal-graded 2D mesh to extend the resolved-MUG Hartmann leg to Ha=20–50.
2. MUG Hunt: implement the out-of-plane axial-flow duct geometry (`vely(x,z)`) with a conducting
   Hartmann-wall EM BC (`ψ` Dirichlet) and compare to the analytic side-jet profile.
3. Fix or report the OpenFOAM v1912 sampling `sha1` bug upstream (workaround in place).
4. (Separate track) the reduced wall-drag closure remains validated only by `test_drag_decay`;
   calibrating its `drag_coeff` against these resolved Hartmann results is a natural follow-up.

---

_Artifacts: `cases/hartmann/` (analytic module, OpenFOAM template, extraction + compare + analyze
scripts, `results/`), `cases/hunt/`, running log `docs/dev/lm_mhd_verification_overnight_log.md`._
