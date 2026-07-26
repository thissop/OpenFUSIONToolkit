# MUG 3D LM-MHD — Scoped Plan (2026-07-26)

**Goal:** extend MUG (OpenFUSIONToolkit) from the validated 2.5D LM-MHD solver (`xmhd_2d`) to a **3D**
LM-MHD capability using the 3D solver `xmhd.F90`, culminating in a **G1-scale high-Ha 3D run** on Ginsburg
that hedges against a COMSOL-3D effort. **Hard robustness gate:** no expensive high-Ha run until a
mesh-convergence study proves the Hartmann layer (δ_Ha = a/Ha) is resolved at the target Ha.

## The core challenge (why this is staged)
The Hartmann layer is δ_Ha = a/Ha thick. **Volume-integrated loads** average over it and are forgiving
(validated to a few % at Ha=2646 with only ~1.6 cells). **Pointwise near-wall quantities** — velocity
peak, and the *timing/overshoot* gates — need ~10 cells across δ_Ha and are the real constraint.

**Honest current 2.5D pointwise ceiling** (from real mesh data, packed-cube mesh):
- 256²/pack48 (hmin=2.33e-4): 10-cell converged to **~Ha 428**; 384²/pack64 (hmin=1.40e-4): **~Ha 714**.
- Reaching disruption-Ha (2646) pointwise needs **~6× finer near-wall** → requires *grading*, not uniform
  refinement (which OOMs). Levers, in order: **(1) geometric boundary-layer grading**, **(2) higher-order
  elements** (blocked at order=1 in 2.5D conjugate-wall, but native in 3D H(curl)), **(3) anisotropic
  wall-normal refinement**, (4) AMR.

## Key scoping fact
`xmhd.F90` is OFT's **plasma** MHD solver: mixed FEM (B in H(curl)+H(grad), v/n/T in Lagrange), and its
only tests are Alfvén/sound waves in a periodic box. **There is no duct, Hartmann, forcing, or
conjugate-wall machinery** — my 2.5D scalar-flux LM-MHD code does not port directly. So this is a real
build, formulation-by-formulation.

---

## STAGE 0 — Foundation (understand the 3D solver)
- 0a. Read `test_alfven.F90` / `test_sound.F90` fully as the canonical 3D-MHD driver template (setup of
  the FE spaces, time loop, forcing, output).
- 0b. Map `xmhd.F90`'s interface: the `oft_xmhd_driver` forcing class, the BC flags (b1_bc/b2_bc/v_bc/
  n_bc/t_bc), how a static/applied field and a sustained drive are (or are not) supported.
- 0c. Determine the 3D mesh path: can we build a graded duct box (streamwise × spanwise × wall-normal)?
  What grading does the mesher support?
- **Exit:** a written note on how to stand up a driven, walled 3D MHD case — the prerequisites for Stage 1.

## STAGE 1 — Low-Ha 3D Hartmann duct feasibility spike  *(de-risk the path, cheap)*
- 1a. 3D duct mesh: a box with insulating side walls + Hartmann walls; modest resolution (low Ha ⇒ thick
  layer, easy to resolve).
- 1b. Axial drive: a uniform body force / pressure gradient along the duct (add to `xmhd.F90` if absent,
  mirroring the `body_force` I added to `xmhd_2d`), plus a transverse uniform B0.
- 1c. Insulating no-slip BCs (v=0 all walls; insulating B condition in the H(curl) formulation).
- 1d. Run at **Ha≈10–20** → extract the axial-velocity profile → **validate against the analytic Hartmann
  solution** (core flattening + Hartmann-layer shape) to a few %.
- **Exit:** MUG's 3D solver reproduces a Hartmann duct at low Ha ⇒ the 3D LM-MHD path is real.

## STAGE 2 — Step Ha up + mesh-convergence study  *(THE robustness gate)*
- 2a. Add the **conjugate (conducting) wall BC** in the 3D H(curl)/H(grad) B-formulation (the LM-MHD c=1
  case). This is the hard formulation work; my 2.5D thin-wall Robin was on the scalar `by`.
- 2b. Sweep Ha (100 → 300 → 500 → 1000 → …), at each running a convergence study with **geometric BL
  grading + higher-order** elements; track loads AND pointwise quantities.
- 2c. **MEASURE the 3D pointwise Ha ceiling** — the max Ha where δ_Ha is resolved and the timing/overshoot
  gates converge. This number IS the go/no-go for Stage 3.
- **Exit:** a measured, documented Ha ceiling; confidence the target Ha is reachable (or an honest cap).

## STAGE 3 — G1-scale high-Ha 3D run  *(only after Stage 2 clears the gate)*
- 3a. Estimate mesh size + Ginsburg cost (multi-node MPI) at the target Ha.
- 3b. Run the G1-scale case; validate against COMSOL-3D / 2.5D / analytic.
- **Do NOT start until Stage 2 proves resolution at the target Ha.**

---

## Parallel / non-blocking
- **2.5D cross-validation still closing** (P1 anchor job 9194431, ILU+pre_freq=8) — independent of 3D.
- **Paper-2** gains a 3D-capability section if Stages 1–2 succeed.

## Tracking
Tasks #1 (close 2.5D), #2 (Stage 1), #3 (Stage 2), #4 (Stage 3). This plan lives at
`OpenFUSIONToolkit/docs/dev/MUG_3D_LM_MHD_PLAN.md`.
