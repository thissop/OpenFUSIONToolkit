# Hartmann wall-drag closure — next steps (staged for Ginsburg/COMSOL)

Context: Guizzo et al. 2026 name effective Hartmann wall-drag/wall-function models as their #1
future-work item. We already had the `use_wall_drag` term in `xmhd_2d.F90` but it used the
CHANNEL closure (over-brakes a duct by ~Ha). **DONE:** added a `drag_closure` selector
(raw|channel|duct_insulating|duct_conducting) computing the correct `alpha_eff=(nu/a^2)f(Ha,c)`
(OFT lm-mhd-upgrades `bf7d446`), and locally validated the insulating closure to machine precision
(resolved duct core u0*Ha=1.0000; figure `F_dragclosure`). Remaining work below is HPC-gated.

## 1. Wall-function validation run (Ginsburg) — coarse mesh + drag = resolved
**Goal:** prove the drag term reproduces the resolved core velocity WITHOUT resolving the Hartmann
layer (that is what "wall-function" means, and what Guizzo needs).
- **Driver:** `test_shercliff_mug` (insulating) and `test_hunt_mug` (conducting), with a drag-on
  path. NOTE: both currently hard-code `use_wall_drag=.FALSE.`. Add a namelist/driver switch:
  `use_wall_drag=T`, `drag_closure='duct_insulating'` (or `'duct_conducting'`, `drag_c=c`),
  `drag_Ha=Ha`, `drag_bhat=[1,0,0]` (perp = cross-field), `nu`, `a_half` set as usual.
- **Mesh:** COARSE, Hartmann layer UNRESOLVED (e.g. ni=16-32, packing=1, NO boundary-layer grading)
  — the whole point is no Hartmann resolution.
- **Check:** core velocity u0 should hit `1/Ha` (insulating) / `(1+1/c)/Ha^2` (conducting) to a few
  %, matching the RESOLVED reference. Control: same coarse mesh with `drag_closure='channel'` should
  be wrong by ~Ha; with drag OFF, unresolved/garbage.
- **Sweep:** Ha in {20,100,1000}, insulating + conducting c=1. ~minutes each (small 2D).

## 2. Apply to Guizzo's device model (Ginsburg) — the "we improve her result" demo
- **Case:** her PR #320 `src/examples/TokaMaker/AdvancedWorkflows/Blanket/vde_blanket.ipynb`
  (coupled `mugtok_td`, Pb-Li blanket VDE). Her headline blanket velocities (1.4-3.5 m/s) are
  computed with NO Hartmann drag.
- **Do:** turn on `use_wall_drag=T` + `drag_closure='duct_conducting'` (her blanket is toroidally
  conducting) with her `Ha`, `c`, `nu`, `a`. `drag_bhat` = the poloidal-field direction.
- **Quantify:** how far the corrected Hartmann drag pulls her peak velocity down from 1.4-3.5 m/s.
  This is the concrete two-paper-compose deliverable. (Her PR #320 is OPEN/unmerged — coordinate with
  Hansen; our contribution touches the same `xmhd_2d.F90`.)

## 3. COMSOL cross-checks (remote desktop, later)
- Optional: a COMSOL coarse-mesh-with-effective-drag vs resolved comparison as an independent second
  code for the wall-function (mirrors the duct two-coding). Low priority — the resolved-vs-closure
  validation is already analytic-exact for insulating.

## Runner-up gaps (Guizzo future-work we don't yet have)
### R1. Open / non-reflecting outflow BC — needed for open/separated geometries
- Paper-2 flagged this for sudden-expansion / backward-facing-step manifolds. Current duct BCs are
  no-slip on ALL walls (`setup_bc` defaults every velocity mask to `oft_blagrange%global%gbe`).
- **KEY FINDING (scoped this session):** the first-cut open outflow needs NO change to `xmhd_2d.F90`.
  `setup_bc` only defaults an UNassociated mask to `gbe`; a driver that allocates its own
  `velx_bc/vely_bc/velz_bc` and sets them `.FALSE.` on the outlet nodes (`.TRUE.` elsewhere) leaves
  the outlet velocity FREE — the momentum weak form's natural condition (zero viscous traction,
  `d(u)/dn=0`) then applies there. This is the standard "do-nothing" outflow and is adequate for
  steady / mildly-unsteady separated flow. Same override pattern the hunt driver already uses for
  `by_bc`. => R1 first cut is a NEW DRIVER + mesh, not shared-physics surgery (cannot break flagship).
- **Plan:**
  1. New driver `test_bfs_mug` (backward-facing step) or reuse a rectangular mesh with an inlet
     Poiseuille profile (Dirichlet inlet), no-slip walls, and a FREE outlet mask. Drive with a fixed
     inlet, not a body force.
  2. Validate: primary-bubble reattachment length x1/S vs Re against Armaly et al. 1983 (expansion
     ratio 1.94). Use the 2D-clean range only: Re=100 -> x1/S ~= 3.0, Re=200 -> ~5.4, Re=389 -> ~8.
     (Above Re~400 the flow is 3D and 2D under-predicts x1 -- do not validate a 2D code there.)
     Confirm the outlet passes the exiting shear layer with no spurious upstream reflection (steady,
     monotone convergence to the benchmark).
- **Refinement (only if the natural BC reflects):** an explicit convective/Sommerfeld update on the
  outlet nodes, `u^{n+1}|_out = u^n|_out - dt * U_n * (du/dn)^n` (advect the boundary value out; `U_n`
  = mean outlet normal speed). Implementable as a per-step Dirichlet value set from the adjacent
  interior gradient — still no weak-form surgery. Add behind an `outflow_convective` flag, default off.
- Scope: driver + meshing + validation (local authorship; Ginsburg compile/run). NO physics-module
  edit for the first cut — this is the safe path.

### R2. Segmented / channelized blanket geometry
- Guizzo: "practical blankets will be segmented into multiple poloidal channels separated by
  structural walls." Her axisymmetric model is a continuous tank.
- **Plan:** this is a MESHING + multi-region effort (poloidal channels + inter-channel walls as
  conjugate regions). Her multi-domain psi/F framework can carry it; the new physics we add is the
  toroidal-Hartmann-wall DRAG on each channel (our closure, item 1, applied per channel). Larger;
  scope only for now.

### R3. Adopt our thin-wall Robin BC for the induced field in her device model
- She uses domain-resolved conductors; our thin-wall Robin BC (`c_wall`) gives the induced-field wall
  coupling without resolving the wall — cheaper. Drop-in for her conducting-blanket boundaries.
- Scope: small (she can enable `c_wall` on the blanket boundary); validate vs her resolved-wall case.

## Priority order when HPC returns
1 (wall-function validation, closes our loop) -> 2 (Guizzo demo, the headline) -> R1 (outflow BC,
the next real capability) -> R3 (Robin adoption) -> R2 (segmented, largest).
