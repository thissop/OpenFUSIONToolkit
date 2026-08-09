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
### R1. Non-reflecting / convective outflow BC (xmhd_2d) — needed for open/separated geometries
- Paper-2 flagged this as the missing piece for sudden-expansion / backward-facing-step manifolds.
  Guizzo also implies open-geometry extensions. Current BCs are Dirichlet/no-slip; an open outlet
  reflects.
- **Plan:** implement a convective (Sommerfeld-type) outflow BC `dt(phi)+U_n dn(phi)=0` on a tagged
  outlet boundary for velocity (and a zero-gradient / do-nothing for pressure), selectable per
  boundary like the Robin BC. Validate on a channel with an open outlet (steady Poiseuille passes
  through cleanly) then the backward-facing step (recirculation length vs Re, Armaly benchmark).
- Scope: moderate (a boundary-integral term + a tagged-outlet mask). Local code; Ginsburg validate.

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
