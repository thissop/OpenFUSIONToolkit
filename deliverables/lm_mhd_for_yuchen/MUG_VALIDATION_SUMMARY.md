# MUG LM-MHD Validation Summary (Paper 2)

Consolidated results of the MUG (`xmhd_2d`) validation campaign against the COMSOL
full-induction reference and analytic benchmarks. Case contract:
`lm-mhd-comsol/VALIDATION_CASES.md`. Cross-agent record: `docs/dev/oft_to_comsol.md`
(MUG→COMSOL) and `lm-mhd-comsol/SYNC/comsol_to_oft.md` (COMSOL→MUG).

All MUG runs: `xmhd_2d` 2.5D visco-resistive MHD, drag OFF, native full induction (psi, by).
Drivers in `src/tests/physics/`, scorers + profiles in `cases/shercliff/`. Nondim per contract:
velocity in unit-forcing units, profiles vs y/a ∈ [−1,1], Hartmann layers ~1/Ha.

## Closed-form vs non-closed-form (important scoping distinction)

Of the closed rungs, **only V4a is genuinely non-closed-form.** Shercliff (S1/S1a/S1-1323) and
Hunt (H1 perfect *and* finite-c) all have analytic series (Shercliff 1953, Hunt 1965), so those
are **solver-correctness** checks — necessary and rigorous, but not validation against physics
with no closed form. V4a (trapezoid) is non-analytic only via geometry non-separability. A
*fully-developed* finite-Rm conjugate duct (old "V1") is also secretly linear/analytic: in the
reduced out-of-plane-flow/in-plane-B₀ model the (u, by) system is linear and Rm only rescales by,
so Hunt is exact at any Rm. Genuine finite-Rm nonlinearity requires **cross-stream flow**
(magnetic obstacle / sudden expansion) or **transient** conjugate coupling. This is why the
load-bearing non-analytic case is V4 (expansion) or V2 (transient conjugate), not V1.

## Status board

| Rung | Case | Reference | MUG result | Agreement | Status |
|---|---|---|---|---|---|
| Analytic 1 | Shercliff duct, Ha=20 | series (double-sourced) | U0 core | 0.004% core, relL2(u) 5.6e-4 | ✅ |
| Analytic 2 | Shercliff S1a, Ha=100 | series + COMSOL | U0 9.99999956e-3 | 4e-9 core, relL2 1.5e-4 | ✅ |
| Analytic 3 | Shercliff S1, **Ha=1323** (canonical) | series (2.5e-10 dual) + COMSOL | U0 7.55857866e-4 | **3e-9 core**, relL2 4.3e-5 | ✅ |
| Non-analytic | V4a trapezoid duct, Ha=100 | COMSOL (mesh-conv 6e-9) | field | **relL2(u) 7.8e-5**, relL2(b) 1.3e-4, Q 0.006% | ✅ |
| Conjugate | H1 Hunt perfect-conductor, Ha=20 | FD referee + COMSOL | u0 1.8205e-3 | core ~1–2%, **max\|b\| 0.13%**, jets 1.7% | ✅ |
| Phase-4 | H1 Hunt finite-c **c=1**, Ha=20 | COMSOL Tier3 conjugate | u0 4.2188e-3 | **core 0.005%**, field relL2 1.2% | ✅ |
| Phase-4 limit | Robin BC c→0 (insulating) | Shercliff series | 4.9831e-2 | 0.18% | ✅ |
| Phase-4 limit | Robin BC c→∞ (perfect) | H1-perfect same mesh | 1.8081e-3 | 0.68% | ✅ |
| **★ Transient capstone** | **V2 transient conjugate, Ha=100, c=1, Pm=1** | COMSOL thin-wall transient ref | v(t,y), by(t,y) @ x=0 | **relL2_v 0.049, relL2_by 0.027, Alfvén overshoot timing err 0.003 τ_A**, by/v amp ratio 3.8% | ✅ |
| Developing | V3 entrance flow, Ha=100 | COMSOL (L_e=0.301) | — | dropped (non-capstone; ~50h grind, both agents concur) | ⏸ |
| Conjugate hi-Ha | H1 finite-c, Ha=1323 | — | — | deprioritized (polish; both agents concur) | ⏸ |

**★ V2 is the load-bearing non-analytic capstone (genuinely finite-Rm, transient conjugate).** MUG
Phase-4 thin-wall vs COMSOL thin-wall reference, rest IC + step-on `−∂p/∂x`, observed 4·τ_A. Match
holds through 2 full Alfvén oscillation periods (period 2τ_A, phase cycle-for-cycle). The damping
log-decrement agrees three ways: MUG 0.10 ≡ COMSOL thin-wall 0.105 ≡ analytic (ν+η)k²/2·2τ_A ≈ 0.10.
**Physics result:** the thin-wall Robin BC's *steady* conductance-equivalence (Hunt c=1 to 0.005%,
above) does **not** extend to the *transient* — a finite-thickness wall (tw=0.2) adds magnetic-diffusion
dissipation (τ_wall≈20 ≫ T_osc=4) that raised the thick-wall relL2 to 0.20; the thin-wall limit (matching
MUG's BC) collapses it to <5%. Both wall limits now cross-validated MUG↔COMSOL. Tooling:
`cases/shercliff/v2_waveform.py` (extract) + `v2_compare.py` (relL2 + overshoot-timing metric).

## Code changes made to MUG during validation (candidate upstream PR)

Two latent bugs in the `by` induction equation's field-line stretching source, both invisible to
every prior use of the code (no test exercised the duct polarization), found by the validation
ladder:

1. **Missing background-field stretch term** `(B_0·∇)v_y` (`xmhd_2d.F90`, `nlfun_apply` +
   Jacobian). Identically zero for the code's original out-of-plane-B_0 use cases, but the *sole*
   driver of the LM-MHD duct ladder (out-of-plane flow, in-plane B_0). Without it every duct case
   silently degenerates to hydrodynamics. Verified against Shercliff (0.004%).
2. **Sign error in the pre-existing psi-based stretch term.** Unexercised by any prior test
   (`test_alfven_2d` covers only the velx↔psi polarization). New `test_alfven_by2d` (by↔vely
   shear-Alfvén): old sign grows exponentially at k·v_A; corrected sign oscillates at the analytic
   period (0.1004 vs 0.1000). Both fixes regression-clean.

**Phase-4 thin-wall/conjugate Robin BC** (`by + c_wall·∂by/∂n = 0`): implemented in residual and
Jacobian via the exact order-1 boundary-edge mass matrix, driver-controlled per wall through the
`by_bc` mask; order/cylindrical guards in `setup_bc`. Was a default-off stub at the start of the
collaboration; now validated across the full c-range (table above).

## Solver notes (for reproducibility)

- **Adaptive dt is capped at the initial namelist dt and only shrinks** (`xmhd_2d.F90:366`,
  `dt = ittarget·avg(dt)/avg(lits)`, then `min(dt, dtin)`). There is no ramp-up: set the namelist
  dt to the *largest step-1-survivable* value. Small "startup" dt stays small forever.
- High-Ha conjugate-wall cases settle much slower than insulating (Hunt ~4× Shercliff at Ha=20);
  high-Ha finite-c from rest needs restart-chaining — the reason Ha=1323 finite-c is deprioritized.
- High-Ha meshes: built-in generator `packing=px,py,pz` gives per-direction boundary-layer grading
  (S1-1323 used 256²/packing 32, ~3.4 cells/Hartmann layer).

## Headline for Paper 2

An independent COMSOL/analytic validation ladder established MUG as a credible LM-MHD full-induction
solver — matching to 10⁻⁹ on analytic anchors and 10⁻⁴–10⁻² on non-analytic and conjugate-wall
references — and in the process found and fixed two latent full-induction bugs and delivered the
conjugate-wall (Phase-4) capability that reaches the finite-c corner where quasi-static surrogates
fail (Paper 1's regime).
