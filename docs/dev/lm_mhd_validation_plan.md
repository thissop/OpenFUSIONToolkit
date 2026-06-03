# LM-MHD Validation Plan

**Branch:** `lm-mhd-upgrades`  
**Date:** 2026-06-02

---

## Implemented and Tested (Phase 1)

### Test 1: Drag decay (`test_drag_decay`) — ✅ PASSING (verified 2026-06-02)

**Purpose:** Verify that the reduced Hartmann/wall drag source term decays velocity
at the correct backward-Euler rate.

**Setup:**
- Domain: 2D periodic mesh, 4 cells in x, `ref_per=T,T,F`
- Fields: uniform velx=1.0, vely=0.75, n=n0, T=t0
- Physics: near-zero η, ν, χ, D → drag term dominates
- drag_bhat = [0,1,0] (y/toroidal direction)
- Natural (free-slip) velocity BCs; Dirichlet-fixed n and T
- Timestep pinned via large `ittarget` so the adaptive controller holds dt = dt_initial
  (required for the constant-dt analytic prediction to hold at any FE order)

**Expected:**
- velx decays: `velx(n) = velx0 / (1 + α·dt)^n`  
- vely is unchanged: `|vely_final - vely0| / vely0 < 1e-6`

**Tolerances:**
- velx relative error vs backward-Euler: < 2%
- vely relative error: < 1e-6

**Measured results (Ubuntu 25.10 aarch64, gfortran 15.2.0, MPICH 4.2.3, native LA backend):**

| Case | FE order | α | dt | nsteps | velx rel. err (⊥, damped) | vely rel. err (∥, undamped) | Result |
|------|----------|---|----|--------|---------------------------|------------------------------|--------|
| `test_drag_p2_r1` | 2 | 2.0 | 0.05 | 10 | 3.52e-5 | 4.44e-16 | PASS |
| `test_drag_p2_r1_high_alpha` | 2 | 5.0 | 0.05 | 10 | 3.42e-5 | 4.44e-16 | PASS |
| `test_drag_p3_r1` | 3 | 2.0 | 0.05 | 10 | 5.17e-5 | 5.33e-13 | PASS |

The perpendicular component decays at the backward-Euler rate to ~5e-5 (limited by the
nonlinear solve tolerance, not the drag model); the parallel component is undamped to
machine precision, confirming the perpendicular projection `u - (u·b̂)b̂` is correct.

**Files:**
- `src/tests/physics/test_drag_decay.F90`
- `src/tests/physics/test_drag_decay.py`

### Test 2 (regression): existing 2D MUG tests with drag default-off — ✅ PASSING

`test_alfven_2d` (9 passed) and `test_sound_2d` (9 passed) were run against the
drag-enabled build. Because `use_wall_drag=.FALSE.` is the default and these tests do
not set it, results are unchanged — confirming the new code path is inert when disabled.
Full suite: **21 passed, 36 skipped** (skips are MPI/coverage parametrizations).

---

## Planned but Not Yet Implemented

### Test 2: No-drag regression

**Purpose:** Verify that setting `use_wall_drag=.FALSE.` gives bitwise-identical
results to the original code path.

**How:** Run test_alfven_2d with default parameters (no drag) and confirm expected
tolerances from the original test still hold.

**Status:** Implicitly covered by existing `test_alfven_2d` and `test_sound_2d` tests
since `use_wall_drag=.FALSE.` is the default and the existing tests do not set it.

### Test 3: Anisotropic conductivity manufactured solution (Phase 2)

**Purpose:** Validate tensor resistivity machinery when implemented.

**Setup:** Eigenmode of ψ diffusion equation:
```
ψ(x,z,t) = sin(k_x x) sin(k_z z) exp[-(η_x k_x² + η_z k_z²) t]
```
For isotropic η: rate = η (k_x² + k_z²)  
For anisotropic η_x ≠ η_z: rate = η_x k_x² + η_z k_z²

**How:** Run 3D xmhd_2d in Cartesian mode with a single spatial Fourier mode, compare
decay rate to analytic formula.

**Status:** Not yet implemented (Phase 2).

### Test 4: Hartmann channel flow benchmark (Phase 3)

**Purpose:** Validate the incompressible LM-MHD solver (when implemented) against the
classical Hartmann flow solution.

**Analytic solution (1D Hartmann flow between insulating plates, distance 2a, B=B₀ ẑ, pressure gradient f_x):**
```
u_x(y) = (f_x a²)/(ν Ha/tanh(Ha)) * [1 - cosh(Ha y/a)/cosh(Ha)]
Ha = B₀ a sqrt(σ/(ρ ν))
Q / Q₀ = tanh(Ha)/Ha   (flow rate ratio vs no-field)
```

**How:** Set up narrow 1D-equivalent channel domain, impose pressure-gradient body force,
run to steady state, compare u_x profile and total flow rate.

**Status:** Not yet implemented. Requires incompressible solver (Phase 3).

### Test 5: Shercliff / side-layer benchmark (Phase 4)

**Purpose:** Validate current-closure wall models and side-layer physics.

**Status:** Not yet implemented. Requires conducting-wall BC (Phase 4).

---

## Canonical Benchmark References

1. **Hartmann (1937)** — original MHD channel flow solution
2. **Shercliff (1953)** — square duct with conducting walls
3. **Hunt (1965)** — rectangular duct with conducting/insulating walls
4. **Smolentsev et al. (2010)** — MHD pressure drop benchmark for fusion blankets
   (Fusion Engineering and Design, 85(7-9), 1111–1118)
5. **Bühler & Mistrangelo (2004)** — 3D effects in blanket ducts

---

## Physics Regime Guide

| Solver | Valid regime | Key limitations |
|--------|-------------|-----------------|
| `xmhd_2d` (compressible, Phase 1) | Rm ~ 0.1–10, Ma ~ 0.01–0.3, LM disruption transients | Not strictly incompressible; Hartmann drag is reduced model |
| `xmhd_2d` + drag (Phase 1) | Qualitative high-Ha duct flow, disruption studies | Drag coefficient α must be calibrated to Ha and geometry |
| Future incompressible solver (Phase 3) | Full LM blanket MHD, Ha ~ 100–10000 | Needs saddle-point solver; anisotropic conductivity |
| Inductionless low-Rm solver (Phase 3b) | Li/PbLi blanket flow, Rm ≪ 1 | B must be imposed externally (from TokaMaker) |

---

## How to Run Tests

Once built (requires CMake build with `OFT_BUILD_TESTS=ON`):

```bash
cd build/tests/physics
pytest test_drag_decay.py -v
pytest test_alfven_2d.py -v   # regression: drag=off by default
pytest test_sound_2d.py -v    # regression: drag=off by default
```
