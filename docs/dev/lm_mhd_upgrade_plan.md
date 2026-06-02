# LM-MHD Upgrade Plan for `xmhd_2d.F90`

**Branch:** `lm-mhd-upgrades`  
**Date:** 2026-06-02  
**Author:** Repository-grounded analysis (Claude Code)

---

## A. Current Architecture Summary

`xmhd_2d.F90` is a self-contained 2D/2.5D visco-resistive **compressible** MHD module
using a scalar Lagrange FE basis (`oft_scalar_bfem`) on 2D surface meshes (`oft_bmesh`).
It is entirely independent of the 3D `xmhd.F90` / `xmhd_lag.F90` modules — different
mesh types, FE spaces, and solver structures.

The main simulation state lives in `oft_xmhd_2d_sim` (TYPE, lines 69–124). The nonlinear
residual is assembled in `nlfun_apply` (lines 717–1006). The approximate Jacobian is
assembled in `build_approx_jacobian` (lines 1010–1516). Both routines loop over cells,
then quadrature points, then test functions in the standard Galerkin FE pattern.

The module exposes `run_simulation` (nonlinear backward-Euler NK loop) and
`run_lin_simulation` (linearized advance around a fixed equilibrium). Build registration
is in `src/physics/CMakeLists.txt` (line 7). Tests use `oft_add_test()` in
`src/tests/physics/CMakeLists.txt` and are driven by pytest through `run_OFT()`.

---

## B. Exact PDEs Currently Solved by `xmhd_2d.F90`

Seven coupled scalar Lagrange fields:

| Field index | Tag    | Physical quantity           |
|-------------|--------|----------------------------|
| 1           | `n`    | number density              |
| 2           | `velx` | velocity component 1 (x/R) |
| 3           | `vely` | velocity component 2 (y/φ) |
| 4           | `velz` | velocity component 3 (z/Z) |
| 5           | `T`    | temperature                 |
| 6           | `psi`  | poloidal flux function      |
| 7           | `by`   | out-of-plane field          |

Reconstructed field: `B = ∇ψ × ŷ + B_y ŷ + B₀` (Cartesian), or `∇ψ/(R+ε) × ŷ + …` (cylindrical).

### Weak-form residuals assembled in `nlfun_apply`

**Continuity** (lines 861–866):
```
F_n += φ n + dt [φ (v·∇n + n∇·v) + ∇φ · D∇n]
```
Fully compressible. `∇·v ≠ 0` term explicitly present.

**Momentum** (lines 868–906):
```
F_v_k += φ v_k + dt [ (B·∇φ)B_k/μ₀mn - B²∇φ_k/2μ₀mn   (Lorentz)
                     + φ 2k_B∇(nT)_k/m                    (pressure gradient from EOS)
                     + φ (v·∇v_k)                          (advection)
                     + ν ∇φ·∇v_k ]                         (viscosity)
```
Pressure is the thermodynamic quantity `p = 2nkT/m` (no separate DOF).

**Temperature** (lines 916–921):
```
F_T += φ T/(γ-1) + dt [ φ v·∇T/(γ-1) + φ T∇·v + χ ∇φ·∇T ]
```
Compressible adiabatic energy equation.

**Induction for ψ** (lines 932–936):
```
F_ψ += φ ψ + dt [ φ (v·∇ψ + (v×B₀)_y) + η ∇φ·∇ψ ]
```
The `v·∇ψ` term is the self-field motional EMF — fully nonlinear, not lagged.

**Induction for By** (lines 950–955):
```
F_By += φ B_y + dt [ -φ (∇ψ×∇v_y)_y + φ v·∇B_y + φ B_y∇·v + η ∇φ·∇B_y ]
```

Sign convention: `F(u_new) = F_old = M·u_old`. Implicit terms appear on the LEFT
(positive sign in F). Dissipative terms (viscosity, drag) are POSITIVE in F.

---

## C. FE Representation and Field Layout

- FE space: `oft_scalar_bfem` (2D scalar Lagrange on `oft_bmesh` surface mesh)
- Composite: `oft_fem_comp_type` with `nfields=7`, all in the same Lagrange space
  (setup, lines 1568–1587)
- Quadrature: `oft_blagrange%quad` (shared quad object)
- Test functions: `oft_blag_eval` / `oft_blag_geval` (basis value and gradient)
- Geometry: `mesh%jacobian(i,pt,jac_mat,jac_det)`, `mesh%log2phys(i,pt)` for coords

Cylindrical mode: `cyl_flag=.TRUE.` adds `coords(1)` (= R) weighting to all
volume integrals and modifies `div_vel` to include the `vel(1)/R` geometric term.
Axis regularization: `gs_epsilon` from `USE oft_gs`.

---

## D. Residual / Jacobian Structure

- **Residual** `nlfun_apply(self, a, b)`: evaluates `F(u_new)` at current Newton
  iterate `a`. Two nested DO loops (cell, quadrature), one reconstruction DO jr
  loop, one residual accumulation DO jr loop.
- **Approximate Jacobian** `build_approx_jacobian(self, a)`: assembles
  `∂F/∂u` at `a`. Dense 7×7 block local matrix assembled per cell, scattered
  atomically to global. Uses `dt_fac = self%jac_dt` (= `dt` for BDF1, `dt/2` for CN).
- **Newton-Krylov driver**: `oft_nksolver` with GMRES inner solver. Jacobian
  rebuilt every `pre_freq` steps. Matrix-free option via `oft_mf_matrix`.

---

## E. BC Structure

Velocity, density, temperature, psi, and by BCs are boolean arrays:
```fortran
velx_bc(node) = .TRUE.  → homogeneous Dirichlet (zero)
velx_bc(node) = .FALSE. → natural BC (zero traction from weak form)
```
Default (`setup_bc`, lines 1602–1610): all fields get Dirichlet on all exterior
nodes (`oft_blagrange%global%gbe`). No Robin, partial-slip, or traction BCs exist.

---

## F. Input / Config Parsing

Parameters are set on `oft_xmhd_2d_sim` by the calling Fortran driver program.
Conventions (from `test_alfven_2d.F90`, `harris_sheet.F90`):
- Driver reads a Fortran NAMELIST (e.g., `xmhd_options`) from `oft.in`
- Driver assigns `mhd_sim%eta = eta`, `mhd_sim%nu = nu`, etc.
- Solver XML config for preconditioner via `oft_in.xml` element `<xmhd2d>`

New parameters should follow the same pattern: add to `oft_xmhd_2d_sim`, document
in this file, have example drivers add to their NAMELIST.

---

## G. Existing Test Infrastructure

Tests in `src/tests/physics/`:
- `test_alfven_2d.F90` / `test_alfven_2d.py` — Alfvén wave round-trip accuracy
- `test_sound_2d.F90` / `test_sound_2d.py` — compressible sound wave accuracy
- Registered in `CMakeLists.txt` via `oft_add_test(test_alfven_2d.F90 test_phys_helpers)`
- Python driver uses `run_OFT("./test_alfven_2d", nproc, timeout)` and reads a
  `.results` file written by the Fortran program

New tests must:
1. Add `oft_add_test(test_NEW.F90)` to `CMakeLists.txt`
2. Add `configure_file(test_NEW.py ...)` and any data files
3. Write a Fortran driver that writes a `.results` file
4. Write a pytest file that calls `run_OFT` and validates `.results`

---

## H. Recommended Implementation Path

### Phase 1 — Reduced Hartmann / wall drag (LOW RISK, HIGH VALUE)

Add optional volumetric momentum sink `S_u = -α_drag · u_⊥` to the existing
compressible solver. Controlled by `use_wall_drag = .FALSE.` (default-off).
This is the minimum viable LM-MHD physics extension for the public solver.

**Files:** `xmhd_2d.F90`, `test_drag_decay.F90`, `test_drag_decay.py`, `CMakeLists.txt`

### Phase 2 — Scalar-to-tensor conductivity infrastructure (MEDIUM RISK)

Add optional per-region scalar coefficients and a conductivity tensor hook.
Reuse the `xmhd.F90` pattern of `eta_reg(mesh%nreg)` for region-specific η.
Then generalize the ψ diffusion term to support a 2×2 symmetric η-tensor.

**Files:** `xmhd_2d.F90` (larger edit), possibly new `lm_material_utils.F90`

### Phase 3 — Incompressible LM-MHD formulation (HIGH RISK)

The TOFE abstract specifies **incompressible** MHD for the liquid-metal region.
This is architecturally incompatible with the current compressible solver.

Two sub-options:
- **3a (recommended):** New module `src/physics/xmhd_2d_lm.F90` — incompressible
  velocity-pressure formulation with Taylor-Hood P2/P1 or stabilized P1/P1.
- **3b (investigation needed):** Add incompressible mode to existing module by
  adding a pressure DOF (8th field) and a div-u=0 constraint residual. Risky
  because it changes the FE block structure.

**For inductionless (low-Rm) formulation:** evolve electric potential φ instead
of ψ/B directly. Much more appropriate for Li/PbLi blanket flows where Rm ≪ 1.

### Phase 4 — Wall / current-closure BCs

Add insulating and conducting-wall electromagnetic BCs for the inductionless solver.
Insulating wall: `J·n = 0` (naturally satisfied if ∂φ/∂n = -(u×B)·n at wall).
Conducting wall: thin-wall conductance ratio C = σ_w t_w / (σ_LM a).

### Phase 5 — TokaMaker bidirectional coupling

Python-level operator-split coupling loop. Requires:
- API to read TokaMaker equilibrium into MUG background field
- API to export LM current distribution from MUG
- Update to `gs_save_mug` or companion `gs_load_mug_currents`

---

## I. Risk Assessment

| Phase | Risk | Rationale |
|-------|------|-----------|
| 1 (wall drag) | Low | Pure residual/Jacobian addition; default-off flag; no structural change |
| 2 (tensor η) | Medium | Requires generalizing the ψ diffusion term; affects Jacobian block structure |
| 3a (new incompressible module) | Medium-High | New FE formulation; saddle-point solver; no mixed-space infrastructure in 2D yet |
| 3b (incompressible mode in existing) | High | Changes field count and block structure; risks breaking existing tests |
| 4 (wall BCs) | Medium | Depends on Phase 3 being complete first |
| 5 (TokaMaker coupling) | High | Requires Python coupling loop and new HDF5 I/O contracts |

---

## J. Minimum Viable Staged Upgrade Plan

**PR 1 (this branch, Phase 1):**
- Add `use_wall_drag`, `drag_coeff`, `drag_bhat(3)` to `oft_xmhd_2d_sim`
- Add drag residual and Jacobian to `nlfun_apply` / `build_approx_jacobian`
- Add `test_drag_decay` Fortran+Python test
- Add this plan document and a validation plan document

**PR 2 (Phase 2, future):**
- Add optional per-region η support reusing 3D MUG pattern
- Add 2×2 η-tensor to ψ diffusion: `∇φ^T · ETA · ∇ψ`
- Add manufactured anisotropic diffusion test

**PR 3 (Phase 3a, future):**
- New `xmhd_2d_lm.F90` with incompressible velocity-pressure fields
- Inductionless formulation with φ potential
- Hartmann channel flow benchmark example

**PR 4 (Phase 4+5, future):**
- Wall BC hooks for current closure
- TokaMaker coupling API and example script

---

## Appendix: Key Line Numbers in `xmhd_2d.F90`

| Item | Lines |
|------|-------|
| `oft_xmhd_2d_sim` type | 69–124 |
| Scalar parameters (η, ν, χ, D) | 82–92 |
| `nlfun_apply` start | 717 |
| Local variable declarations | 722–726 |
| Parameters copied from parent_sim | 748–755 |
| BLOCK + OMP private clause | 769–779 |
| Quadrature loop start | 798 |
| Reconstruction DO jr | 814–829 |
| B-field / div_vel reconstruction | 830–836 |
| Residual DO jr start | 852 |
| Continuity residual | 861–866 |
| Momentum residual (cyl) | 868–890 |
| Momentum residual (Cartesian) | 892–906 |
| Temperature residual | 907–921 |
| ψ residual | 923–936 |
| B_y residual | 938–955 |
| `build_approx_jacobian` start | 1010 |
| Local var decls (Jacobian) | 1013–1017 |
| Parameters from self | 1040–1048 |
| Jacobian BLOCK + OMP private | 1057–1070 |
| Jacobian vel,vel diagonal | 1232–1265 |
| `setup` (field registration) | 1537–1597 |
| `setup_bc` (default BCs) | 1602–1611 |
