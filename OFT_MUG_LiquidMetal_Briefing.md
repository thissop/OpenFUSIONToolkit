# OFT / MUG Liquid-Metal MHD: Full Technical Briefing

**Prepared for:** Claude on the web  
**Project:** OpenFUSIONToolkit (OFT), MUG physics solver  
**Repo:** https://github.com/OpenFUSIONToolkit/OpenFUSIONToolkit  
**Context date:** June 2026  
**Prepared by:** Repository inspection via Claude Code (5 parallel agents, ~12,000 lines of Fortran read)

---

## Who is involved

- **Sophia Guizzo** — PhD student, Columbia University. Author of the TOFE 2025 and APS DPP abstracts on axisymmetric MUG + liquid-metal blanket coupling. Email: s.guizzo@columbia.edu
- **Chris Hansen** — Columbia, OFT principal developer
- **S. Freiberger, C. Paz-Soldan** — Columbia collaborators on the blanket project
- **Yuchen / Sergey** — Other OFT maintainers (referenced in internal context)
- **Thaddaeus Kiker** — Columbia, writing this briefing (tjk2147@columbia.edu)

---

## What OFT/MUG is

The **OpenFUSIONToolkit** is a research-grade finite-element MHD code suite written mostly in Fortran 90 with Python wrappers. It contains:

- **MUG** — the MHD physics solver. Solves extended/visco-resistive MHD on unstructured meshes.
- **TokaMaker** — a 2D time-dependent plasma equilibrium solver (Grad-Shafranov based).
- **ThinCurr** — a thin-wall eddy-current solver (for passive conducting structures).

MUG has three Fortran physics modules:

| File | Lines | Basis | Description |
|---|---|---|---|
| `src/physics/xmhd.F90` | 5077 | H(Curl)+H1+Lagrange (mixed) | 3D extended-MHD, curl-conforming B |
| `src/physics/xmhd_lag.F90` | 4929 | All-Lagrange | 3D extended-MHD, Lagrange B with divergence cleaning |
| `src/physics/xmhd_2d.F90` | 1811 | All-Lagrange scalar | **New** 2D/2.5D MHD (commit 5faae9a) |

**Important:** "lag" in `xmhd_lag.F90` means **Lagrange basis**, NOT time-lagged. Both 3D modules are physically equivalent.

---

## The Sophia Guizzo TOFE 2025 abstract (most recent, most authoritative)

> **"Flexible, Axisymmetric MHD Solver for Coupled Plasma and Liquid Metal Blanket Dynamics"**  
> S. Guizzo, C. Hansen, S. Freiberger, C. Paz-Soldan — Columbia University

Key claims:
1. A **new capability within OFT** for simulating axisymmetric plasma evolution coupled to static and/or flowing conducting regions.
2. The **flowing conductor regions** (liquid lithium, PbLi) are treated with **time-dependent visco-resistive incompressible MHD equations**.
3. The **plasma region** neglects inertia and solves an **MHD equilibrium at each timestep** (TokaMaker-style), allowing simulation on the liquid-metal flow timescale.
4. Supports **anisotropic conductivity**.
5. Geometry features: **poloidal ducts**, **toroidally isolated blanket segments**.
6. Applied to **disruption scenarios**: predicts induced velocities, pressure variations, and blanket feedback on plasma evolution.

**Critical word: INCOMPRESSIBLE.** This is fundamental — see Section 4 below.

---

## What is actually in the public main branch (code-confirmed)

### 1. `xmhd_2d.F90` — what it really solves

This is a **compressible** 2D/2.5D MHD solver. It evolves **seven coupled scalar fields**, all in scalar Lagrange FE space:

| Field | Variable | Physical quantity |
|---|---|---|
| 1 | `n` | number density |
| 2 | `vel(1)` | velocity x (or R in cylindrical) |
| 3 | `vel(2)` | velocity y (toroidal/out-of-plane) |
| 4 | `vel(3)` | velocity z |
| 5 | `T` | temperature |
| 6 | `psi` | poloidal magnetic flux function |
| 7 | `by` | out-of-plane magnetic field |

The **magnetic field** is reconstructed as:
```
B = ∇ψ × ŷ + By ŷ + B₀           [Cartesian]
B = ∇ψ/(R+ε) × ŷ + By ŷ/(R+ε) + B₀  [Cylindrical, cyl_flag=T]
```

### The actual PDEs (from weak-form residual `nlfun_apply`, lines 717–1006):

**Continuity:**
```
∂n/∂t + v·∇n + n∇·v = D∇²n
```
Fully compressible — density evolves. `∇·v ≠ 0`.

**Momentum:**
```
ρ(∂v/∂t + v·∇v) = J×B − ∇p + ν∇²v
```
where `p = 2nkT/m` (factor 2 for two-species pressure p_i + p_e). **There is no separate pressure DOF.** No saddle-point system. No `∇·v = 0` constraint.

**Temperature:**
```
(∂T/∂t + v·∇T)/(γ−1) + T∇·v = χ∇·(∇T/n)·n
```
Compressible energy equation. No internal energy / enthalpy formulation.

**Induction for ψ (poloidal flux):**
```
∂ψ/∂t + v·∇ψ + (v×B₀)_y = η∇²ψ
```
The `v·∇ψ` term IS the self-field motional EMF (equivalent to `(v×B_in-plane)_y`). This term is **fully nonlinear** in the Newton residual — it is NOT lagged or linearized.

**Induction for By (toroidal field):**
```
∂By/∂t + v·∇By + By∇·v − (∇ψ×∇vy)_y = η∇²By
```

### Key physics flags (lines 69–91):

```fortran
LOGICAL :: cyl_flag = .FALSE.    ! T = R-Z cylindrical, F = Cartesian
LOGICAL :: linear   = .FALSE.    ! T = linearized perturbation mode
LOGICAL :: timestep_cn = .FALSE. ! T = Crank-Nicolson, F = backward Euler
LOGICAL :: mfnk     = .FALSE.    ! T = matrix-free Newton-Krylov
REAL(r8) :: nu      = ...        ! kinematic viscosity [m²/s], uniform scalar
REAL(r8) :: eta     = ...        ! resistivity / μ₀ [m²/s], uniform scalar
REAL(r8) :: chi     = ...        ! thermal diffusivity, uniform scalar
REAL(r8) :: D_diff  = ...        ! density diffusion coefficient
REAL(r8) :: B_0(3)  = [0,0,0]   ! static background magnetic field [T]
```

All material coefficients are **uniform scalars** — no spatial variability.

### Nonlinear solver:

Newton-Krylov with GMRES inner solver and approximate block-Jacobian preconditioner. Jacobian rebuilt every 4 Newton steps (line 242: `nksolver%up_freq=4`). Adaptive timestep control targeting 40 GMRES iterations/step.

### Boundary conditions for velocity:

**Only homogeneous Dirichlet (no-slip) or natural (zero-traction).** Applied via boolean arrays `velx_bc`, `vely_bc`, `velz_bc`. Default (`setup_bc`, lines 1602–1607): all boundary nodes get zero-velocity Dirichlet (no-slip). No Robin, no traction, no partial-slip BC exists.

### Critical limitation — cylindrical linear solve:

```fortran
! xmhd_2d.F90:378–379
IF(self%cyl_flag)CALL oft_abort( &
    "Cylindrical mode not currently supported for linear solve","run_lin_simulation",__FILE__)
```
Linear perturbation mode aborts in cylindrical coordinates. Only nonlinear advance supports R-Z.

---

### 2. `xmhd.F90` — the 3D module (summary)

The 3D code is much more capable:

- Mixed H(Curl)+H1+Lagrange FE basis for B (automatically divergence-free)
- Full extended-MHD: Hall term (off by default), two-temperature Ti/Te model, Braginskii transport
- Three viscosity types: kinematic (`kin`), isotropic Newtonian (`iso`), anisotropic (`ani`)
- Temperature-dependent Spitzer resistivity: `η ∝ T^{-3/2}`
- Region-specific resistivity via XML `<region>` blocks (`xmhd_setup_regions`, lines 2591–2675)
- `solid_cell` / `xmhd_rw` mechanism: marks type-2 mesh regions as zero-velocity conducting solids (resistive wall mode)
- `vbc='all'/'norm'/'tang'/'none'` velocity BC options (normal or tangential Dirichlet component selection)
- Gravity term (`g_accel`)

**The `solid_cell` wall model in 3D:** Type-2 regions have velocity forced to zero and skip the momentum solve entirely. Only resistive B-diffusion (`dB/dt = -∇×(ηJ)`) is solved in solid cells. This is a **passive conducting wall** — not a dynamic liquid-metal fluid region.

**No multi-region viscosity, no multi-region density, no interface conditions** between fluid regions.

---

### 3. TokaMaker coupling — what actually exists

One-way HDF5 export only:
- `gs_save_mug()` at `src/physics/grad_shaf.F90:4851` — writes BR, BT, BZ, pressure, PSI to HDF5
- Python wrapper: `TokaMakerCore.save_mug()` at `src/python/OpenFUSIONToolkit/TokaMaker/_core.py:3012`
- `xmhd_2d.F90` can ingest the background field via `B_0` parameter

**No reverse coupling** (MUG → TokaMaker). No coupled timestep loop. No blanket current feeding back into plasma equilibrium. This is the core missing piece of the TOFE capability.

---

### 4. The central discrepancy

| | TOFE abstract claims | What's in `xmhd_2d.F90` |
|---|---|---|
| **Fluid model** | Incompressible MHD | **Compressible** MHD |
| **Pressure** | Separate pressure DOF (required for incompress.) | Thermodynamic `p=2nkT/m`, no DOF |
| **∇·v=0 constraint** | Required | Absent |
| **Conductivity** | Anisotropic tensor | Scalar uniform η |
| **TokaMaker coupling** | Bidirectional, plasma evolves on LM timescale | One-way file export |
| **Multi-region** | LM region + plasma region | Uniform material only |
| **Code status** | "new capability" at time of abstract | Not in public repo |

**Conclusion:** The Guizzo TOFE blanket solver is a **separate development** from `xmhd_2d.F90`. It is not on any public branch. It presumably lives on a private branch or is under active development.

---

## What IS validated and committed

### Tests using `xmhd_2d.F90`:

1. **2D Alfvén wave test** (`src/tests/physics/test_alfven_2d.F90`)  
   — Propagates Alfvénic perturbation against guide field `B₀`  
   — n and T frozen via Dirichlet; small but nonzero η, ν, χ  
   — Checks waveform returns to initial values (round-trip accuracy)  
   — Tolerance ~3% error at polynomial orders 2–4

2. **2D sound wave test** (`src/tests/physics/test_sound_2d.F90`)  
   — Propagates compressible sound wave with near-ideal parameters  
   — No magnetic field  
   — Checks n, v, T waveform errors

3. **Harris sheet / tearing mode** (`src/examples/MUG/harris_sheet/harris_sheet.F90`)  
   — Classic GEM-challenge 2D tearing mode problem  
   — Uses `xmhd_2d` with full η, ν, χ, D_diff  
   — No guide field, periodic-x domain, no-slip top/bottom walls  
   — Nonlinear, fully compressible

### 3D MUG examples (use `xmhd.F90`, not `xmhd_2d`):

- `Spheromak_tilt` — tilt instability in spheromak geometry
- `Spheromak_heating` — Ohmic + viscous heating in spheromak
- `slab_recon` — slab reconnection

---

## Complete inventory of missing LM-MHD physics

### Definitively absent (zero hits in all `.F90` files):

- `liquid`, `lithium`, `blanket`, `PbLi`
- `Hartmann`, `Shercliff`, `Hunt flow`, `duct`
- `drag`, `friction`, `wall_drag`
- `incompress`, `div.*free`, `pressure.*solve`, `saddle.*point`
- `anisotropic.*conduct`, `sigma.*tensor`, `conductivity.*tensor`
- `conducting.*wall`, `thin.*wall`, `wall.*sigma`
- Any Hartmann flow, Shercliff flow, or duct-flow benchmark

### Present but simplified:

- `eta`, `nu`, `chi` — scalar uniform, no spatial variation
- Velocity BCs — Dirichlet only, no traction/Robin
- Cylindrical mode — nonlinear OK, linear aborts
- Cylindrical viscous stress — geometrically simplified (comment at `xmhd_2d.F90:628`)
- `solid_cell` wall — zero-velocity EM wall only, not a dynamic fluid

---

## Physics model hierarchy for liquid-metal MHD

For reference, here is the correct model hierarchy for LM-MHD in fusion blanket contexts:

### Tier 1: What `xmhd_2d.F90` can do now (with caveats)

The code solves visco-resistive compressible MHD. For liquid metals (Ma ≪ 1), running compressible is wasteful (must resolve acoustics) but may give qualitative answers if the Mach number is low enough. The `u×B` EMF is fully present and nonlinear. The J×B Lorentz force is correct.

**Use case:** Rough qualitative study of LM flow induced by transient magnetic fields, if Mach number constraints are satisfied.

### Tier 2: Add incompressible solver

For liquid metals, the proper formulation is:
```
∂u/∂t + u·∇u = (1/ρ)(J×B) − (1/ρ)∇p + ν∇²u
∇·u = 0
J = σ(E + u×B)
∇·J = 0  (charge conservation / current closure)
```
This requires a **pressure DOF** in the FE space (H1 elements for pressure or mixed Taylor-Hood `P2/P1`) and a saddle-point solver. The current `xmhd_2d.F90` lacks all of this.

For the inductionless (low-Rm) approximation, B is externally imposed and only the electric potential φ is evolved:
```
J = σ(-∇φ + u×B)
∇·J = ∇·(σ(-∇φ + u×B)) = 0  → elliptic equation for φ
```

### Tier 3: Add Hartmann friction (reduced model)

For unresolved Hartmann layers (Ha ≫ 1), add a volumetric momentum sink:
```
S_u = -α_H · σ · B₀² · u_⊥ / (ρ · L_z)
```
where `L_z` is the duct half-width in the B direction and `α_H ~ 1` is an O(1) coefficient. This is equivalent to adding `-α_drag · u` to the momentum residual. In normalized units, `α_H = Ha / (√3 L_z)` for a Hartmann duct.

Implementation in `xmhd_2d.F90:nlfun_apply` (lines 892–906): add ~10 lines to the momentum residual and ~10 lines to the Jacobian block in `build_approx_jacobian`.

### Tier 4: Anisotropic conductivity / wall conductance

The generalized Ohm's law in the presence of a strong B-field takes the form:
```
J = σ_⊥(E + u×B) + σ_H(E + u×B)×b̂/B + (σ_‖ - σ_⊥)(E·b̂)b̂
```
where `b̂ = B/|B|` is the field direction. In the inductionless limit, this becomes an anisotropic conductivity tensor in the potential equation.

### Tier 5: Conducting-wall electromagnetic coupling

At the LM-wall interface, the boundary conditions depend on whether the wall is:
- **Insulating:** `J·n̂ = 0` at wall (no current passes through)
- **Conducting (thin wall):** Current passes through via wall conductance ratio `C = σ_w t_w / (σ_LM a)` where `t_w` is wall thickness, `a` is duct half-width

Conducting-wall BCs create electromagnetic coupling between the LM flow and the wall that qualitatively changes the flow profile (Hunt flow vs. Shercliff flow vs. Hartmann flow).

### Tier 6: Full TokaMaker coupling (the Guizzo TOFE capability)

The plasma region treats MHD as quasi-static equilibrium (∂/∂t terms dropped, only Grad-Shafranov balance). The LM region evolves on its flow timescale. At each LM timestep:
1. TokaMaker updates plasma equilibrium given current plasma shape and external coils
2. Updated B-field from TokaMaker is passed as boundary condition or source to LM solver
3. LM solver advances one timestep (or several) 
4. LM currents are fed back to TokaMaker (if feedback is included)
5. Repeat

This requires a coupled-solver API that does not currently exist in the public repo.

---

## Specific code locations for potential changes

### Adding Hartmann drag source term to `xmhd_2d.F90`

**Where to add it:** `nlfun_apply` subroutine, lines 892–906 (momentum residual loop)

**What to add (pseudocode):**
```fortran
! Existing momentum residual assembly ends around line 906
! Add after the viscous term:
IF(self%use_wall_drag) THEN
  drag_fac = self%drag_coeff * int_factor * self%dt
  DO k = 1, 3
    IF(k /= self%drag_dir) THEN   ! apply only to components perpendicular to B
      res_loc(jr, k+1) = res_loc(jr, k+1) + drag_fac * vel(k)
    END IF
  END DO
END IF
```

**Corresponding Jacobian entry** in `build_approx_jacobian` around lines 1192–1260:
```fortran
! In the vel_k / vel_k diagonal block:
IF(self%use_wall_drag .AND. k /= self%drag_dir) THEN
  jac_loc(k+1,k+1)%m(jr,jc) = jac_loc(k+1,k+1)%m(jr,jc) &
    + dt_fac * drag_fac * basis_vals(jc) * int_factor
END IF
```

**New parameters needed on the `oft_xmhd_2d_sim` type:**
```fortran
REAL(r8) :: drag_coeff = 0.d0      ! α_H σB²/(ρ L) in code units
INTEGER :: drag_dir = 2             ! direction of B (1=x, 2=y/toroidal, 3=z)
LOGICAL :: use_wall_drag = .FALSE.  ! enable/disable
```

### Adding incompressible pressure DOF (larger change)

This requires adding an 8th field to the composite FE representation:
- `setup` subroutine (lines 1537–1597): add `p` Lagrange DOF
- `nlfun_apply`: replace `p = 2nkT/m` with actual pressure variable; add `∇·v = 0` constraint equation (mass residual)
- `build_approx_jacobian`: add pressure-velocity coupling blocks (saddle point)
- FE spaces: velocity in `P2` (quadratic), pressure in `P1` (linear) for Taylor-Hood stability, OR use stabilized equal-order formulation with pressure stabilization

This is a significant structural change that would likely be a new module.

---

## File locations reference

```
src/physics/xmhd_2d.F90         — 2D/2.5D MHD module (1811 lines)
src/physics/xmhd.F90             — 3D extended-MHD, H(Curl) basis (5077 lines)
src/physics/xmhd_lag.F90         — 3D extended-MHD, Lagrange basis (4929 lines)
src/physics/mhd_utils.F90        — Physical constants, Braginskii, Spitzer (198 lines)
src/physics/grad_shaf.F90        — Grad-Shafranov (TokaMaker), includes gs_save_mug at line 4851
src/physics/thin_wall.F90        — ThinCurr passive conductor
src/examples/MUG/harris_sheet/
  harris_sheet.F90               — 2D tearing mode example using xmhd_2d
  oft.in                         — Input parameters (eta=742.6, nu=7424.22, chi=2943.4)
  oft_in.xml                     — Solver config (GMRES + block-Jacobi + LU)
src/tests/physics/
  test_alfven_2d.F90             — Alfvén wave test for xmhd_2d
  test_sound_2d.F90              — Sound wave test for xmhd_2d
src/python/OpenFUSIONToolkit/TokaMaker/
  _core.py                       — Python TokaMaker driver, save_mug() at line 3012
src/docs/
  doc_mug_main.md                — MUG documentation
  doc_tokamaker_main.md          — TokaMaker documentation
```

---

## Key line numbers in `xmhd_2d.F90` for quick navigation

| Location | Line(s) | Content |
|---|---|---|
| Module parameters | 69–95 | All physics flags and material coefficients |
| `cyl_flag` | 70 | Cartesian/cylindrical switch |
| `run_simulation` | 136–346 | Nonlinear time advance driver |
| `run_lin_simulation` | 350–526 | Linear perturbation advance |
| Cylindrical abort | 378–379 | Linear solve blocked for cyl_flag=T |
| `nlfun_apply` (core residual) | 717–1006 | All PDE weak forms |
| Continuity residual | 861–866 | `∂n/∂t + v·∇n + n∇·v = D∇²n` |
| Momentum residual | 868–906 | Lorentz + pressure + advection + viscosity |
| Temperature residual | 916–921 | Compressible energy equation |
| ψ induction residual | 932–936 | `∂ψ/∂t + v·∇ψ + (v×B₀)_y = η∇²ψ` |
| By induction residual | 950–955 | Out-of-plane field induction |
| Dirichlet BC application | 987–993 | Velocity BCs applied (Dirichlet only) |
| `build_approx_jacobian` | 1010–1516 | Jacobian assembly |
| `setup` | 1537–1597 | FE representation setup |
| `setup_bc` | 1602–1611 | Default BC assignment (`gbe` = all boundary nodes) |

---

## Summary: what to ask Claude on the web

Given this briefing, productive questions for Claude on the web include:

### Physics/formulation questions:
- "Given the compressible xmhd_2d.F90 formulation, under what conditions does using it for liquid-metal MHD give valid results? What Mach number criterion must be satisfied?"
- "Design an incompressible pressure-velocity FE formulation for liquid-metal MHD that could be added as an 8th field (pressure) to xmhd_2d.F90. What FE spaces ensure inf-sup stability?"
- "What is the correct Hartmann friction coefficient for a rectangular duct with conducting walls (wall conductance ratio C) in terms of Ha, C, and duct geometry?"
- "Write the weak form for adding a reduced Hartmann drag source term `S_u = -α n u_⊥` to the xmhd_2d momentum residual."

### Implementation questions:
- "Given the Fortran structure of `nlfun_apply` in xmhd_2d.F90 [paste relevant lines], show exactly how to add a volumetric drag term `S = -alpha * vel(k) * n` with its Jacobian contribution."
- "What is the correct anisotropic Ohm's law tensor for a liquid-metal in a strong magnetic field, and how should it be discretized in a finite-element context?"

### Coupling / architecture questions:
- "Design a Python-level coupling loop between TokaMaker (plasma equilibrium) and xmhd_2d (liquid-metal MHD). What is the minimum API needed? What information crosses the interface at each step?"
- "How should the boundary condition between the plasma region (TokaMaker, no inertia) and the liquid-metal region (xmhd_2d, full MHD) be formulated? What is the electromagnetic interface condition?"

### Validation / benchmarking questions:
- "What is the analytical solution for Hartmann channel flow (velocity profile, pressure drop, current distribution) that could be used as a benchmark for xmhd_2d with wall drag added?"
- "How can I verify an incompressible MHD solver using the Hartmann flow benchmark? What are the key non-dimensional parameters and what should the velocity profile look like?"

---

## What NOT to assume

1. Do NOT assume the Guizzo TOFE/APS blanket solver is in the public repo — it is not.
2. Do NOT assume "xmhd_lag" means time-lagged physics — it means Lagrange basis.
3. Do NOT assume xmhd_2d is incompressible — it is explicitly compressible.
4. Do NOT assume wall drag exists anywhere in the codebase — it does not.
5. Do NOT assume TokaMaker coupling is bidirectional — it is one-way file export only.
6. Do NOT assume there are LM-MHD benchmarks in the repo — there are none.

---

## Unanswered questions (for Sophia / maintainers)

1. Is the Guizzo TOFE capability on a private branch? When will it be upstreamed?
2. Is the liquid-metal solver genuinely incompressible, or is it a low-Mach compressible run?
3. What FE spaces are used for pressure in the incompressible formulation?
4. How is the anisotropic conductivity tensor implemented? What variables does it depend on?
5. How does the TokaMaker feedback coupling work at the code level (inner loop, outer loop, operator-split)?
6. What wall BC is used — insulating, conducting, or mixed?
7. Has the cylindrical linear solve limitation (line 378) been fixed in the private branch?
8. Are there any Hartmann flow or duct-flow benchmarks in the private branch?
