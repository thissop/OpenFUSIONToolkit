# LM-MHD Upgrade Lima Handoff

**Created:** 2026-06-02  
**Branch:** `lm-mhd-upgrades`  
**Repo:** OpenFUSIONToolkit (local checkout at `/Users/tkiker/Documents/GitHub/OpenFUSIONToolkit`)

---

## 1. Current Branch and Status

```
Branch:   lm-mhd-upgrades
Base:     main (commit 1a4bc3c "tokamaker: add handling for out of bounds index")
HEAD:     6c80aaf "tests: add drag decay regression test for xmhd_2d wall drag"
```

**What is on this branch (3 commits ahead of main):**
- `9453740` — docs: add LM-MHD upgrade plan and validation plan
- `0dbdf81` — physics: add reduced Hartmann/wall-drag source to xmhd_2d
- `6c80aaf` — tests: add drag decay regression test for xmhd_2d wall drag

**Phase 1 is complete at the code level.** It adds an optional reduced Hartmann/wall-drag
source term to `src/physics/xmhd_2d.F90`. The feature is default-off and should not
affect any existing MUG 2D tests when disabled (`use_wall_drag = .FALSE.` is the default).

**The full runtime build and test has NOT yet been completed.** This work was done in a
macOS environment without a full OFT build directory. Only a partial `gfortran -fsyntax-only`
check was possible (it confirmed macro and module-path issues that are expected without
a full build — these are not code errors). Everything must be compiled and tested in
the Lima Ubuntu VM.

---

## 2. Why This Work Exists

The public `xmhd_2d.F90` is a compressible 2D/2.5D visco-resistive MHD solver. It has
plasma-MHD validation tests (Alfvén wave, sound wave, Harris sheet tearing mode) but no
committed liquid-metal MHD validation cases.

Sophia Guizzo's TOFE 2025 abstract describes a capability within OFT for
"axisymmetric plasma evolution coupled to static and/or flowing conducting regions"
with "time-dependent visco-resistive **incompressible** MHD equations" and
"anisotropic conductivity." This capability is **not present in the public main branch**.
The TOFE capability is on a private/development branch.

This Phase 1 patch adds a reduced Hartmann/wall-drag term as the **first low-risk,
LM-MHD-relevant upgrade** to the public standalone solver. It is:
- NOT full wall-resolved LM blanket MHD
- NOT incompressible (the solver remains compressible)
- NOT anisotropic conductivity
- NOT bidirectional TokaMaker coupling
- A physically motivated, architecturally clean first step

What remains as future work:
- Incompressibility (`∇·u = 0` constraint, pressure DOF) — high-risk new module
- Anisotropic/tensor conductivity — medium risk
- Wall current-closure and conducting-wall BCs — depends on incompressible solver
- Bidirectional TokaMaker coupling — API design + Python loop
- LM duct benchmarks (Hartmann channel flow, Shercliff, Hunt flow)

---

## 3. Implemented Feature: Reduced Hartmann/Wall Drag

### New fields on `oft_xmhd_2d_sim` (xmhd_2d.F90, lines 93–96)

```fortran
!--- Reduced Hartmann / wall-drag parameters (LM-MHD extension, default off)
LOGICAL :: use_wall_drag = .FALSE. ! Enable reduced wall/Hartmann drag source term
REAL(r8) :: drag_coeff = 0.d0     ! Drag coefficient alpha [1/time in problem units]
REAL(r8) :: drag_bhat(3) = [0.d0,1.d0,0.d0] ! Applied-field direction; drag NOT along this
```

**Set these from your driver program:**
```fortran
mhd_sim%use_wall_drag = .TRUE.
mhd_sim%drag_coeff    = alpha       ! [s^{-1} in problem units]
mhd_sim%drag_bhat     = [0.d0, 1.d0, 0.d0]  ! must be a unit vector (not normalized by code)
```

### Physics

The term models Hartmann-layer friction: velocity perpendicular to the applied field
is damped; velocity parallel to the field is not.

```
u_perp = u - (u · b̂) b̂
dv/dt += -alpha * u_perp
```

For uniform flow and pure drag (all other terms near-zero), this gives:
```
v_perp(t) = v_perp(0) * exp(-alpha * t)
```
Backward-Euler discrete:
```
v_perp_{n+1} = v_perp_n / (1 + alpha * dt)
```

Physical mapping (approximate, for reference):
```
alpha_Hartmann ~ sigma * B0^2 / (rho * L_Ha)
             ~ Ha * nu / L^2
```
where Ha = B₀ L √(σ/ρν) is the Hartmann number and L is the duct half-width
perpendicular to B. **This mapping is not in the code — document it externally.**

### Exact residual expression (xmhd_2d.F90, lines 919–934)

Residual sign convention: `F(u_new) = F_old = M·u_old`. Dissipative terms appear
with a **positive** sign in F (on the left side of `F = F_old`). The drag is a
dissipative term.

**Residual added per test function `jr`, per velocity component `k`:**

Cartesian (line 930–931):
```fortran
res_loc(jr,k+1) = res_loc(jr,k+1) &
  + drag_coeff*self%dt*basis_vals(jr)*(vel(k)-vdotb_drag*drag_bhat(k))*int_factor
```

Cylindrical (line 927–928):
```fortran
res_loc(jr,k+1) = res_loc(jr,k+1) &
  + drag_coeff*self%dt*basis_vals(jr)*(vel(k)-vdotb_drag*drag_bhat(k))*int_factor*coords(1)
```

where `vdotb_drag = DOT_PRODUCT(vel, drag_bhat)` is precomputed once per quadrature
point outside the test-function loop (line 862).

### Exact Jacobian expression (xmhd_2d.F90, lines 1288–1319)

```
dF_vel_k / dv_l = alpha * dt_fac * phi_jr * phi_jc * (delta_kl - bhat_k * bhat_l)
```

Cartesian (lines 1311–1318):
```fortran
IF(use_wall_drag) THEN
  DO k=1,3
    DO l=1,3
      jac_loc(k+1,l+1)%m(jr,jc) = jac_loc(k+1,l+1)%m(jr,jc) &
        + dt_fac*drag_coeff*basis_vals(jr)*basis_vals(jc) &
          *(merge(1.d0,0.d0,k==l) - drag_bhat(k)*drag_bhat(l))*int_factor
    END DO
  END DO
END IF
```

Cylindrical (lines 1288–1295): same but multiplied by `*int_factor*coords(1)`.

`dt_fac = self%jac_dt` which equals `dt` for backward Euler and `dt/2` for
Crank-Nicolson.

### Important caveats / TODOs

- ~~**`drag_bhat` is NOT normalized by the code.**~~ **RESOLVED 2026-06-02.** `drag_bhat`
  is now normalized once in `run_simulation`/`run_lin_simulation` (with a zero-vector abort).
  Users may supply any non-zero vector; see §6 fix 1.
- Drag is currently applied uniformly over the entire domain. A per-region drag
  (applying only in the liquid-metal region, not the plasma region) is future work.
- No input namelist parsing is added for the drag parameters. The user must set them
  programmatically in the Fortran driver. This is consistent with how `eta`, `nu`,
  etc. are handled, but should be documented clearly.

---

## 4. Changed Files

All 6 files exist and are committed on `lm-mhd-upgrades`:

| File | What changed |
|------|-------------|
| `src/physics/xmhd_2d.F90` | Added `use_wall_drag`, `drag_coeff`, `drag_bhat` to type; residual and Jacobian drag terms in both Cartesian and cylindrical paths |
| `src/tests/physics/test_drag_decay.F90` | New 223-line Fortran test program |
| `src/tests/physics/test_drag_decay.py` | New 122-line pytest driver |
| `src/tests/physics/CMakeLists.txt` | 3 lines added to register the new test |
| `docs/dev/lm_mhd_upgrade_plan.md` | 264-line architecture and phased implementation plan |
| `docs/dev/lm_mhd_validation_plan.md` | 121-line validation roadmap with analytic benchmarks |
| `docs/dev/lm_mhd_lima_handoff.md` | This handoff file |

---

## 5. Exact Code Locations

### New type fields (`xmhd_2d.F90`)

```
Line 94:  LOGICAL :: use_wall_drag = .FALSE.
Line 95:  REAL(r8) :: drag_coeff = 0.d0
Line 96:  REAL(r8) :: drag_bhat(3) = [0.d0,1.d0,0.d0]
```

### No input/namelist parsing for drag parameters

The drag parameters are set programmatically (like `eta`, `nu`, `chi`).
Verify by running:
```
grep -n "drag" src/physics/xmhd_2d.F90 | grep -i "namelist\|read\|input\|nml"
```
Expected: zero hits. **Action needed:** add parsing to example drivers and document.

### Local variable copies in `nlfun_apply`

```
Lines 731–732:  LOGICAL :: use_wall_drag / REAL(r8) :: drag_coeff, drag_bhat(3)
Lines 762–764:  Assigned from self%parent_sim%*
Line  784:      REAL(r8) :: vdotb_drag  (in BLOCK, listed in OMP private clause)
Line  789:      vdotb_drag in !$omp parallel private(...) clause
Line  862:      IF(use_wall_drag) vdotb_drag = DOT_PRODUCT(vel, drag_bhat)
Lines 919–934:  IF(use_wall_drag) drag residual loop
```

### Local variable copies in `build_approx_jacobian`

```
Lines 1042, 1047:  LOGICAL :: use_wall_drag / REAL(r8) :: drag_coeff, drag_bhat(3)
Lines 1080–1082:   Assigned from self%*
Lines 1288–1295:   Cylindrical drag Jacobian
Lines 1311–1318:   Cartesian drag Jacobian
```

### CMakeLists registration

```
src/tests/physics/CMakeLists.txt:30  oft_add_test( test_drag_decay.F90 )
src/tests/physics/CMakeLists.txt:31  configure_file( test_drag_decay.py ... COPYONLY )
```

### Test program key lines (`test_drag_decay.F90`)

```
Line  64: NAMELIST/drag_test_options/  (includes alpha_drag)
Line  79: READ(io_unit, drag_test_options, ...)
Lines 104-108: mhd_sim%use_wall_drag / drag_coeff / drag_bhat set
Lines 113-131: Velocity BC override (natural BCs) before run_simulation
Lines 144-170: Initial conditions set uniformly
Line  197: run_simulation() called
Lines 187-195: Final velocities read from mhd_sim%u after run
Line  201: velx_expected = velx0 / (1+alpha_drag*dt)^nsteps
Lines 204-205: velx_relerr, vely_relerr computed
Lines 208-211: Results written to drag_decay.results
Lines 213-219: PASS/FAIL + oft_abort on failure
```

### Pytest assertions (`test_drag_decay.py`)

```
Line  77: def validate_drag_result(velx_err_tol=0.02, vely_err_tol=1.e-6):
Line  80: reads drag_decay.results (2 lines: velx_err, vely_err)
Line  87: if velx_err > velx_err_tol → FAIL
Line  91: if vely_err > vely_err_tol → FAIL
Lines 101-104: test_drag_p2_r1() — alpha=2, order=2, dt=0.05, nsteps=10
Lines 110-113: test_drag_p2_r1_high_alpha() — alpha=5
Lines 117-119: test_drag_p3_r1() — alpha=2, order=3
```

---

## 6. Build/Test Status — ✅ COMPLETED (2026-06-02, Lima Ubuntu)

**Build: COMPLETED successfully.** Full CMake build in a Linux (Lima) environment.

**Environment:**
- OS: Ubuntu 25.10, aarch64
- Compilers: gcc/g++/gfortran 15.2.0
- CMake 3.31.6, Python 3.13 (venv at `oft_venv/`)
- External libs built from source via `build_libs.py`: MPICH 4.2.3, HDF5 1.14.6,
  OpenBLAS 0.3.30, METIS 5.1.0, SuperLU 7.0.0, UMFPACK 6.3.5, ARPACK-ng 3.9.1, libxml2 2.15.2
- LA backend: native; OpenMP on (6 threads)

**Exact build commands used:**
```bash
# Prerequisites (Ubuntu apt)
sudo apt-get install -y python3-venv python3-dev build-essential m4 autoconf automake libtool

# Python venv + deps
python3 -m venv oft_venv
echo "source $(pwd)/oft_venv/bin/activate" > setup_env.sh
source setup_env.sh
python -m pip install pytest numpy scipy h5py matplotlib xarray

# Stage 1: external libraries (~20 min, mirrors CI copilot-setup-steps.yml)
mkdir -p builds && cd builds
export CC=gcc CXX=g++ FC=gfortran
python ../src/utilities/build_libs.py --build_umfpack=1 --build_superlu=1 \
  --no_dl_progress --nthread=6 --build_arpack=1 --oft_build_tests=1 --build_mpich=1

# Stage 2: configure + build OFT
bash config_cmake.sh                 # generates + runs cmake into builds/build_release
cd build_release && make -j6         # NOTE: source the venv with an ABSOLUTE path first;
                                     # a relative `source setup_env.sh` fails in a fresh shell

# Run tests
cd tests/physics
python -m pytest test_alfven_2d.py test_sound_2d.py test_drag_decay.py -v
```

**Test results — ALL PASS:**

| Test | Purpose | Result |
|------|---------|--------|
| `test_alfven_2d` | Alfvén wave round-trip accuracy (drag default-off) | 9 passed |
| `test_sound_2d` | Compressible sound wave accuracy (drag default-off) | 9 passed |
| `test_drag_decay` | Backward-Euler drag decay rate (3 cases) | 3 passed |
| **Full suite** | `test_alfven_2d + test_sound_2d + test_drag_decay` | **21 passed, 36 skipped** |

(36 skips are MPI/coverage parametrizations, not failures.)

**`test_drag_decay` measured values:**

| Case | order | α | velx rel err (⊥, damped) | vely rel err (∥, undamped) |
|------|-------|---|--------------------------|-----------------------------|
| p2, α=2 | 2 | 2.0 | 3.52e-5 | 4.44e-16 |
| p2, α=5 | 2 | 5.0 | 3.42e-5 | 4.44e-16 |
| p3, α=2 | 3 | 2.0 | 5.17e-5 | 5.33e-13 |

All within tolerance (velx < 2%, vely < 1e-6). Perpendicular decay matches the
backward-Euler analytic rate; parallel component undamped to machine precision.

### Fixes made during build/test (2026-06-02)

1. **`drag_bhat` normalization (handoff §8 item 5, §11 known gap).** Added a one-time
   normalization with a zero-vector abort guard in both `run_simulation` and
   `run_lin_simulation` (`xmhd_2d.F90`, just after `CALL self%setup_bc()`):
   ```fortran
   IF(self%use_wall_drag)THEN
     IF(NORM2(self%drag_bhat)<=1.d-12)THEN
       CALL oft_abort('use_wall_drag=.TRUE. but drag_bhat is a zero vector', &
         'run_simulation',__FILE__)
     END IF
     self%drag_bhat=self%drag_bhat/NORM2(self%drag_bhat)
   END IF
   ```
   Users may now supply a non-unit `drag_bhat`; the perpendicular projection stays correct.

2. **Test segfault fix (`test_drag_decay.F90`).** `vec_vals` was declared
   `REAL(r8), POINTER :: vec_vals(:)` without `=> NULL()`. `vec_get_local` allocates only
   when `.NOT.ASSOCIATED(array)`, so the undefined pointer association caused a dangling-read
   SIGSEGV inside `vec_restore_local` on the first IC assignment. Fixed by initializing the
   pointer: `REAL(r8), POINTER :: vec_vals(:) => NULL()`. (Pure test-harness bug; no physics change.)

3. **Test timestep pinning (`test_drag_decay.F90`).** The order-3 case failed at 10.25%
   velx error because the adaptive timestep controller (`ittarget=40` default,
   `dt = ittarget·Σdt/Σlits` capped at `dt_initial`) shrank dt below 0.05 when the linear
   iteration count rose at higher FE order — so the constant-dt analytic prediction no longer
   held. Fixed by setting `mhd_sim%ittarget = 1000000`, which makes the controller cap dt at
   its initial value every step. p3 error dropped 10.25% → 5.17e-5. (Test-only change; the
   drag physics was always correct.)

**Files modified this session (uncommitted):**
- `src/physics/xmhd_2d.F90` — drag_bhat normalization (fix 1)
- `src/tests/physics/test_drag_decay.F90` — fixes 2 and 3

### Harris sheet example

`src/examples/MUG/harris_sheet/` (added by commit 5faae9a) exists but was **not built**:
the config used `-DOFT_BUILD_EXAMPLES:BOOL=FALSE`. It is a longer nonlinear tearing-mode run
and is not required for Phase 1 drag validation. To run it, reconfigure with
`OFT_BUILD_EXAMPLES=ON` (note `config_cmake.sh` wipes `build_release` on rerun).

---

## 7. Lima Setup Instructions

**Do NOT start or stop Lima VMs automatically.** You may have multiple Lima
instances running for different tasks. Before doing anything:

```bash
limactl list        # see what VMs are running
limactl shell <name>  # enter the correct one
```

Ask the user which Lima VM to use if there are multiple.

**Once inside the correct Lima VM:**

```bash
# Check environment
cmake --version           # need ≥ 3.18
gfortran --version        # need ≥ 9 for OpenMP and Fortran 2008+ features
python3 --version         # need 3.8+
mpirun --version || true  # OFT can build without MPI but it's preferred
h5cc -showconfig || true  # check HDF5 with parallel support

# Check repo state
cd /path/to/OpenFUSIONToolkit   # adjust to actual mount path
git status
git branch --show-current       # must be lm-mhd-upgrades
git log --oneline -5
```

**Read the build documentation FIRST:**
```bash
cat CONTRIBUTING.md
cat README.md
ls docs/
```

**Potential Ubuntu apt dependencies (confirm from OFT docs before installing):**
```bash
sudo apt-get update
sudo apt-get install -y \
  build-essential \
  gfortran \
  cmake \
  ninja-build \
  python3 \
  python3-pip \
  python3-pytest \
  libopenmpi-dev \
  openmpi-bin \
  libhdf5-openmpi-dev \
  libblas-dev \
  liblapack-dev \
  pkg-config
```

**CMake configure (example — adjust paths to match your HDF5 install):**
```bash
mkdir -p build && cd build
cmake .. \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DOFT_BUILD_TESTS=ON \
  -DOFT_BUILD_EXAMPLES=OFF \
  -DOFT_USE_MPI=ON \
  -G Ninja
ninja -j$(nproc)
```

**Do not overwrite user changes.** Run `git status` before any edit. The branch
`lm-mhd-upgrades` has 3 commits ahead of main — do not merge or reset.

---

## 8. Exact Next Actions for Fresh Claude Code in Lima

1. **Read this entire handoff file** (`docs/dev/lm_mhd_lima_handoff.md`).

2. **Confirm branch and state:**
   ```bash
   git status
   git branch --show-current    # expect: lm-mhd-upgrades
   git log --oneline -5
   ```

3. **Inspect the 6 changed files** listed in Section 4. Confirm they exist and look
   correct before building anything.

4. **Check drag parameter parsing:** The drag parameters are NOT in any input namelist
   in the physics module. They must be set programmatically in driver programs.
   Decide whether to add namelist parsing to `test_drag_decay.F90` (already reads
   `alpha_drag` from `drag_test_options` namelist — that's fine for the test) or
   whether to add general parsing guidance to example drivers. No action required
   for build, but document this limitation.

5. **Confirm `drag_bhat` normalization:** The code does NOT normalize `drag_bhat`.
   Add a normalization step or `CALL oft_abort` if `|drag_bhat| ≠ 1.0` to within
   tolerance. This is a small fix needed before claiming the feature is production-ready.

6. **Configure and build OFT:**
   ```bash
   mkdir -p build && cd build
   cmake .. -DOFT_BUILD_TESTS=ON -G Ninja
   ninja -j$(nproc)
   ```
   Fix any compilation errors before proceeding. If errors are in new code,
   inspect the sign convention and variable scope carefully.

7. **Run existing regression tests first:**
   ```bash
   cd build/tests/physics
   pytest test_alfven_2d.py -v
   pytest test_sound_2d.py -v
   ```
   Both must pass before running the new test. If they fail, the new code broke
   something despite the `use_wall_drag=.FALSE.` guard — investigate immediately.

8. **Run the new drag decay test:**
   ```bash
   pytest test_drag_decay.py -v
   ```
   Expected output: all three cases (`test_drag_p2_r1`, `test_drag_p2_r1_high_alpha`,
   `test_drag_p3_r1`) should PASS.

9. **If the drag test fails**, debug in this order:
   - Check `drag_decay.results` for the actual errors
   - Verify the sign convention: `F += dt*alpha*phi*(vel_k - vdotb*bhat_k)` means
     `(1 + alpha*dt)*v_perp_new = v_perp_old` (backward Euler)
   - Verify the OMP private clause includes `vdotb_drag`
   - Check whether `vdotb_drag` is initialized before use (the `IF(use_wall_drag)`
     guard means it's only assigned when drag is on, which is correct for the test)
   - Check the BC setup: velocity BCs are `.FALSE.` (natural), n and T are `.TRUE.`
     (Dirichlet fixed) — this is essential for the uniform-field analytic solution

10. **Add a Python postprocessing script if helpful:**
    ```python
    # Quick check of drag_decay.results
    with open('drag_decay.results') as f:
        velx_err = float(f.readline())
        vely_err = float(f.readline())
    print(f"velx relative error: {velx_err:.4e} (expect < 0.02)")
    print(f"vely relative error: {vely_err:.4e} (expect < 1e-6)")
    ```

11. **Update `docs/dev/lm_mhd_lima_handoff.md`** with test results (actual numbers,
    any fixes applied, compiler version used).

12. **Update `docs/dev/lm_mhd_validation_plan.md`** — move `test_drag_decay` from
    "Planned" to "Completed" with the actual measured error values.

13. **Only then begin OpenFOAM/FreeMHD validation** — see Section 9 below.

---

## 9. OpenFOAM / FreeMHD Validation Plan

> **A detailed, runnable cross-validation plan now lives in
> [`lm_mhd_openfoam_validation.md`](lm_mhd_openfoam_validation.md)** — including the
> analytic Hartmann reference, exact apt/`mhdFoam` install + run commands, the `Ha` sweep,
> and (importantly) an honest framing of what the reduced-drag model can and cannot be
> compared against (bulk damping yes; resolved `cosh` boundary-layer profile no).
> OpenFOAM was confirmed installable via the Ubuntu `openfoam` apt package but was **not
> installed/run this session** — the primary OFT Phase 1 goal (build + tests) is complete.
> The summary below is retained for context.

**Goal:** Use external MHD codes as cross-checks for LM-MHD-type behavior,
starting simple and escalating only if time allows.

**Order of priority:**

1. **Internal OFT drag-only decay** — validates the reduced drag model (Phase 1 test above)
2. **Internal OFT reduced Hartmann friction scaling** — vary α and dt, check that
   the decay rate scales correctly with both
3. **Analytic Hartmann channel-flow comparison** — classic 1D exact solution
4. **OpenFOAM `mhdFoam` Hartmann channel** — external cross-check
5. **FreeMHD** — only if OpenFOAM install succeeds and time permits

**Start with OpenFOAM `mhdFoam`, not FreeMHD.** OpenFOAM is likely easier to install
on Ubuntu and has a well-documented `mhdFoam` solver for incompressible MHD.

**Do NOT start this section until OFT builds and `test_drag_decay` passes.**

### OpenFOAM Hartmann flow setup (for reference)

Analytic solution — 1D Hartmann flow between insulating plates separated by 2a,
transverse field B₀, pressure gradient G = -dp/dx:

```
u_x(y) = (G a²)/(ν Ha) * [Ha/tanh(Ha) - cosh(Ha·y/a)/sinh(Ha)]
Ha = B₀ a sqrt(σ/(ρ ν))
Q / Q_Poiseuille = tanh(Ha) / Ha * 3/Ha²  (flow rate ratio)
```

For `mhdFoam` in OpenFOAM, set up a 2D channel with:
- `blockMesh` domain 0 < y < 2a, periodic in x
- Boundary conditions: no-slip at y=0 and y=2a, insulating walls (dJ/dn=0)
- Body force in x via `fvOptions` or `buoyantBoussinesqSimpleFoam` analogy
- Run `mhdFoam` to steady state
- Compare to analytic u_x(y) profile

**OpenFOAM installation (Ubuntu):**
```bash
# Check if already installed
mhdFoam -help 2>/dev/null || echo "not installed"

# Install OpenFOAM 10 or 11 from official repo
# See https://openfoam.org/download/ubuntu/
```

**FreeMHD:**
- Repo: typically GitHub or SourceForge
- Usually requires compilation from source; may need additional dependencies
- Attempt only after OpenFOAM is working

---

## 10. Notes for Future Mentor Report

**Do not email yet — notes for later mentor report**

When writing to Sophia / Chris / advisor, the following careful phrasing is appropriate:

- "I implemented an optional default-off reduced Hartmann/wall-drag source term in the standalone 2D OFT/MUG solver (`xmhd_2d.F90`)."
- "The term damps velocity perpendicular to a user-specified applied-field direction, modeling Hartmann-layer friction without resolving the thin wall layer on the mesh."
- "I included the corresponding Jacobian contribution, consistent with the backward-Euler implicit time-stepping used by the existing solver."
- "I added a drag-decay validation test comparing the numerical solution against the backward-Euler analytic prediction. It was built and run on Linux (Ubuntu/gfortran 15.2, MPICH): the perpendicular component decays at the analytic backward-Euler rate (rel. err ~5e-5 at FE orders 2 and 3) and the parallel component is undamped to machine precision. The existing 2D Alfvén and sound-wave tests still pass unchanged with drag default-off (21 passed, 36 skipped overall)."
- "This validates only the reduced drag model. It does not claim incompressible MHD, anisotropic conductivity, current-closure wall BCs, or full LM blanket validity."
- "The larger incompressible LM-MHD formulation described in the TOFE abstract should be developed as a separate module (`xmhd_2d_lm.F90`) rather than added to the compressible solver, to avoid compromising plasma-MHD use cases."
- "Next physics additions in order of risk and value: (1) region-specific material coefficients, (2) anisotropic conductivity tensor, (3) incompressible pressure-velocity formulation, (4) wall current-closure BCs, (5) bidirectional TokaMaker coupling."

---

## 11. Risks and Things to Verify

| Risk | Severity | Status |
|------|----------|--------|
| Sign convention in residual | High | ✅ RESOLVED — drag_decay test passes (velx 3.5e-5) |
| `drag_bhat` not normalized | Medium | ✅ RESOLVED — normalization + abort added (§6 fix 1) |
| Cylindrical velocity component treatment | Medium | Code written; Cartesian path tested. Cyl path NOT yet exercised by a test — still open |
| Jacobian matches residual exactly | High | ✅ RESOLVED — NK solver converges; drag_decay passes at p2 and p3 |
| OMP `private(vdotb_drag)` is correct | Medium | ✅ RESOLVED — ran with 6 OpenMP threads, results correct |
| Default-off leaves existing tests unchanged | High | ✅ RESOLVED — alfven_2d (9) + sound_2d (9) pass against drag build |
| No input parser exposes drag options | Low | Programmatic setting is consistent with existing style |
| `velx` vs `vely` direction in test | Medium | Confirm `vec_vals(1)` is a meaningful representative DOF for uniform field |
| `alpha_drag * dt` must be not too large | Low | `alpha=2, dt=0.05 → alpha*dt=0.1` — should be fine for BE stability |
| Physical mapping `alpha → Ha` undocumented | Low | Document in validation plan before mentor report |
| `mhd_sim%u%get_local` after run_simulation | Medium | Confirmed by reading `run_simulation` line 320: `CALL self%u%add(0,1,u)` at each step |

---

## 12. Final Instruction for Future Claude Code

**Fresh Claude Code: read this entire file before editing. First build and test the
existing Phase 1 wall-drag patch. Do not start Lima or stop any Lima VM automatically.
Do not proceed to OpenFOAM/FreeMHD until OFT builds and the existing/new tests run.
Update this handoff file with every major result.**
