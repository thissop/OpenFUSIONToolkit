# OFT / MUG LM-MHD — Phase 1 Build & Validation Run Summary

**Run date:** 2026-06-02
**Branch:** `lm-mhd-upgrades`
**Repo:** `/Users/tkiker/Documents/GitHub/OpenFUSIONToolkit` (macOS host, mounted into a Lima
Ubuntu VM via virtiofs — this file is on the Mac filesystem)
**Environment where the build/tests ran:** Lima Ubuntu 25.10, aarch64

> This is a self-contained record of the Phase 1 build-and-test run. It does not require
> reading any other file. Companion documents (same `docs/dev/` folder, also on the Mac):
> `lm_mhd_lima_handoff.md` (full handoff), `lm_mhd_upgrade_plan.md`, `lm_mhd_validation_plan.md`,
> `lm_mhd_openfoam_validation.md`.

---

## 1. One-paragraph summary

Phase 1 adds an optional, **default-off** reduced Hartmann / wall-drag source term to the
standalone 2D MUG solver (`src/physics/xmhd_2d.F90`). This run compiled the full Open FUSION
Toolkit from source on Linux, ran the existing 2D MUG regression tests (Alfvén wave, sound
wave) and the new drag-decay test, and confirmed **all pass (21 passed, 36 skipped)**. The
reduced-drag model damps velocity perpendicular to a user-specified field direction at the
correct backward-Euler rate while leaving the parallel component untouched; existing
plasma-MHD results are unchanged when the feature is disabled. Three issues found during the
run were fixed (one source hardening + two test-harness bugs). This validates the **reduced
drag model only** — not full wall-resolved liquid-metal blanket MHD.

---

## 2. What Phase 1 implements

In `src/physics/xmhd_2d.F90`, three new fields on `oft_xmhd_2d_sim` (default-off):

```fortran
LOGICAL  :: use_wall_drag = .FALSE.            ! enable the reduced drag source
REAL(r8) :: drag_coeff    = 0.d0              ! drag coefficient alpha [1/time]
REAL(r8) :: drag_bhat(3)  = [0.d0,1.d0,0.d0]  ! applied-field direction (drag NOT along this)
```

Physics — a volumetric momentum sink modelling unresolved Hartmann-layer friction:

```
S_u = -alpha * u_perp ,     u_perp = u - (u . b_hat) b_hat
```

For pure drag this gives `u_perp(t) = u_perp(0) exp(-alpha t)`, with backward-Euler discrete
form `u_perp_{n+1} = u_perp_n / (1 + alpha*dt)`. The term is added to both the nonlinear
residual (`nlfun_apply`) and the analytic Jacobian (`build_approx_jacobian`), in both the
Cartesian and cylindrical paths. When `use_wall_drag = .FALSE.` (the default), none of the
new code executes.

---

## 3. Build (Linux, from source)

**Toolchain:** gcc/g++/gfortran 15.2.0, CMake 3.31.6, Python 3.13 (venv), OpenMP (6 threads),
native LA backend. External libraries built from source by `build_libs.py`: MPICH 4.2.3,
HDF5 1.14.6, OpenBLAS 0.3.30, METIS 5.1.0, SuperLU 7.0.0, UMFPACK 6.3.5, ARPACK-ng 3.9.1,
libxml2 2.15.2.

**Commands (exact):**
```bash
# Ubuntu prerequisites
sudo apt-get install -y python3-venv python3-dev build-essential m4 autoconf automake libtool

# Python venv + deps
python3 -m venv oft_venv
echo "source $(pwd)/oft_venv/bin/activate" > setup_env.sh
source setup_env.sh
python -m pip install pytest numpy scipy h5py matplotlib xarray

# Stage 1 — external libraries (~20 min; mirrors CI copilot-setup-steps.yml)
mkdir -p builds && cd builds
export CC=gcc CXX=g++ FC=gfortran
python ../src/utilities/build_libs.py --build_umfpack=1 --build_superlu=1 \
  --no_dl_progress --nthread=6 --build_arpack=1 --oft_build_tests=1 --build_mpich=1

# Stage 2 — configure + build OFT
bash config_cmake.sh                 # generates + runs cmake into builds/build_release
cd build_release
#   IMPORTANT: source the venv with an ABSOLUTE path in a fresh shell, e.g.
#   source /Users/tkiker/Documents/GitHub/OpenFUSIONToolkit/oft_venv/bin/activate
make -j6
```

**Result:** OFT built with **zero compile errors**. `xmhd_2d.F90` (with the Phase 1 drag
changes) compiled cleanly; the three test executables built:
`tests/physics/test_alfven_2d`, `test_sound_2d`, `test_drag_decay`.

---

## 4. Tests — ALL PASS

```bash
source /Users/tkiker/Documents/GitHub/OpenFUSIONToolkit/oft_venv/bin/activate
cd builds/build_release/tests/physics
python -m pytest test_alfven_2d.py test_sound_2d.py test_drag_decay.py -v
```

| Test | Purpose | Result |
|------|---------|--------|
| `test_alfven_2d` | Alfvén wave round-trip accuracy (drag default-off) | 9 passed |
| `test_sound_2d`  | Compressible sound wave accuracy (drag default-off) | 9 passed |
| `test_drag_decay`| Backward-Euler drag decay rate (3 cases) | 3 passed |
| **Full suite**   | all three together | **21 passed, 36 skipped** |

(The 36 skips are MPI/coverage parametrizations, not failures.)

### `test_drag_decay` measured values

| Case | FE order | alpha | dt | nsteps | velx rel err (perp, damped) | vely rel err (parallel, undamped) |
|------|----------|-------|----|--------|------------------------------|------------------------------------|
| `test_drag_p2_r1`            | 2 | 2.0 | 0.05 | 10 | 3.52e-5 | 4.44e-16 |
| `test_drag_p2_r1_high_alpha` | 2 | 5.0 | 0.05 | 10 | 3.42e-5 | 4.44e-16 |
| `test_drag_p3_r1`            | 3 | 2.0 | 0.05 | 10 | 5.17e-5 | 5.33e-13 |

Interpretation: the perpendicular component decays at the analytic backward-Euler rate to
~5e-5 (limited by the nonlinear-solver tolerance, not the drag model); the parallel component
is undamped to machine precision — confirming the projection `u - (u.b_hat)b_hat` is correct.
Because `use_wall_drag` defaults to `.FALSE.`, the Alfvén and sound tests are unchanged,
confirming the new path is inert when disabled.

---

## 5. Issues found and fixed during this run

All three were found by actually building and running; the drag residual/Jacobian math was
correct as written. Fixes 2 and 3 are test-harness only.

1. **`drag_bhat` normalization (source hardening).** `drag_bhat` was not normalized, so a
   non-unit vector would not give a true perpendicular projection. Added a one-time
   normalization with a zero-vector abort in `run_simulation` and `run_lin_simulation`
   (`xmhd_2d.F90`, just after `CALL self%setup_bc()`):
   ```fortran
   IF(self%use_wall_drag)THEN
     IF(NORM2(self%drag_bhat)<=1.d-12)THEN
       CALL oft_abort('use_wall_drag=.TRUE. but drag_bhat is a zero vector', &
         'run_simulation',__FILE__)
     END IF
     self%drag_bhat=self%drag_bhat/NORM2(self%drag_bhat)
   END IF
   ```

2. **Segfault in the test (`test_drag_decay.F90`).** `vec_vals` was declared
   `REAL(r8), POINTER :: vec_vals(:)` without `=> NULL()`. `vec_get_local` allocates only
   when `.NOT.ASSOCIATED(array)`, so an undefined pointer association caused a dangling-read
   SIGSEGV in `vec_restore_local` on the first initial-condition assignment. Fixed:
   `REAL(r8), POINTER :: vec_vals(:) => NULL()`.

3. **Order-3 tolerance failure (`test_drag_decay.F90`).** The p3 case failed at 10.25% velx
   error because the adaptive timestep controller (`ittarget=40` default;
   `dt = ittarget*Sum(dt)/Sum(lits)`, capped at the initial dt) shrank dt below 0.05 when the
   linear-iteration count rose at higher FE order — breaking the constant-dt analytic
   prediction. Fixed by pinning dt with `mhd_sim%ittarget = 1000000`, which holds dt at its
   initial value every step. p3 error dropped 10.25% -> 5.17e-5.

---

## 6. Files changed this run (all on the Mac mount; uncommitted)

| File | Status | Change |
|------|--------|--------|
| `src/physics/xmhd_2d.F90` | modified | drag_bhat normalization + zero-vector abort (fix 1) |
| `src/tests/physics/test_drag_decay.F90` | modified | pointer init (fix 2) + ittarget dt-pin (fix 3) |
| `docs/dev/lm_mhd_lima_handoff.md` | modified | build/test results, fixes, resolved-risk table |
| `docs/dev/lm_mhd_validation_plan.md` | modified | drag test marked PASSING + measured values |
| `docs/dev/lm_mhd_openfoam_validation.md` | new (untracked) | OpenFOAM Hartmann cross-validation plan |
| `OFT_MUG_Phase1_Run_Summary.md` | new (this file) | standalone run summary |

> "Uncommitted" = saved to disk on the Mac, not yet committed to git history. To preserve in
> git: `git add -A && git commit` on the `lm-mhd-upgrades` branch.

Build artifacts (`builds/`, `oft_venv/`) were produced inside Linux and are aarch64-Linux
binaries — not usable from macOS natively, and not part of the deliverable. The source and
docs above are plain text and fully portable.

---

## 7. Honest scope — what is and is NOT validated

**Validated:**
- The reduced drag term integrates at the correct backward-Euler decay rate.
- The perpendicular/parallel projection is correct.
- The feature is inert when disabled (no regression in existing 2D MUG tests).
- Residual and Jacobian are consistent (the nonlinear solver converges).

**NOT validated / NOT claimed (future work):**
- Full wall-resolved liquid-metal blanket MHD (the Hartmann boundary layer is modelled, not
  resolved).
- Incompressible MHD (the solver remains compressible).
- Anisotropic / tensor conductivity.
- Wall current-closure / conducting-wall electromagnetic BCs.
- Bidirectional TokaMaker coupling.
- LM duct benchmarks (Hartmann channel flow, Shercliff, Hunt).
- The **cylindrical** drag path is implemented but not yet exercised by a test (the one
  remaining open code item).

---

## 8. Next steps (recommended order)

1. (Optional) Commit the two source files + docs to `lm-mhd-upgrades`.
2. OpenFOAM `mhdFoam` Hartmann cross-check — see `docs/dev/lm_mhd_openfoam_validation.md`
   (OpenFOAM is installable via the Ubuntu `openfoam` apt package; not run this session). Use
   it to calibrate `drag_coeff` to a Hartmann number and check the bulk flow-rate reduction.
3. Add a test for the cylindrical drag path.
4. Only after the above: begin the larger LM-MHD work (region-specific material coefficients,
   then anisotropic conductivity, then a separate incompressible module `xmhd_2d_lm.F90`).

---

*Generated by Claude Code during the 2026-06-02 Lima build/test session. No Lima VM was
started, stopped, or modified. No git commit was made.*
