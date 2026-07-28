# 3D xMHD high-Ha solver + transient-overshoot notes (2026-07-28)

Ginsburg build. Companion to `src/tests/physics/test_hartmann3d.F90` (the 3D Hartmann-duct /
disruption driver: `ramp_driver` Sdrv source via the `xmhd_run` `driver=` hook, smooth cosine
turn-on, per-step velocity probe).

## 3D solver: unblock + high-Ha robustness (SOLVED)
The 3D `xmhd` solver was blocked (native solver heap-corrupts; PETSc uncompiled), then walled out at
high Ha (stiff Hartmann layer) with the default DIAGONAL preconditioner.

1. **PETSc build.** `build_libs.py --build_petsc=1 --petsc_superlu_dist=1 --petsc_version=3.23`
   against the bundled MPICH-4.2.3 (the Ginsburg petsc module is OpenMPI = incompatible). superlu_dist
   is REQUIRED — PETSc's built-in LU cannot factor OFT's MPI/block matrices (`MatGetFactor` fails).
2. **Runtime env (order matters):** `LD_LIBRARY_PATH = gcc-11.2/lib64 : mpich/lib : petsc/lib :
   anaconda/lib` (gcc libgfortran before anaconda for GFORTRAN_10; mpich before anaconda for the
   weak `MPI_Comm_c2f`). `use_petsc=T`. Run `mpirun -np 1 ./bin oft.in oft_in.xml`.
3. **Robust preconditioner = the solver XML passed as argv[2]** (`oft_in.xml`). The 3D solver routes
   the PC through multigrid (`create_ml_xml`) even at mesh nlevels=1, so use the MG-smoother format
   with an LU-direct (superlu_dist) smoother+coarse — see `pc_ludirect.xml` in this dir:
   ```xml
   <oft><xmhd><pre>
     <smoother direction="both"><solver type="gmres"><its>2</its><nrits>2</nrits>
       <pre type="lu"></pre></solver></smoother>
     <coarse><solver type="gmres"><its>4</its><nrits>4</nrits>
       <pre type="lu"></pre></solver></coarse>
   </pre></xmhd></oft>
   ```
   **Result: NL converges in 1 iteration, ~1-2 s/step at ANY Ha** (factorization cost ~ mesh size,
   not Ha). `ilu` aborts under PETSc; use `lu` (superlu_dist). **MUST use FE order>=2** (order-1 heap-
   corrupts). packing sets near-wall grading at NO extra factorization cost (matrix size = cell count).

## Transient overshoot: regime map (for reproducing the paper sec:threeD A_os)
Impulsive/ramp-driven Hartmann duct, core v_z(t) overshoot A_os = peak/plateau.
- **Pm = nu/eta matters:** Pm=1 (nu=eta=1) is viscously OVER-damped -> monotonic rise, NO overshoot.
  The standing-Alfven overshoot only appears at LOW Pm (~1e-6, the liquid-metal regime).
- **A_os-1 scales ~linearly with Ha** (matches box-1's amplitude law alpha~1.0): small at low Ha,
  ~1.9 at box-1's Ha=1323.
- **Resolve the Hartmann layer** (delta = a/Ha): under-resolved -> braking missing -> the low-Pm
  oscillation GROWS/blows up instead of settling to a plateau. High packing resolves it cheaply.
- **Finite-L mechanism** (bz=0 axial caps, `ref_per=F,F,F`): caps suppress the short-duct flow ->
  A_os rises with L/a toward the 2.5D (periodic) limit. Demonstrated qualitatively at Ha=18.

## TODO to match box-1's A_os(L/a) = 1.24/1.73/1.91 (PbLi re-calibration)
Switch from demo units (den_scale=5.98e26, rho~2) to box-1 PbLi: **rho=9486 => den_scale~2.84e30**,
eta(=lambda)=1.137 (SI eta=1.4286e-6 Ohm.m), Ha=2645.7*B0y with **B0y=0.5T => Ha=1323**, c=600,
tau_A=2a/U_A=0.0437s. Then: LU-direct solver XML + resolved Hartmann layer (high packing) + SHORT runs
(~2-3 tau_A, first overshoot before instability) + the bz=0-cap L/a={1,2,4} sweep. Box-1 driver =
constant Sdrv=100 T/s in the induction eq (the `ramp_driver`).
