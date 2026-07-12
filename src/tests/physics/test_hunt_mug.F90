!---------------------------------------------------------------------------
!> @file test_hunt_mug.F90
!
!> Hunt duct (VALIDATION_CASES.md case H1): conducting Hartmann walls,
!> insulating side walls, out-of-plane flow — same 2.5D machinery as
!> test_shercliff_mug.F90, only the by wall BC changes:
!>   - side walls x=+/-a: insulating, Dirichlet by=0
!>   - Hartmann walls z=+/-a: c_wall=0 -> natural d(by)/dn=0 (perfect
!>     conductor, Hunt 1965); c_wall>0 -> Phase-4 thin-wall Robin BC
!>     by + c_wall*d(by)/dn = 0 (finite wall conductance)
!> Signature physics vs Shercliff: high-velocity side-wall jets.
!> Reference: cases/hunt/hunt_duct.py FD referee (perfect-conductor case).
!---------------------------------------------------------------------------
PROGRAM test_hunt_mug
USE oft_base
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
USE oft_la_base, ONLY: oft_vector
USE mhd_utils, ONLY: mu0, proton_mass
USE xmhd_2d
IMPLICIT NONE

INTEGER(i4) :: io_unit, ierr, i, np
REAL(r8), POINTER :: vec_vals(:) => NULL()

!---Runtime options (namelist)
INTEGER(i4) :: order    = 1
INTEGER(i4) :: nsteps   = 400
INTEGER(i4) :: rst_freq = 100000
INTEGER(i4) :: ittarget = 1000000   !< large => pin dt at its initial value
REAL(r8) :: dt          = 2.d-2
REAL(r8) :: a_half      = 1.d0       !< duct half-width (domain must match in oft.in)
REAL(r8) :: B0          = 5.6d-2     !< applied field (in-plane z); Ha=B0*a/sqrt(mu0*eta*rho*nu)
REAL(r8) :: eta         = 1.d0       !< magnetic diffusivity (resistivity in mu0 units)
REAL(r8) :: nu          = 1.d0       !< kinematic viscosity
REAL(r8) :: fy          = 1.d-2      !< body force per unit mass (axial, out-of-plane)
REAL(r8) :: n0          = 1.d0       !< density field value (internal n = n0*den_scale)
REAL(r8) :: t0          = 1.d0       !< temperature (frozen)
REAL(r8) :: chi         = 1.d-30
REAL(r8) :: D_diff      = 1.d-30
REAL(r8) :: gamma       = 1.67d0
REAL(r8) :: den_scale   = 1.d0       !< set so rho = m_i*n0*den_scale is as intended
REAL(r8) :: lin_tol     = 1.d-9
REAL(r8) :: nl_tol      = 1.d-7
REAL(r8) :: c_wall      = 0.d0       !< Hartmann-wall conductance ratio: 0 = perfect
                                     !  conductor (natural d(by)/dn=0, Hunt's case);
                                     !  >0 = thin-wall Robin (Phase-4 finite c)
CHARACTER(LEN=256) :: restart_file = '' !< if set, continue from this .rst instead of
                                        !  building ICs from scratch (restart-chaining
                                        !  for slow-settling high-Ha conjugate cases; the
                                        !  vector, t, and dt are all restored from the file)
INTEGER(i4) :: rst_base = 0             !< output restart index offset (continue numbering)
LOGICAL :: pm = .FALSE.
NAMELIST/hunt_options/ order, nsteps, rst_freq, ittarget, dt, a_half, B0, eta, &
  nu, fy, n0, t0, chi, D_diff, gamma, den_scale, lin_tol, nl_tol, c_wall, &
  restart_file, rst_base, pm

TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
CLASS(oft_vector), POINTER :: u => NULL()
REAL(r8) :: Ha_pred, rho

CALL oft_init
OPEN(NEWUNIT=io_unit, FILE=oft_env%ifile)
READ(io_unit, hunt_options, IOSTAT=ierr)
CLOSE(io_unit)
oft_env%pm = pm

CALL multigrid_construct_surf(mg_mesh)
CALL mhd_sim%setup(mg_mesh, order)

!---Physics
mhd_sim%dt        = dt
mhd_sim%nsteps    = nsteps
mhd_sim%rst_freq  = rst_freq
mhd_sim%ittarget  = ittarget
mhd_sim%chi       = chi
mhd_sim%eta       = eta
mhd_sim%nu        = nu
mhd_sim%D_diff    = D_diff
mhd_sim%gamma     = gamma
mhd_sim%den_scale = den_scale
mhd_sim%lin_tol   = lin_tol
mhd_sim%nl_tol    = nl_tol
mhd_sim%mfnk      = .FALSE.
mhd_sim%B_0       = [0.d0, 0.d0, B0]    ! Hartmann direction (in-plane z)
mhd_sim%use_body_force = .TRUE.
mhd_sim%body_force = [0.d0, fy, 0.d0]   ! axial (out-of-plane) drive
! drag stays OFF (default use_wall_drag=.FALSE.)

rho = proton_mass * n0 * den_scale
Ha_pred = B0 * a_half / SQRT(mu0 * eta * rho * nu)
IF(oft_env%head_proc) WRITE(*,'(A,ES12.4)') 'Predicted Ha = ', Ha_pred

!---Boundary conditions
! n, T frozen everywhere (incompressible, no pressure force):
ALLOCATE(mhd_sim%n_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%n_bc = .TRUE.
mhd_sim%t_bc => mhd_sim%n_bc
! psi: natural BC (in-plane induced field; stays zero for Shercliff flow):
ALLOCATE(mhd_sim%psi_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%psi_bc = .FALSE.
! velx, vely, velz: leave unset -> setup_bc uses gbe (no-slip Dirichlet, all walls).
! by (Hunt's case): Dirichlet by=0 ONLY on the insulating SIDE walls x=+/-a
! (corners included); FREE on the Hartmann walls z=+/-a. With c_wall=0 the free
! Hartmann walls take the natural condition d(by)/dn=0 = PERFECT conductor
! (Hunt 1965); with c_wall>0 the Phase-4 thin-wall Robin term applies there.
mhd_sim%c_wall = c_wall
ALLOCATE(mhd_sim%by_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%by_bc = .FALSE.
BLOCK
INTEGER(i4) :: ip
DO ip=1,mg_mesh%smesh%np
  IF(ABS(ABS(mg_mesh%smesh%r(1,ip)) - a_half) < 1.d-8) mhd_sim%by_bc(ip) = .TRUE.
END DO
END BLOCK

!---Initial conditions
CALL ML_oft_blagrange%vec_create(u)
mhd_sim%rst_base = rst_base
IF(LEN_TRIM(restart_file) > 0)THEN
  !---Restart-chaining: continue from a saved state (vector + t + dt restored).
  ! Enables slow-settling high-Ha conjugate cases to accumulate toward steady
  ! across multiple walltime-limited jobs, and captures transient waveforms.
  CALL mhd_sim%rst_load(mhd_sim%u, TRIM(restart_file), 'U', mhd_sim%t, mhd_sim%dt)
  IF(oft_env%head_proc) WRITE(*,'(A,A,A,ES12.4,A,ES12.4)') &
    'Restarted from ', TRIM(restart_file), ' at t=', mhd_sim%t, ' dt=', mhd_sim%dt
ELSE
  !---Fresh start: fluid at rest, uniform n and T, no field perturbation
  CALL u%set(n0);    CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,1)
  DEALLOCATE(vec_vals); NULLIFY(vec_vals)
  CALL u%set(0.d0);  CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,2)
  DEALLOCATE(vec_vals); NULLIFY(vec_vals)
  CALL u%set(0.d0);  CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,3)
  DEALLOCATE(vec_vals); NULLIFY(vec_vals)
  CALL u%set(0.d0);  CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,4)
  DEALLOCATE(vec_vals); NULLIFY(vec_vals)
  CALL u%set(t0);    CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,5)
  DEALLOCATE(vec_vals); NULLIFY(vec_vals)
  CALL u%set(0.d0);  CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,6)
  DEALLOCATE(vec_vals); NULLIFY(vec_vals)
  CALL u%set(0.d0);  CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,7)
  DEALLOCATE(vec_vals); NULLIFY(vec_vals)
END IF

!---March (from rest or from the restart state)
CALL mhd_sim%run_simulation()

!---Extract vely and by over all mesh points (order=1: DOF=vertex)
BLOCK
REAL(r8), POINTER :: by_vals(:)
NULLIFY(by_vals)
CALL mhd_sim%u%get_local(vec_vals, 3)
CALL mhd_sim%u%get_local(by_vals, 7)
np = mg_mesh%smesh%np
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit, FILE='hunt_mug.profile')
  WRITE(io_unit,'(A,ES16.8)') '# Ha_pred = ', Ha_pred
  WRITE(io_unit,'(A,ES16.8)') '# fy = ', fy
  WRITE(io_unit,'(A,ES16.8)') '# c_wall = ', c_wall
  WRITE(io_unit,'(A)') '# x  z  vely  by   (u_nondim = vely*nu/(fy*a^2))'
  DO i=1,np
    WRITE(io_unit,'(4ES20.10)') mg_mesh%smesh%r(1,i), mg_mesh%smesh%r(2,i), vec_vals(i), by_vals(i)
  END DO
  CLOSE(io_unit)
  WRITE(*,'(A,I0,A)') 'Wrote hunt_mug.profile (', np, ' points)'
END IF
DEALLOCATE(vec_vals, by_vals); NULLIFY(vec_vals, by_vals)
END BLOCK

CALL oft_finalize()
END PROGRAM test_hunt_mug
