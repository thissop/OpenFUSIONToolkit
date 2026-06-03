!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file test_drag_decay.F90
!
!> Regression test for the reduced wall-drag (Hartmann friction) source term
!> added to xmhd_2d for liquid-metal MHD capability.
!>
!> Physics:
!>   With use_wall_drag=T and all other physics nearly inactive, the momentum
!>   equation for a spatially uniform field reduces (per DOF) to:
!>     (1 + dt * drag_coeff) * v_perp_new = v_perp_old
!>   giving backward-Euler discrete decay:
!>     v_perp_{n+1} = v_perp_n / (1 + drag_coeff * dt)
!>
!>   drag_bhat = [0,1,0] means the y-component (toroidal/out-of-plane) is
!>   parallel to the applied field, so:
!>     - velx (x-component) IS damped (perpendicular to bhat)
!>     - vely (y-component) is NOT damped (parallel to bhat)
!>
!> Checks (written to drag_decay.results):
!>   Line 1: relative error of velx vs backward-Euler prediction (expect < 0.02)
!>   Line 2: relative error of vely vs initial value (expect < 1e-6)
!---------------------------------------------------------------------------
PROGRAM test_drag_decay
!---Runtime
USE oft_base
!---Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
!---LA
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!---FE
USE oft_blag_operators, ONLY: oft_blag_getmop
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: elec_charge, proton_mass, mu0
USE xmhd_2d
IMPLICIT NONE

INTEGER(i4) :: io_unit, ierr
REAL(r8), POINTER :: vec_vals(:) => NULL()

!---Runtime options (read from namelist)
INTEGER(i4) :: order    = 2
INTEGER(i4) :: nsteps   = 10
INTEGER(i4) :: rst_freq = 1000
REAL(r8) :: dt          = 0.05d0
REAL(r8) :: n0          = 1.d19
REAL(r8) :: velx0       = 1.d0    !< initial velx (perp to drag_bhat=[0,1,0])
REAL(r8) :: vely0       = 0.75d0  !< initial vely (par to drag_bhat, should not decay)
REAL(r8) :: t0          = 1.d3
REAL(r8) :: chi         = 1.d-30
REAL(r8) :: eta         = 1.d-30
REAL(r8) :: nu          = 1.d-30
REAL(r8) :: D_diff      = 1.d-30
REAL(r8) :: gamma       = 1.67d0
REAL(r8) :: den_scale   = 1.d19
REAL(r8) :: alpha_drag  = 2.d0    !< drag coefficient alpha [1/time]
REAL(r8) :: lin_tol     = 1.d-10
REAL(r8) :: nl_tol      = 1.d-8
LOGICAL :: pm = .FALSE.
NAMELIST/drag_test_options/ order, nsteps, rst_freq, dt, n0, velx0, vely0, t0, &
  chi, eta, nu, D_diff, gamma, den_scale, alpha_drag, lin_tol, nl_tol, pm

!---Simulation objects
TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
CLASS(oft_vector), POINTER :: u => NULL()
REAL(r8) :: velx_final, vely_final
REAL(r8) :: velx_expected, velx_relerr, vely_relerr

!==============================================================================
! Initialize and read options
!==============================================================================
CALL oft_init
OPEN(NEWUNIT=io_unit, FILE=oft_env%ifile)
READ(io_unit, drag_test_options, IOSTAT=ierr)
CLOSE(io_unit)

!==============================================================================
! Set up mesh and FE space
!==============================================================================
CALL multigrid_construct_surf(mg_mesh)
CALL mhd_sim%setup(mg_mesh, order)

!==============================================================================
! Configure physics
!==============================================================================
mhd_sim%dt        = dt
mhd_sim%nsteps    = nsteps
mhd_sim%rst_freq  = rst_freq
mhd_sim%chi       = chi
mhd_sim%eta       = eta
mhd_sim%nu        = nu
mhd_sim%D_diff    = D_diff
mhd_sim%gamma     = gamma
mhd_sim%den_scale = den_scale
mhd_sim%lin_tol   = lin_tol
mhd_sim%nl_tol    = nl_tol
mhd_sim%mfnk      = .FALSE.
!---Pin the timestep: a large ittarget makes the adaptive controller cap dt at
!   its initial value every step (see run_simulation: dt = ittarget*Sum(dt)/Sum(lits),
!   then capped at dtin). A fixed dt is required for the analytic backward-Euler
!   decay prediction velx0/(1+alpha*dt)^nsteps to be exact, independent of FE order.
mhd_sim%ittarget  = 1000000

!---Enable wall drag perpendicular to y (bhat=[0,1,0])
!   velx and velz are damped; vely is not.
mhd_sim%use_wall_drag = .TRUE.
mhd_sim%drag_coeff    = alpha_drag
mhd_sim%drag_bhat     = [0.d0, 1.d0, 0.d0]

!==============================================================================
! Set boundary conditions BEFORE run_simulation (which calls setup_bc internally).
! setup_bc skips any pointer that is already associated.
!
! Velocity: natural (free-slip) at all boundaries → .FALSE. everywhere.
!   This avoids Dirichlet-zero walls conflicting with the uniform-decay solution.
! n, T: Dirichlet-fixed at all nodes → .TRUE. everywhere.
!   This prevents n and T from evolving and coupling back to momentum.
! psi, by: use default (Dirichlet zero at boundary from gbe) — fine since B≈0.
!==============================================================================
ALLOCATE(mhd_sim%velx_bc(ML_oft_blagrange%current_level%ne))
ALLOCATE(mhd_sim%vely_bc(ML_oft_blagrange%current_level%ne))
ALLOCATE(mhd_sim%velz_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%velx_bc = .FALSE.
mhd_sim%vely_bc = .FALSE.
mhd_sim%velz_bc = .FALSE.

ALLOCATE(mhd_sim%n_bc(ML_oft_blagrange%current_level%ne))
ALLOCATE(mhd_sim%T_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%n_bc = .TRUE.   ! Dirichlet: fix n = n0 everywhere
mhd_sim%T_bc = .TRUE.   ! Dirichlet: fix T = t0 everywhere

!==============================================================================
! Set initial conditions (uniform fields)
!==============================================================================
CALL ML_oft_blagrange%vec_create(u)

!---n = n0/den_scale (uniform)
CALL u%set(n0 / den_scale)
CALL u%get_local(vec_vals)
CALL mhd_sim%u%restore_local(vec_vals, 1)
DEALLOCATE(vec_vals)
NULLIFY(vec_vals)

!---velx = velx0 (uniform, perpendicular to drag_bhat → will decay)
CALL u%set(velx0)
CALL u%get_local(vec_vals)
CALL mhd_sim%u%restore_local(vec_vals, 2)
DEALLOCATE(vec_vals)
NULLIFY(vec_vals)

!---vely = vely0 (uniform, parallel to drag_bhat → should NOT decay)
CALL u%set(vely0)
CALL u%get_local(vec_vals)
CALL mhd_sim%u%restore_local(vec_vals, 3)
DEALLOCATE(vec_vals)
NULLIFY(vec_vals)

!---velz = 0 (uniform)
CALL u%set(0.d0)
CALL u%get_local(vec_vals)
CALL mhd_sim%u%restore_local(vec_vals, 4)
DEALLOCATE(vec_vals)
NULLIFY(vec_vals)

!---T = t0 (uniform)
CALL u%set(t0)
CALL u%get_local(vec_vals)
CALL mhd_sim%u%restore_local(vec_vals, 5)
DEALLOCATE(vec_vals)
NULLIFY(vec_vals)

!---psi = 0, by = 0 (default zero)

!==============================================================================
! Run simulation
! After completion, mhd_sim%u holds the final state (updated each step at line
! ~320 of run_simulation: CALL self%u%add(0.d0,1.d0,u)).
!==============================================================================
CALL mhd_sim%run_simulation()

!==============================================================================
! Read final velocities from mhd_sim%u
! For a uniform initial condition with drag only, all DOF values are equal.
! We read node/edge index 1 as representative.
!==============================================================================
CALL mhd_sim%u%get_local(vec_vals, 2)   ! velx field
velx_final = vec_vals(1)
DEALLOCATE(vec_vals)
NULLIFY(vec_vals)

CALL mhd_sim%u%get_local(vec_vals, 3)   ! vely field
vely_final = vec_vals(1)
DEALLOCATE(vec_vals)
NULLIFY(vec_vals)

!==============================================================================
! Backward-Euler analytic prediction after nsteps:
!   v_perp_{n+1} = v_perp_n / (1 + alpha * dt)
!   => after nsteps: velx = velx0 / (1 + alpha*dt)^nsteps
!==============================================================================
velx_expected = velx0 / (1.d0 + alpha_drag*dt)**nsteps
velx_relerr   = ABS(velx_final - velx_expected) / ABS(velx_expected)
vely_relerr   = ABS(vely_final - vely0) / ABS(vely0)

!==============================================================================
! Write results for pytest
!==============================================================================
OPEN(NEWUNIT=io_unit, FILE='drag_decay.results')
WRITE(io_unit, *) velx_relerr
WRITE(io_unit, *) vely_relerr
CLOSE(io_unit)

IF(velx_relerr < 0.02d0 .AND. vely_relerr < 1.d-6) THEN
  WRITE(*,'(A)') 'PASS: drag decay test'
ELSE
  WRITE(*,'(A,2(A,ES12.4))') 'FAIL:', '  velx_relerr=', velx_relerr, &
    '  vely_relerr=', vely_relerr
  CALL oft_abort('Drag decay test failed', 'test_drag_decay', __FILE__)
END IF

CALL oft_finalize()

END PROGRAM test_drag_decay
