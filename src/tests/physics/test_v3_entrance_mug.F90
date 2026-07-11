!---------------------------------------------------------------------------
!> @file test_v3_entrance_mug.F90
!
!> V3: developing / entrance MHD channel flow (non-analytic rung, no Phase-4).
!> Formulation per docs/dev/oft_to_comsol.md 2026-07-11 proposal (reduced
!> entrance problem — PIN WITH COMSOL SIDE BEFORE PRODUCTION RUNS):
!>   - in-plane machinery of test_hartmann_mug.F90: streamwise velx(x,z),
!>     wall-normal mesh z, applied field B_0 = (0,0,B0), induced field via psi
!>   - domain x in [0,L], walls z = +/-1; NO periodicity
!>   - inlet x=0:  velx = u_in (slug, = developed core velocity), psi = 0
!>   - outlet x=L: natural (du/dx = dpsi/dx = 0)
!>   - walls: no-slip, insulating (psi natural)
!>   - constant body force fx sustains the developed region; the slug's
!>     Hartmann layers grow over the MHD entry length L_e(Ha)
!> March to steady state; the steady field IS the entrance solution.
!>
!> Output: v3_entrance_mug.profile with "x z velx psi" over mesh points.
!---------------------------------------------------------------------------
PROGRAM test_v3_entrance_mug
USE oft_base
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
USE oft_la_base, ONLY: oft_vector
USE mhd_utils, ONLY: mu0, proton_mass
USE xmhd_2d
IMPLICIT NONE

INTEGER(i4) :: io_unit, ierr, i, np
REAL(r8), POINTER :: vec_vals(:) => NULL(), psi_vals(:) => NULL()

!---Runtime options (namelist)
INTEGER(i4) :: order    = 1
INTEGER(i4) :: nsteps   = 800
INTEGER(i4) :: rst_freq = 100000
INTEGER(i4) :: ittarget = 40
REAL(r8) :: dt          = 1.d-3
REAL(r8) :: a_half      = 1.d0       !< channel half-height (mesh z in [-a,a])
REAL(r8) :: xlen        = 4.d0       !< channel length (mesh x in [0,L])
REAL(r8) :: B0          = 1.121125979562d-1  !< Ha = 100 default
REAL(r8) :: eta         = 1.d0
REAL(r8) :: nu          = 1.d0
REAL(r8) :: fx          = 1.d-2      !< body force per unit mass (streamwise)
REAL(r8) :: n0          = 1.d0
REAL(r8) :: t0          = 1.d0
REAL(r8) :: chi         = 1.d-30
REAL(r8) :: D_diff      = 1.d-30
REAL(r8) :: gamma       = 1.67d0
REAL(r8) :: den_scale   = 5.98d26
REAL(r8) :: lin_tol     = 1.d-9
REAL(r8) :: nl_tol      = 1.d-7
LOGICAL :: pm = .FALSE.
NAMELIST/v3_options/ order, nsteps, rst_freq, ittarget, dt, a_half, xlen, B0, eta, &
  nu, fx, n0, t0, chi, D_diff, gamma, den_scale, lin_tol, nl_tol, pm

TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
CLASS(oft_vector), POINTER :: u => NULL()
REAL(r8) :: Ha_pred, rho, u_in
REAL(r8), PARAMETER :: btol = 1.d-8

CALL oft_init
OPEN(NEWUNIT=io_unit, FILE=oft_env%ifile)
READ(io_unit, v3_options, IOSTAT=ierr)
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
mhd_sim%B_0       = [0.d0, 0.d0, B0]    ! wall-normal applied field (in-plane z)
mhd_sim%use_body_force = .TRUE.
mhd_sim%body_force = [fx, 0.d0, 0.d0]   ! streamwise drive sustains developed region

rho = proton_mass * n0 * den_scale
Ha_pred = B0 * a_half / SQRT(mu0 * eta * rho * nu)
! developed core velocity of the insulating Hartmann channel (high-Ha limit
! G/(sigma*B0^2) with the cosh correction), used as the slug inlet speed
u_in = fx * mu0 * eta * rho / B0**2 * (1.d0 - 1.d0/COSH(Ha_pred))
IF(oft_env%head_proc)THEN
  WRITE(*,'(A,ES12.4)') 'Predicted Ha = ', Ha_pred
  WRITE(*,'(A,ES12.4)') 'Slug inlet u_in = ', u_in
END IF

!---Boundary conditions (coordinate-based masks; order=1 => DOF = vertex)
! n, T, vely, velz frozen everywhere:
ALLOCATE(mhd_sim%n_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%n_bc = .TRUE.
mhd_sim%t_bc => mhd_sim%n_bc
mhd_sim%vely_bc => mhd_sim%n_bc
mhd_sim%velz_bc => mhd_sim%n_bc
! velx: Dirichlet at walls (no-slip) and inlet (slug); FREE at outlet:
ALLOCATE(mhd_sim%velx_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%velx_bc = .FALSE.
! psi: Dirichlet at inlet only (field-free entry); natural at walls (insulating)
! and outlet (outflow):
ALLOCATE(mhd_sim%psi_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%psi_bc = .FALSE.
! by: natural everywhere (identically zero for this in-plane configuration):
ALLOCATE(mhd_sim%by_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%by_bc = .FALSE.
DO i=1,mg_mesh%smesh%np
  IF(ABS(ABS(mg_mesh%smesh%r(2,i)) - a_half) < btol) mhd_sim%velx_bc(i) = .TRUE.  ! walls
  IF(ABS(mg_mesh%smesh%r(1,i)) < btol)THEN                                        ! inlet
    mhd_sim%velx_bc(i) = .TRUE.
    mhd_sim%psi_bc(i) = .TRUE.
  END IF
END DO

!---Initial conditions: slug flow everywhere (walls pinned at 0 by BC init below)
CALL ML_oft_blagrange%vec_create(u)
CALL u%set(n0);    CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,1)
DEALLOCATE(vec_vals); NULLIFY(vec_vals)
! velx = u_in everywhere except exactly 0 at walls (Dirichlet holds initial values)
CALL u%set(u_in);  CALL u%get_local(vec_vals)
DO i=1,mg_mesh%smesh%np
  IF(ABS(ABS(mg_mesh%smesh%r(2,i)) - a_half) < btol) vec_vals(i) = 0.d0
END DO
CALL mhd_sim%u%restore_local(vec_vals,2)
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

!---March to steady state (the steady field is the entrance solution)
CALL mhd_sim%run_simulation()

!---Extract velx and psi over all mesh points
CALL mhd_sim%u%get_local(vec_vals, 2)
CALL mhd_sim%u%get_local(psi_vals, 6)
np = mg_mesh%smesh%np
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit, FILE='v3_entrance_mug.profile')
  WRITE(io_unit,'(A,ES16.8)') '# Ha_pred = ', Ha_pred
  WRITE(io_unit,'(A,ES16.8)') '# fx = ', fx
  WRITE(io_unit,'(A,ES16.8)') '# u_in = ', u_in
  WRITE(io_unit,'(A)') '# x  z  velx  psi'
  DO i=1,np
    WRITE(io_unit,'(4ES20.10)') mg_mesh%smesh%r(1,i), mg_mesh%smesh%r(2,i), &
      vec_vals(i), psi_vals(i)
  END DO
  CLOSE(io_unit)
  WRITE(*,'(A,I0,A)') 'Wrote v3_entrance_mug.profile (', np, ' points)'
END IF
DEALLOCATE(vec_vals, psi_vals); NULLIFY(vec_vals, psi_vals)

CALL oft_finalize()
END PROGRAM test_v3_entrance_mug
