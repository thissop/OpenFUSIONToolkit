!---------------------------------------------------------------------------
!> @file test_hartmann_mug.F90
!
!> Driven Hartmann channel in the resolved 2D MUG solver (xmhd_2d), drag OFF.
!>
!> PERSONAL VERIFICATION exercise (not a MUG capability claim): reproduce the
!> classic insulating-wall Hartmann velocity profile by resolving the Hartmann
!> layers directly, using the default-off uniform body-force momentum source to
!> drive a fully-developed channel.
!>
!> Geometry mapping (see xmhd_2d.F90: mesh plane is x-z, y is the out-of-plane
!> invariant direction; basis_grads(2)=0, div_vel=dvx/dx+dvz/dz):
!>   - streamwise flow  : velx, function of the wall-normal coordinate
!>   - wall-normal coord : mesh second coordinate (smesh%r(2,:)), walls at +/-a
!>   - applied field     : B_0 = (0,0,B0)  (wall-normal, in-plane z-component)
!>   - induced field     : b_x = -dpsi/d(wall-normal); Lorentz braking on velx
!> Hartmann number:  Ha = B0 * a / sqrt(mu0 * eta * rho * nu),  rho = m_i * n.
!>
!> Setup:
!>   - periodic in x (oft.in: ref_per=T,F,F), no-slip walls (default gbe BC)
!>   - n, T frozen (Dirichlet everywhere) -> incompressible, no pressure force
!>   - psi: natural BC (insulating walls: dpsi/dn=0 => b_x=0 at walls)
!>   - body_force = (fx,0,0) drives the flow; march to steady state
!>
!> Output: hartmann_mug.profile with rows "x  y  velx" over all mesh points
!> (order=1 => DOFs at vertices). A Python script bins velx by y and compares
!> to the analytic Hartmann profile (fits Ha) — see cases/hartmann/mug/.
!---------------------------------------------------------------------------
PROGRAM test_hartmann_mug
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
REAL(r8) :: a_half      = 1.d0       !< channel half-width (domain must match in oft.in)
REAL(r8) :: B0          = 5.6d-3     !< applied wall-normal field; Ha=B0*a/sqrt(mu0*eta*rho*nu)
REAL(r8) :: eta         = 1.d0       !< magnetic diffusivity (resistivity in mu0 units)
REAL(r8) :: nu          = 1.d0       !< kinematic viscosity
REAL(r8) :: fx          = 1.d-2      !< body force per unit mass (streamwise)
REAL(r8) :: n0          = 1.d0       !< density field value (internal n = n0*den_scale)
REAL(r8) :: t0          = 1.d0       !< temperature (frozen)
REAL(r8) :: chi         = 1.d-30
REAL(r8) :: D_diff      = 1.d-30
REAL(r8) :: gamma       = 1.67d0
REAL(r8) :: den_scale   = 1.d0       !< set so rho = m_i*n0*den_scale is as intended
REAL(r8) :: lin_tol     = 1.d-9
REAL(r8) :: nl_tol      = 1.d-7
LOGICAL :: pm = .FALSE.
NAMELIST/hartmann_options/ order, nsteps, rst_freq, ittarget, dt, a_half, B0, eta, &
  nu, fx, n0, t0, chi, D_diff, gamma, den_scale, lin_tol, nl_tol, pm

TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
CLASS(oft_vector), POINTER :: u => NULL()
REAL(r8) :: Ha_pred, rho

CALL oft_init
OPEN(NEWUNIT=io_unit, FILE=oft_env%ifile)
READ(io_unit, hartmann_options, IOSTAT=ierr)
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
mhd_sim%body_force = [fx, 0.d0, 0.d0]   ! streamwise drive
! drag stays OFF (default use_wall_drag=.FALSE.)

rho = proton_mass * n0 * den_scale
Ha_pred = B0 * a_half / SQRT(mu0 * eta * rho * nu)
IF(oft_env%head_proc) WRITE(*,'(A,ES12.4)') 'Predicted Ha = ', Ha_pred

!---Boundary conditions
! n, T frozen everywhere (incompressible, no pressure force):
ALLOCATE(mhd_sim%n_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%n_bc = .TRUE.
mhd_sim%t_bc => mhd_sim%n_bc
! psi: natural BC at walls (insulating: dpsi/dn=0 => induced b_x=0 at walls):
ALLOCATE(mhd_sim%psi_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%psi_bc = .FALSE.
! velx, vely, velz, by: leave unset -> setup_bc uses gbe (Dirichlet at walls):
!   x is periodic so gbe = wall nodes only => no-slip walls. by=0 (Dirichlet).

!---Initial conditions: fluid at rest, uniform n and T, no field perturbation
CALL ML_oft_blagrange%vec_create(u)
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

!---March to steady state
CALL mhd_sim%run_simulation()

!---Extract streamwise velocity profile velx over all mesh points (order=1: DOF=vertex)
CALL mhd_sim%u%get_local(vec_vals, 2)
np = mg_mesh%smesh%np
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit, FILE='hartmann_mug.profile')
  WRITE(io_unit,'(A,ES16.8)') '# Ha_pred = ', Ha_pred
  WRITE(io_unit,'(A)') '# x  y  velx'
  DO i=1,np
    WRITE(io_unit,'(3ES20.10)') mg_mesh%smesh%r(1,i), mg_mesh%smesh%r(2,i), vec_vals(i)
  END DO
  CLOSE(io_unit)
  WRITE(*,'(A,I0,A)') 'Wrote hartmann_mug.profile (', np, ' points)'
END IF
DEALLOCATE(vec_vals); NULLIFY(vec_vals)

CALL oft_finalize()
END PROGRAM test_hartmann_mug
