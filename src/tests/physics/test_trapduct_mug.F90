!---------------------------------------------------------------------------
!> @file test_trapduct_mug.F90
!
!> V4a: trapezoidal insulating duct, Ha=100, full induction (non-analytic rung).
!> COMSOL reference: lm-mhd-comsol data/results/v4a_trapduct_Ha100.csv;
!> mesh-converged flow-rate target Q = 2.76810e-2 (forcing -1 nondim units).
!>
!> Same 2.5D machinery as test_shercliff_mug.F90 (axial vely, induced by,
!> B_0 in-plane z, insulating by=0 Dirichlet). The only geometric change: the
!> rectangle mesh is tapered in x by a linear-in-z factor (width base_w at
!> z=-1 to top_w at z=+1) before FEM setup. Compare flow rate Q = int(u)dA
!> and the full (u,b) field, NOT u(0,0) (Ha*u0=1 for any shape).
!---------------------------------------------------------------------------
PROGRAM test_trapduct_mug
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
REAL(r8) :: base_w      = 2.d0       !< duct width at z = -1 (VALIDATION V4a: 2.0)
REAL(r8) :: top_w       = 1.2d0      !< duct width at z = +1 (VALIDATION V4a: 1.2)
LOGICAL :: pm = .FALSE.
NAMELIST/trapduct_options/ order, nsteps, rst_freq, ittarget, dt, a_half, B0, eta, &
  nu, fy, n0, t0, chi, D_diff, gamma, den_scale, lin_tol, nl_tol, base_w, top_w, pm

TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
CLASS(oft_vector), POINTER :: u => NULL()
REAL(r8) :: Ha_pred, rho

CALL oft_init
OPEN(NEWUNIT=io_unit, FILE=oft_env%ifile)
READ(io_unit, trapduct_options, IOSTAT=ierr)
CLOSE(io_unit)
oft_env%pm = pm

CALL multigrid_construct_surf(mg_mesh)
!---Taper the rectangle into the V4a trapezoid: width base_w at z=-1 -> top_w at z=+1.
! Linear-in-z scaling of x; straight triangles + on-the-fly jacobians make the warped
! mesh exact for order-1 elements. Boundary topology/flags are unchanged.
BLOCK
  INTEGER(i4) :: ip
  REAL(r8) :: sfac
  DO ip=1,mg_mesh%smesh%np
    sfac = (base_w + (top_w-base_w)*(mg_mesh%smesh%r(2,ip)+1.d0)/2.d0)/2.d0
    mg_mesh%smesh%r(1,ip) = mg_mesh%smesh%r(1,ip)*sfac
  END DO
END BLOCK
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
! velx, vely, velz, by: leave unset -> setup_bc uses gbe (Dirichlet at walls):
!   no periodicity => all four walls. No-slip velocity; by=0 = insulating
!   (Shercliff). Conducting/thin-wall variants (Hunt H1, Tc) will replace the
!   by Dirichlet flag with the Phase-4 Robin boundary term.

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

!---Extract vely and by over all mesh points (order=1: DOF=vertex)
BLOCK
REAL(r8), POINTER :: by_vals(:)
NULLIFY(by_vals)
CALL mhd_sim%u%get_local(vec_vals, 3)
CALL mhd_sim%u%get_local(by_vals, 7)
np = mg_mesh%smesh%np
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit, FILE='trapduct_mug.profile')
  WRITE(io_unit,'(A,ES16.8)') '# Ha_pred = ', Ha_pred
  WRITE(io_unit,'(A,ES16.8)') '# fy = ', fy
  WRITE(io_unit,'(A)') '# x  z  vely  by   (u_nondim = vely*nu/(fy*a^2))'
  DO i=1,np
    WRITE(io_unit,'(4ES20.10)') mg_mesh%smesh%r(1,i), mg_mesh%smesh%r(2,i), vec_vals(i), by_vals(i)
  END DO
  CLOSE(io_unit)
  WRITE(*,'(A,I0,A)') 'Wrote trapduct_mug.profile (', np, ' points)'
END IF
DEALLOCATE(vec_vals, by_vals); NULLIFY(vec_vals)
END BLOCK

CALL oft_finalize()
END PROGRAM test_trapduct_mug
