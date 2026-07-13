!---------------------------------------------------------------------------
!> @file test_cavity_mug.F90
!
!> Lid-driven cavity — DECISIVE go/no-go for whether MUG's true-pressure
!> low-Mach mode can produce incompressible RECIRCULATION (the capability V4
!> BFS needs; see docs/dev/oft_to_comsol.md V4 correction). If a primary vortex
!> forms and the center matches the Ghia et al. (1982) benchmark, V4 is
!> feasible; if the low-Mach flow won't recirculate/stabilize, V4 is off.
!>
!> Config: unit square, top lid (z = a) moves at velx = U_lid, other 3 walls
!> no-slip. Unfrozen n (=> grad(nT) pressure force active, low-Mach ~
!> incompressible); n pinned at one corner as the pressure reference; T frozen;
!> B0 = 0 (pure hydrodynamic recirculation test first — Ha added later).
!> In-plane flow velx (streamwise, mesh r1) + velz (wall-normal, mesh r2).
!> Re = U_lid * a / nu.
!>
!> Output: cavity_mug.profile "x z velx velz" — post-process for the vortex
!> center (velx=velz=0 interior point) and compare to Ghia Re=100/400.
!---------------------------------------------------------------------------
PROGRAM test_cavity_mug
USE oft_base
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
USE oft_la_base, ONLY: oft_vector
USE mhd_utils, ONLY: mu0, proton_mass
USE xmhd_2d
IMPLICIT NONE

INTEGER(i4) :: io_unit, ierr, i, np
REAL(r8), POINTER :: vec_vals(:) => NULL(), vz_vals(:) => NULL()

INTEGER(i4) :: order    = 1
INTEGER(i4) :: nsteps   = 3000
INTEGER(i4) :: rst_freq = 100000
INTEGER(i4) :: ittarget = 30
REAL(r8) :: dt          = 1.d-5      !< acoustic-CFL limited (low-Mach); dt-cap => this is the max
REAL(r8) :: a_len       = 1.d0       !< cavity side
REAL(r8) :: u_lid       = 1.d-2      !< lid speed (small => low Mach)
REAL(r8) :: nu          = 1.d-4      !< Re = u_lid*a/nu = 100 with these defaults
REAL(r8) :: eta         = 1.d0
REAL(r8) :: B0          = 0.d0       !< Ha=0 (pure hydro recirculation test)
REAL(r8) :: n0          = 1.d0
REAL(r8) :: t0          = 1.d0       !< sets sound speed; low-Mach => u_lid << sqrt(T)
REAL(r8) :: chi         = 1.d-30
REAL(r8) :: D_diff      = 1.d-6
REAL(r8) :: gamma       = 1.67d0
REAL(r8) :: den_scale   = 5.98d26
REAL(r8) :: lin_tol     = 1.d-9
REAL(r8) :: nl_tol      = 1.d-7
LOGICAL :: pm = .FALSE.
NAMELIST/cavity_options/ order, nsteps, rst_freq, ittarget, dt, a_len, u_lid, nu, &
  eta, B0, n0, t0, chi, D_diff, gamma, den_scale, lin_tol, nl_tol, pm

TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
CLASS(oft_vector), POINTER :: u => NULL()
REAL(r8) :: rho
REAL(r8), PARAMETER :: btol = 1.d-8

CALL oft_init
OPEN(NEWUNIT=io_unit, FILE=oft_env%ifile)
READ(io_unit, cavity_options, IOSTAT=ierr)
CLOSE(io_unit)
oft_env%pm = pm

CALL multigrid_construct_surf(mg_mesh)
CALL mhd_sim%setup(mg_mesh, order)

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
mhd_sim%B_0       = [0.d0, 0.d0, B0]

rho = proton_mass * n0 * den_scale
IF(oft_env%head_proc) WRITE(*,'(A,ES12.4)') 'Re = u_lid*a/nu = ', u_lid*a_len/nu

!---BCs: true-pressure mode. n free except pinned at the (0,0) corner (pressure
! reference); T frozen; velx/velz Dirichlet on all walls (lid moves, rest no-slip);
! psi/by natural (B0=0 => no induction dynamics of interest).
ALLOCATE(mhd_sim%n_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%n_bc = .FALSE.
ALLOCATE(mhd_sim%t_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%t_bc = .TRUE.
ALLOCATE(mhd_sim%velx_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%velx_bc = .FALSE.
ALLOCATE(mhd_sim%velz_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%velz_bc = .FALSE.
mhd_sim%vely_bc => mhd_sim%t_bc         ! vely (out-of-plane) frozen 0
ALLOCATE(mhd_sim%psi_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%psi_bc = .TRUE.                 ! psi pinned (no in-plane induction at B0=0)
ALLOCATE(mhd_sim%by_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%by_bc = .TRUE.
DO i=1,mg_mesh%smesh%np
  ! all 4 walls: velx, velz Dirichlet (values set in ICs below)
  IF(ABS(mg_mesh%smesh%r(1,i)) < btol .OR. ABS(mg_mesh%smesh%r(1,i)-a_len) < btol .OR. &
     ABS(mg_mesh%smesh%r(2,i)) < btol .OR. ABS(mg_mesh%smesh%r(2,i)-a_len) < btol)THEN
    mhd_sim%velx_bc(i) = .TRUE.
    mhd_sim%velz_bc(i) = .TRUE.
  END IF
  ! pressure reference: pin n at the (0,0) corner
  IF(ABS(mg_mesh%smesh%r(1,i)) < btol .AND. ABS(mg_mesh%smesh%r(2,i)) < btol) &
    mhd_sim%n_bc(i) = .TRUE.
END DO

!---ICs: fluid at rest, uniform n,T; lid velx = u_lid on the top wall (z=a)
CALL ML_oft_blagrange%vec_create(u)
CALL u%set(n0);  CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,1)
DEALLOCATE(vec_vals); NULLIFY(vec_vals)
CALL u%set(0.d0);  CALL u%get_local(vec_vals)
DO i=1,mg_mesh%smesh%np
  IF(ABS(mg_mesh%smesh%r(2,i)-a_len) < btol) vec_vals(i) = u_lid   ! moving lid
END DO
CALL mhd_sim%u%restore_local(vec_vals,2)
DEALLOCATE(vec_vals); NULLIFY(vec_vals)
CALL u%set(0.d0);  CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,3)
DEALLOCATE(vec_vals); NULLIFY(vec_vals)
CALL u%set(0.d0);  CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,4)
DEALLOCATE(vec_vals); NULLIFY(vec_vals)
CALL u%set(t0);  CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,5)
DEALLOCATE(vec_vals); NULLIFY(vec_vals)
CALL u%set(0.d0);  CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,6)
DEALLOCATE(vec_vals); NULLIFY(vec_vals)
CALL u%set(0.d0);  CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,7)
DEALLOCATE(vec_vals); NULLIFY(vec_vals)

CALL mhd_sim%run_simulation()

CALL mhd_sim%u%get_local(vec_vals, 2)
CALL mhd_sim%u%get_local(vz_vals, 4)
np = mg_mesh%smesh%np
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit, FILE='cavity_mug.profile')
  WRITE(io_unit,'(A,ES16.8)') '# u_lid = ', u_lid
  WRITE(io_unit,'(A)') '# x  z  velx  velz'
  DO i=1,np
    WRITE(io_unit,'(4ES20.10)') mg_mesh%smesh%r(1,i), mg_mesh%smesh%r(2,i), &
      vec_vals(i), vz_vals(i)
  END DO
  CLOSE(io_unit)
  WRITE(*,'(A,I0,A)') 'Wrote cavity_mug.profile (', np, ' points)'
END IF
DEALLOCATE(vec_vals, vz_vals); NULLIFY(vec_vals, vz_vals)

CALL oft_finalize()
END PROGRAM test_cavity_mug
