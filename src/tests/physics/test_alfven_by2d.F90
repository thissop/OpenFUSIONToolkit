!---------------------------------------------------------------------------
!> @file test_alfven_by2d.F90
!
!> Shear-Alfven wave in the by<->vely polarization — sign/consistency test for
!> the psi-generated field-line stretching term in the by induction equation.
!>
!> Motivation (see docs/dev/oft_to_comsol.md, 2026-07-10 "MUG rung 1" entry):
!> the background-field (B_0) stretch term's sign had to be fixed empirically
!> for the Shercliff duct; the PRE-EXISTING dpsi-stretch term in the by
!> equation (`cross_product(dpsi, dvel(2,:))`) is exercised by NO other test —
!> test_alfven_2d covers only the velx<->psi polarization. This test closes
!> that gap: if the dpsi term is sign-consistent with the momentum tension
!> term, a vely perturbation on a psi-generated in-plane field oscillates at
!> omega = k*v_A (decaying slowly via backward-Euler dissipation); if it is
!> inconsistent, kinetic energy grows exponentially.
!>
!> Setup:
!>   - background field FROM PSI (deliberately, B_0 = 0): psi = -B*x gives a
!>     uniform in-plane field along mesh z; psi pinned by default Dirichlet BC
!>   - perturbation: vely = v_delta*sin(2*pi*z/lam), by = 0
!>   - periodic in z (oft.in ref_per=F,T,F); vely/by natural BCs (x-uniform);
!>     velx/velz/n/T frozen Dirichlet everywhere
!>   - eta = nu = 1e-6 (near-ideal); expect period T = lam/v_alf
!>
!> PASS: kinetic energy (column 3 of xmhd_2d.hist diag block / printed KE)
!> oscillates at 2*omega and decays; FAIL: monotonic exponential growth.
!---------------------------------------------------------------------------
MODULE alfven_by2d_init
USE oft_base
IMPLICIT NONE
REAL(r8) :: init_v_delta = 1.d-4  !< velocity perturbation amplitude
REAL(r8) :: init_lam = 1.d0       !< wavelength along z
REAL(r8) :: init_B_bg = 0.d0      !< background field magnitude
CONTAINS
SUBROUTINE vely_wave(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = init_v_delta * SIN(2.d0*pi*pt(2)/init_lam)
END SUBROUTINE vely_wave
SUBROUTINE psi_bg(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = -init_B_bg * pt(1)
END SUBROUTINE psi_bg
END MODULE alfven_by2d_init

PROGRAM test_alfven_by2d
USE oft_base
USE alfven_by2d_init
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_blag_operators, ONLY: oft_blag_getmop, oft_blag_project
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: mu0, proton_mass
USE xmhd_2d
IMPLICIT NONE

INTEGER(i4) :: io_unit, ierr
REAL(r8), POINTER :: vec_vals(:) => NULL()
TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
TYPE(poss_scalar_bfield) :: field_init
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u => NULL(), v => NULL()
REAL(r8) :: B_bg, rho

!---Runtime options (namelist)
INTEGER(i4) :: order    = 1
INTEGER(i4) :: nsteps   = 300
INTEGER(i4) :: rst_freq = 100000
INTEGER(i4) :: ittarget = 1000000
REAL(r8) :: dt          = 1.d-3
REAL(r8) :: v_alf       = 10.d0     !< Alfven speed; B_bg = v_alf*sqrt(mu0*rho)
REAL(r8) :: lam         = 1.d0      !< wavelength along z
REAL(r8) :: v_delta     = 1.d-4     !< velocity perturbation amplitude
REAL(r8) :: eta         = 1.d-6
REAL(r8) :: nu          = 1.d-6
REAL(r8) :: n0          = 1.d0
REAL(r8) :: t0          = 1.d0
REAL(r8) :: chi         = 1.d-30
REAL(r8) :: D_diff      = 1.d-30
REAL(r8) :: gamma       = 1.67d0
REAL(r8) :: den_scale   = 5.98d26
REAL(r8) :: lin_tol     = 1.d-10
REAL(r8) :: nl_tol      = 1.d-8
LOGICAL :: pm = .FALSE.
NAMELIST/alfven_by_options/ order, nsteps, rst_freq, ittarget, dt, v_alf, lam, &
  v_delta, eta, nu, n0, t0, chi, D_diff, gamma, den_scale, lin_tol, nl_tol, pm

CALL oft_init
OPEN(NEWUNIT=io_unit, FILE=oft_env%ifile)
READ(io_unit, alfven_by_options, IOSTAT=ierr)
CLOSE(io_unit)
oft_env%pm = pm

CALL multigrid_construct_surf(mg_mesh)
CALL mhd_sim%setup(mg_mesh, order)

rho = proton_mass * n0 * den_scale
B_bg = v_alf * SQRT(mu0 * rho)
init_v_delta = v_delta
init_lam = lam
init_B_bg = B_bg

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
mhd_sim%B_0       = 0.d0            ! background comes from psi ON PURPOSE

IF(oft_env%head_proc)THEN
  WRITE(*,'(A,ES12.4)') 'B_bg      = ', B_bg
  WRITE(*,'(A,ES12.4)') 'expect T  = ', lam/v_alf
END IF

!---Boundary conditions
! n, T, velx, velz frozen everywhere; psi default (Dirichlet walls, pins background);
! vely, by natural everywhere (x-uniform wave, periodic z)
ALLOCATE(mhd_sim%n_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%n_bc = .TRUE.
mhd_sim%t_bc => mhd_sim%n_bc
mhd_sim%velx_bc => mhd_sim%n_bc
mhd_sim%velz_bc => mhd_sim%n_bc
ALLOCATE(mhd_sim%vely_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%vely_bc = .FALSE.
mhd_sim%by_bc => mhd_sim%vely_bc

!---Initial conditions via projection
NULLIFY(u, v, mop)
CALL oft_blag_getmop(ML_oft_blagrange%current_level, mop)
CALL create_cg_solver(minv)
minv%A => mop
minv%its = -2
CALL create_diag_pre(minv%pre)
CALL ML_oft_blagrange%vec_create(u)
CALL ML_oft_blagrange%vec_create(v)
field_init%mesh => mesh

! n
CALL u%set(n0); CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,1)
DEALLOCATE(vec_vals); NULLIFY(vec_vals)
! velx, velz zero
CALL u%set(0.d0); CALL u%get_local(vec_vals)
CALL mhd_sim%u%restore_local(vec_vals,2); CALL mhd_sim%u%restore_local(vec_vals,4)
DEALLOCATE(vec_vals); NULLIFY(vec_vals)
! vely = v_delta * sin(2*pi*z/lam)
field_init%func => vely_wave
CALL oft_blag_project(ML_oft_blagrange%current_level, field_init, v)
CALL u%set(0.d0); CALL minv%apply(u, v)
CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,3)
DEALLOCATE(vec_vals); NULLIFY(vec_vals)
! T
CALL u%set(t0); CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,5)
DEALLOCATE(vec_vals); NULLIFY(vec_vals)
! psi = -B_bg * x (uniform in-plane background field via the flux function)
field_init%func => psi_bg
CALL oft_blag_project(ML_oft_blagrange%current_level, field_init, v)
CALL u%set(0.d0); CALL minv%apply(u, v)
CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,6)
DEALLOCATE(vec_vals); NULLIFY(vec_vals)
! by = 0
CALL u%set(0.d0); CALL u%get_local(vec_vals); CALL mhd_sim%u%restore_local(vec_vals,7)
DEALLOCATE(vec_vals); NULLIFY(vec_vals)

!---Run; kinetic-energy history lands in xmhd_2d.hist (diag column 2)
CALL mhd_sim%run_simulation()
IF(oft_env%head_proc) WRITE(*,'(A)') 'Run complete; check KE column of xmhd_2d.hist'
CALL oft_finalize()
END PROGRAM test_alfven_by2d
