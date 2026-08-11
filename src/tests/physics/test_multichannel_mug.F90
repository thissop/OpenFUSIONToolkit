!---------------------------------------------------------------------------
!> @file test_multichannel_mug.F90
!
!> Segmented multi-channel duct (LM-MHD flagship item 3): two out-of-plane
!> flow channels (mesh regions 1 and 3) separated by a conducting structural
!> septum (region 2), built via the cube_options septum_lo/septum_hi tags.
!> The septum is solid (velocity pinned to 0 = no-slip wall) but carries the
!> induced field with its own resistivity eta_reg(2)=eta_wall, so induced
!> current can LEAK across it and electromagnetically couple the two channels.
!>   - eta_wall large  -> insulating septum -> decoupled Shercliff/Hunt channels
!>   - eta_wall finite -> conducting septum -> inter-channel leakage cross-talk
!> Drive: steady pumping (fy) and/or a disruption dB/dt ramp (by_source). The
!> per-channel flow rate Q_i(t) and septum leakage current I_leak(t) are
!> reconstructed in post from the (x,z,vely,by) profile.
!> Single-region limit (no septum in the mesh) reduces exactly to test_hunt_mug.
!---------------------------------------------------------------------------
PROGRAM test_multichannel_mug
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
INTEGER(i4) :: ittarget = 1000000
REAL(r8) :: dt          = 2.d-2
REAL(r8) :: a_half      = 1.d0       !< Hartmann-wall half-gap in z (domain must match oft.in)
REAL(r8) :: B0          = 5.6d-2     !< applied field (in-plane z); Ha=B0*a/sqrt(mu0*eta*rho*nu)
REAL(r8) :: eta         = 1.d0       !< fluid magnetic diffusivity (resistivity in mu0 units)
REAL(r8) :: nu          = 1.d0       !< kinematic viscosity
REAL(r8) :: fy          = 1.d-2      !< body force per unit mass (axial pumping; 0 for pure disruption)
REAL(r8) :: n0          = 1.d0
REAL(r8) :: t0          = 1.d0
REAL(r8) :: chi         = 1.d-30
REAL(r8) :: D_diff      = 1.d-30
REAL(r8) :: gamma       = 1.67d0
REAL(r8) :: den_scale   = 1.d0
REAL(r8) :: lin_tol     = 1.d-9
REAL(r8) :: nl_tol      = 1.d-7
REAL(r8) :: c_wall      = 0.d0       !< Hartmann-wall (z=+/-a) conductance ratio (0 = perfect conductor)
REAL(r8) :: eta_wall    = 1.d0       !< septum (region 2) resistivity ratio eta_reg(2); >>1 insulating, <1 conducting
CHARACTER(LEN=256) :: restart_file = ''
INTEGER(i4) :: rst_base = 0
LOGICAL :: pm = .FALSE.
REAL(r8) :: by_source_dB   = 0.d0    !< disruption ramp total field change dB (0 => off)
REAL(r8) :: by_source_tauq = 0.d0    !< disruption ramp duration tau_q (0 => off)
LOGICAL :: use_ilu = .FALSE.
INTEGER(i4) :: pre_freq = 1
NAMELIST/multichannel_options/ order, nsteps, rst_freq, ittarget, dt, a_half, B0, eta, &
  nu, fy, n0, t0, chi, D_diff, gamma, den_scale, lin_tol, nl_tol, c_wall, eta_wall, &
  restart_file, rst_base, pm, by_source_dB, by_source_tauq, use_ilu, pre_freq

TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
CLASS(oft_vector), POINTER :: u => NULL()
REAL(r8) :: Ha_pred, rho, Lx, Lz
INTEGER(i4) :: nreg

CALL oft_init
OPEN(NEWUNIT=io_unit, FILE=oft_env%ifile)
READ(io_unit, multichannel_options, IOSTAT=ierr)
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
mhd_sim%use_ilu   = use_ilu
mhd_sim%pre_freq  = pre_freq
mhd_sim%B_0       = [0.d0, 0.d0, B0]    ! Hartmann direction (in-plane z)
mhd_sim%use_body_force = .TRUE.
mhd_sim%body_force = [0.d0, fy, 0.d0]   ! axial (out-of-plane) pumping (fy=0 for pure disruption)
mhd_sim%use_by_source  = (by_source_tauq > 0.d0)
mhd_sim%by_source_dB   = by_source_dB
mhd_sim%by_source_tauq = by_source_tauq

rho = proton_mass * n0 * den_scale
Ha_pred = B0 * a_half / SQRT(mu0 * eta * rho * nu)
IF(oft_env%head_proc) WRITE(*,'(A,ES12.4)') 'Predicted Ha = ', Ha_pred

!---Region materials: region 2 = structural septum (solid + its own resistivity).
!   Only engaged when the mesh actually carries a septum (nreg>=2, via septum_lo/hi in
!   &cube_options). Single-region mesh => eta_reg unassociated => bit-identical to Hunt.
nreg = mg_mesh%smesh%nreg
IF(nreg >= 2)THEN
  ALLOCATE(mhd_sim%eta_reg(nreg), mhd_sim%solid_reg(nreg))
  mhd_sim%eta_reg   = 1.d0
  mhd_sim%solid_reg = .FALSE.
  mhd_sim%eta_reg(2)   = eta_wall   ! septum resistivity ratio (leakage coupling knob)
  mhd_sim%solid_reg(2) = .TRUE.
  IF(oft_env%head_proc) WRITE(*,'(A,I0,A,ES12.4)') &
    'Multi-channel: nreg = ', nreg, ', septum eta_wall = ', eta_wall
END IF

!---Boundary conditions
! n, T frozen everywhere (incompressible, no pressure force):
ALLOCATE(mhd_sim%n_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%n_bc = .TRUE.
mhd_sim%t_bc => mhd_sim%n_bc
! psi: natural BC:
ALLOCATE(mhd_sim%psi_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%psi_bc = .FALSE.
mhd_sim%c_wall = c_wall
mhd_sim%a_half = a_half

!---Domain extents (for coordinate-based wall detection; order=1 => DOF index == node index)
Lx = MAXVAL(ABS(mg_mesh%smesh%r(1,:)))
Lz = MAXVAL(ABS(mg_mesh%smesh%r(2,:)))

! by BC: insulating OUTER side walls (x=+/-Lx) Dirichlet by=0; Hartmann walls (z=+/-Lz)
! free (natural / c_wall Robin); septum FREE so induced current leaks across it.
ALLOCATE(mhd_sim%by_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%by_bc = .FALSE.
BLOCK
INTEGER(i4) :: ip
DO ip=1,mg_mesh%smesh%np
  IF(ABS(ABS(mg_mesh%smesh%r(1,ip)) - Lx) < 1.d-8) mhd_sim%by_bc(ip) = .TRUE.
END DO
END BLOCK

! Velocity BC: no-slip on every OUTER wall (x=+/-Lx, z=+/-Lz) AND on the septum (region 2),
! which pins velocity=0 in the structural wall. Build one mask, share it across components.
ALLOCATE(mhd_sim%velx_bc(ML_oft_blagrange%current_level%ne))
mhd_sim%velx_bc = .FALSE.
BLOCK
INTEGER(i4) :: ip, ic, j, cnp
DO ip=1,mg_mesh%smesh%np
  IF(ABS(ABS(mg_mesh%smesh%r(1,ip)) - Lx) < 1.d-8) mhd_sim%velx_bc(ip) = .TRUE.
  IF(ABS(ABS(mg_mesh%smesh%r(2,ip)) - Lz) < 1.d-8) mhd_sim%velx_bc(ip) = .TRUE.
END DO
IF(nreg >= 2)THEN
  cnp = SIZE(mg_mesh%smesh%lc, 1)
  DO ic=1,mg_mesh%smesh%nc
    IF(mg_mesh%smesh%reg(ic) == 2)THEN
      DO j=1,cnp
        mhd_sim%velx_bc(mg_mesh%smesh%lc(j,ic)) = .TRUE.
      END DO
    END IF
  END DO
END IF
END BLOCK
mhd_sim%vely_bc => mhd_sim%velx_bc
mhd_sim%velz_bc => mhd_sim%velx_bc

!---Initial conditions
CALL ML_oft_blagrange%vec_create(u)
mhd_sim%rst_base = rst_base
IF(LEN_TRIM(restart_file) > 0)THEN
  CALL mhd_sim%rst_load(mhd_sim%u, TRIM(restart_file), 'U', mhd_sim%t, mhd_sim%dt)
  IF(oft_env%head_proc) WRITE(*,'(A,A,A,ES12.4,A,ES12.4)') &
    'Restarted from ', TRIM(restart_file), ' at t=', mhd_sim%t, ' dt=', mhd_sim%dt
ELSE
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

!---March
CALL mhd_sim%run_simulation()

!---Extract vely and by over all mesh points (order=1: DOF=vertex). Post computes
!   per-channel Q (integrate vely over x<0 vs x>0, excluding the septum band) and
!   the septum leakage current from d(by)/dx across the septum.
BLOCK
REAL(r8), POINTER :: by_vals(:)
NULLIFY(by_vals)
CALL mhd_sim%u%get_local(vec_vals, 3)
CALL mhd_sim%u%get_local(by_vals, 7)
np = mg_mesh%smesh%np
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit, FILE='multichannel_mug.profile')
  WRITE(io_unit,'(A,ES16.8)') '# Ha_pred = ', Ha_pred
  WRITE(io_unit,'(A,ES16.8)') '# fy = ', fy
  WRITE(io_unit,'(A,ES16.8)') '# eta_wall = ', eta_wall
  WRITE(io_unit,'(A,ES16.8)') '# c_wall = ', c_wall
  WRITE(io_unit,'(A)') '# x  z  vely  by   (u_nondim = vely*nu/(fy*a^2))'
  DO i=1,np
    WRITE(io_unit,'(4ES20.10)') mg_mesh%smesh%r(1,i), mg_mesh%smesh%r(2,i), vec_vals(i), by_vals(i)
  END DO
  CLOSE(io_unit)
  WRITE(*,'(A,I0,A)') 'Wrote multichannel_mug.profile (', np, ' points)'
END IF
DEALLOCATE(vec_vals, by_vals); NULLIFY(vec_vals, by_vals)
END BLOCK

CALL oft_finalize()
END PROGRAM test_multichannel_mug
