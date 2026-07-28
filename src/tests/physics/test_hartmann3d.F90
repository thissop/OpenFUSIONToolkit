!------------------------------------------------------------------------------
!> @file test_hartmann3d.F90
!
!> 3D Hartmann duct feasibility spike for the xMHD (3D full-induction) solver.
!> Axial (z) flow in a duct, driven by the built-in uniform body force g_accel
!> (-z), braked by a uniform transverse Hartmann field B0 = B0 xhat. Duct is
!> walled in x (Hartmann walls, x=+/-a) and y (side walls), periodic in z
!> (fully-developed). xmhd_run applies no-slip v + Dirichlet-B on the boundary
!> automatically (gbe). Low Ha => thick Hartmann layer => cheap to resolve.
!> Goal: does the 3D solver reproduce the analytic Hartmann velocity profile,
!> and (transient) does the core velocity overshoot on the standing-Alfven clock?
!> A per-step probe (hartmann_probe) writes the core/centerline axial velocity
!> vs time to hartmann3d.probe -> the overshoot amplitude A_os for paper-1 sec.5.
!> Stage 1 of the MUG 3D LM-MHD plan (docs/dev/MUG_3D_LM_MHD_PLAN.md).
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Per-step probe: extract core-masked max|v_z| (|x|<0.5*a_half) and the
!> centerline v_z (node nearest the x=y=0 axis) from the live velocity field.
!------------------------------------------------------------------------------
MODULE hartmann_probe_mod
USE oft_base
USE oft_la_base, ONLY: oft_vector
USE xmhd, ONLY: oft_xmhd_probe, xmhd_sub_fields
IMPLICIT NONE
PRIVATE
PUBLIC :: hartmann_probe
TYPE, EXTENDS(oft_xmhd_probe) :: hartmann_probe
  REAL(r8), POINTER :: r(:,:) => NULL()     !< mesh vertex coords (=>mg_mesh%mesh%r)
  REAL(r8), POINTER :: vvals(:) => NULL()   !< reused get_local buffer
  INTEGER(i4) :: np = 0                      !< number of mesh vertices
  INTEGER(i4) :: io_unit = 0                 !< open unit (NEWUNIT is negative; 0 = inactive)
  REAL(r8) :: a_half = 0.5d0                 !< duct half-width (core mask |x|<0.5*a_half)
CONTAINS
  PROCEDURE :: apply => hartmann_probe_apply
END TYPE hartmann_probe
CONTAINS
SUBROUTINE hartmann_probe_apply(self,sub_fields,t)
  CLASS(hartmann_probe), INTENT(inout) :: self
  TYPE(xmhd_sub_fields), INTENT(inout) :: sub_fields
  REAL(r8), INTENT(in) :: t
  REAL(r8) :: cmax(3),vcen(3),d,dmin
  INTEGER(i4) :: ib,i,nn
  ! Vector-Lagrange is BLOCKED: get_local(vals,ib) returns component-ib block
  ! (block 1=x,2=y,3=z by convention). Compute per-component core-max |v| over
  ! |x|<0.5*a_half and the value at the node nearest the x=y=0 axis. Writing all
  ! three lets us confirm which block is the axial (drive) velocity.
  dmin=HUGE(1.d0)
  cmax=0.d0; vcen=0.d0
  DO ib=1,3
    CALL sub_fields%V%get_local(self%vvals,ib)   ! reuses buffer after 1st alloc
    nn=MIN(self%np,SIZE(self%vvals))
    DO i=1,nn
      IF(ABS(self%r(1,i)) < 0.5d0*self%a_half)THEN
        IF(ABS(self%vvals(i))>cmax(ib)) cmax(ib)=ABS(self%vvals(i))
      END IF
    END DO
  END DO
  ! centerline node (nearest x=y=0): capture all 3 components there
  DO ib=1,3
    CALL sub_fields%V%get_local(self%vvals,ib)
    dmin=HUGE(1.d0)
    DO i=1,MIN(self%np,SIZE(self%vvals))
      d=self%r(1,i)**2+self%r(2,i)**2
      IF(d<dmin)THEN; dmin=d; vcen(ib)=self%vvals(i); END IF
    END DO
  END DO
  IF(self%io_unit/=0)THEN
    WRITE(self%io_unit,'(7ES18.8)') t,cmax(1),cmax(2),cmax(3),vcen(1),vcen(2),vcen(3)
    FLUSH(self%io_unit)
  END IF
END SUBROUTINE hartmann_probe_apply
END MODULE hartmann_probe_mod

PROGRAM test_hartmann3d
USE oft_base
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_vgetmop
!---H1 (Grad) subspace
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp
!---H(Curl) space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_setup_interp, hcurl_grad_getmop, &
  oft_hcurl_grad_project, oft_hcurl_grad_gzerop
!---Physics
USE oft_vector_inits, ONLY: uniform_field
USE mhd_utils, ONLY: mu0, proton_mass
USE xmhd, ONLY: xmhd_run, xmhd_plot, xmhd_minlev, xmhd_taxis, g_accel, &
  xmhd_sub_fields, xmhd_ML_hcurl, xmhd_ML_H1, xmhd_ML_hcurl_grad, xmhd_ML_H1grad, &
  xmhd_ML_lagrange, xmhd_ML_vlagrange
USE hartmann_probe_mod, ONLY: hartmann_probe
IMPLICIT NONE
!---Full H(Curl) space mass-matrix solver
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
!---Fields
CLASS(oft_vector), POINTER :: b,vel,u,v,be,den,temp
TYPE(xmhd_sub_fields) :: equil_fields
TYPE(uniform_field) :: x_field
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_hcurl_grad_gzerop), TARGET :: hcurl_grad_gzerop
TYPE(hartmann_probe) :: probe
INTEGER(i4) :: io_unit,ierr,i,np,probe_unit
REAL(r8), POINTER :: vel_vals(:) => NULL()
!---Runtime options (namelist)
INTEGER(i4) :: order  = 2
INTEGER(i4) :: minlev = 1
REAL(r8) :: B0    = 1.d-2     !< transverse Hartmann field (x); Ha = B0*a*sqrt(sigma/(rho*nu))
REAL(r8) :: drive = 1.d0      !< body-force acceleration g_accel (-z axial drive)
REAL(r8) :: n0    = 1.d19     !< number density
REAL(r8) :: t0    = 1.d0      !< temperature [eV]
REAL(r8) :: a_half = 0.5d0    !< duct half-width (for the core velocity mask)
NAMELIST/hartmann3d_options/order,minlev,B0,drive,n0,t0,a_half
!------------------------------------------------------------------------------
! Initialize + read options
!------------------------------------------------------------------------------
CALL oft_init
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,hartmann3d_options,IOSTAT=ierr)
CLOSE(io_unit)
!------------------------------------------------------------------------------
! Grid (duct box from &cube_options: walled in x,y; periodic in z)
!------------------------------------------------------------------------------
CALL multigrid_construct(mg_mesh)
!------------------------------------------------------------------------------
! FE structures (mirrors test_alfven)
!------------------------------------------------------------------------------
ALLOCATE(xmhd_ML_lagrange,xmhd_ML_vlagrange)
CALL oft_lag_setup(mg_mesh,order,xmhd_ML_lagrange,ML_vlag_obj=xmhd_ML_vlagrange,minlev=minlev)
CALL lag_setup_interp(xmhd_ML_lagrange)
ALLOCATE(xmhd_ML_H1)
CALL oft_h1_setup(mg_mesh,order+1,xmhd_ML_H1,minlev=minlev)
CALL h1_setup_interp(xmhd_ML_H1)
ALLOCATE(xmhd_ML_hcurl)
CALL oft_hcurl_setup(mg_mesh,order,xmhd_ML_hcurl,minlev=minlev)
CALL hcurl_setup_interp(xmhd_ML_hcurl)
ALLOCATE(xmhd_ML_hcurl_grad,xmhd_ML_H1grad)
CALL oft_hcurl_grad_setup(xmhd_ML_hcurl,xmhd_ML_H1,xmhd_ML_hcurl_grad,xmhd_ML_H1grad,minlev)
CALL hcurl_grad_setup_interp(xmhd_ML_hcurl_grad,xmhd_ML_H1)
hcurl_grad_gzerop%ML_hcurl_grad_rep=>xmhd_ML_hcurl_grad
!------------------------------------------------------------------------------
! Full H(Curl) mass-matrix solver (to project the uniform B0)
!------------------------------------------------------------------------------
NULLIFY(mop)
CALL hcurl_grad_getmop(xmhd_ML_hcurl_grad%current_level,mop,"none")
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
minv%atol=1.d-10
CALL create_diag_pre(minv%pre)
CALL xmhd_ML_hcurl_grad%vec_create(u)
CALL xmhd_ML_hcurl_grad%vec_create(v)
CALL xmhd_ML_hcurl_grad%vec_create(b)
CALL xmhd_ML_hcurl_grad%vec_create(be)
!------------------------------------------------------------------------------
! Uniform transverse Hartmann field B0 = B0 * xhat
!------------------------------------------------------------------------------
x_field%val=(/B0,0.d0,0.d0/)
CALL oft_hcurl_grad_project(xmhd_ML_hcurl_grad%current_level,x_field,v)
CALL hcurl_grad_gzerop%apply(v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL b%add(0.d0,1.d0,u)
CALL be%add(0.d0,1.d0,u)
!------------------------------------------------------------------------------
! Velocity (rest), density, temperature
!------------------------------------------------------------------------------
CALL xmhd_ML_vlagrange%vec_create(vel)
CALL vel%set(0.d0)
CALL xmhd_ML_lagrange%vec_create(den)
CALL xmhd_ML_lagrange%vec_create(temp)
CALL den%set(n0)
CALL temp%set(t0)
!------------------------------------------------------------------------------
! Per-step probe -> hartmann3d.probe (t, core-max|vz|, centerline vz)
!------------------------------------------------------------------------------
probe%r=>mg_mesh%mesh%r
probe%np=mg_mesh%mesh%np
probe%a_half=a_half
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=probe_unit,FILE='hartmann3d.probe')
  WRITE(probe_unit,'(A)') '# t  coremax_vx coremax_vy coremax_vz  cen_vx cen_vy cen_vz'
  probe%io_unit=probe_unit
END IF
!------------------------------------------------------------------------------
! Axial drive via built-in body force (-z), then run (with probe)
!------------------------------------------------------------------------------
g_accel=drive
xmhd_minlev=minlev
xmhd_taxis=2
oft_env%pm=.FALSE.
equil_fields%B=>b
equil_fields%V=>vel
equil_fields%Ne=>den
equil_fields%Ti=>temp
CALL xmhd_run(equil_fields,probes=probe)
IF(probe%io_unit/=0) CLOSE(probe%io_unit)
!------------------------------------------------------------------------------
! Dump solution (xmhd_plot) + raw axial-velocity profile (order-1 vertex sample)
!------------------------------------------------------------------------------
CALL xmhd_plot
CALL vel%get_local(vel_vals)
np=mg_mesh%mesh%np
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='hartmann3d.profile')
  WRITE(io_unit,'(A,ES12.4,A,ES12.4)') '# B0 = ',B0,'  g_accel = ',drive
  WRITE(io_unit,'(A)') '# x  y  z  vel(3*(i-1)+1..3)   (order>=2: vertex DOFs only)'
  DO i=1,MIN(np,SIZE(vel_vals)/3)
    WRITE(io_unit,'(6ES18.8)') mg_mesh%mesh%r(1,i),mg_mesh%mesh%r(2,i),mg_mesh%mesh%r(3,i), &
      vel_vals(3*(i-1)+1),vel_vals(3*(i-1)+2),vel_vals(3*(i-1)+3)
  END DO
  CLOSE(io_unit)
  WRITE(*,'(A,I0,A)') 'Wrote hartmann3d.profile (',np,' points)'
END IF
CALL oft_finalize
END PROGRAM test_hartmann3d
