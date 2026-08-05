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
!> Per-step probe: extract fluid-core max|v_z| (|x|<a_half, i.e. excludes the
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
  REAL(r8) :: a_half = 0.5d0                 !< fluid half-width (probe samples |x|<a_half)
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
      IF(ABS(self%r(1,i)) < self%a_half)THEN
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

!------------------------------------------------------------------------------
!> Constant field-ramp driver (box-1's disruption Sdrv = dB/dt, constant): adds a
!> uniform induction source S = Sdrv*zhat to the B residual each step. The source
!> load b_S = INT(Sdrv*zhat . phi) dV is precomputed (an H(Curl)+Grad(H1) vector);
!> apply() adds scale*dt*b_S to the B blocks (1=H(Curl), 2=Grad(H1)) of residual v.
!> Sign/magnitude via `scale` (calibrate: with no flow, B_z must grow at ~Sdrv).
!------------------------------------------------------------------------------
MODULE ramp_driver_mod
USE oft_base
USE oft_la_base, ONLY: oft_vector
USE xmhd, ONLY: oft_xmhd_driver
IMPLICIT NONE
PRIVATE
PUBLIC :: ramp_driver
TYPE, EXTENDS(oft_xmhd_driver) :: ramp_driver
  CLASS(oft_vector), POINTER :: b_S => NULL()   !< precomputed source load (B-space)
  REAL(r8) :: scale = 1.d0                        !< calibration factor (sign+magnitude)
  REAL(r8) :: t_ramp = 0.d0                       !< smooth turn-on time (0 = instant on)
CONTAINS
  PROCEDURE :: apply => ramp_driver_apply
END TYPE ramp_driver
CONTAINS
SUBROUTINE ramp_driver_apply(self,up,u,v,t,dt)
  CLASS(ramp_driver), INTENT(inout) :: self
  CLASS(oft_vector), INTENT(inout) :: up,u,v
  REAL(r8), INTENT(in) :: t,dt
  REAL(r8), POINTER :: svals(:),vvals(:)
  REAL(r8) :: f
  INTEGER(i4) :: ib
  ! smooth cosine turn-on over t_ramp kills the impulsive-on shock (which otherwise
  ! contaminates the early overshoot at low Ha); f=1 for instant-on (t_ramp<=0)
  IF(self%t_ramp>0.d0 .AND. t<self%t_ramp)THEN
    f=0.5d0*(1.d0-COS(3.14159265358979d0*t/self%t_ramp))
  ELSE
    f=1.d0
  END IF
  DO ib=1,2                                       ! B occupies residual blocks 1,2
    NULLIFY(svals); NULLIFY(vvals)
    CALL self%b_S%get_local(svals,iblock=ib)
    CALL v%get_local(vvals,iblock=ib)
    vvals = vvals + f*self%scale*dt*svals
    CALL v%restore_local(vvals,iblock=ib,wait=.TRUE.)
    DEALLOCATE(svals,vvals)
  END DO
END SUBROUTINE ramp_driver_apply
END MODULE ramp_driver_mod

!------------------------------------------------------------------------------
!> Region-aware disruption source: S = Sdrv*zhat in the FLUID, zero in the
!> conjugate wall (region 2). Projecting THIS (instead of a uniform field) makes
!> the source load b_S fluid-only, so the wall is driven purely by conjugate
!> diffusion of the induced field (box-1's reference model). Reduces to the
!> whole-domain uniform source when no wall is tagged (all cells /= region 2).
!------------------------------------------------------------------------------
MODULE region_source_mod
USE oft_base
USE fem_utils, ONLY: fem_interp
IMPLICIT NONE
PRIVATE
PUBLIC :: region_source_field
TYPE, EXTENDS(fem_interp) :: region_source_field
  INTEGER(i4) :: n = 3                                  !< component count
  REAL(r8) :: sval = 0.d0                               !< Sdrv (axial source amplitude)
  INTEGER(i4), POINTER, DIMENSION(:) :: reg => NULL()   !< => mesh%reg (per-cell region id)
CONTAINS
  PROCEDURE :: interp => region_source_interp
END TYPE region_source_field
CONTAINS
SUBROUTINE region_source_interp(self,cell,f,gop,val)
  CLASS(region_source_field), INTENT(inout) :: self
  INTEGER(i4), INTENT(in) :: cell
  REAL(r8), INTENT(in) :: f(:)
  REAL(r8), INTENT(in) :: gop(3,4)
  REAL(r8), INTENT(out) :: val(:)
  val=0.d0
  IF(self%reg(cell) /= 2) val(3)=self%sval   ! fluid (any non-wall region) gets the source
END SUBROUTINE region_source_interp
END MODULE region_source_mod

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
USE ramp_driver_mod, ONLY: ramp_driver
USE region_source_mod, ONLY: region_source_field
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
TYPE(ramp_driver) :: driver
TYPE(uniform_field) :: sdrv_field
TYPE(region_source_field) :: sdrv_rfield
CLASS(oft_vector), POINTER :: b_S => NULL()
INTEGER(i4) :: io_unit,ierr,i,j,np,probe_unit
REAL(r8) :: cx
REAL(r8), POINTER :: vel_vals(:) => NULL()
!---Runtime options (namelist)
INTEGER(i4) :: order  = 2
INTEGER(i4) :: minlev = 1
REAL(r8) :: B0    = 1.d-2     !< transverse Hartmann field (x); Ha = B0*a*sqrt(sigma/(rho*nu))
REAL(r8) :: drive = 1.d0      !< body-force acceleration g_accel (-z axial drive)
REAL(r8) :: n0    = 1.d19     !< number density
REAL(r8) :: t0    = 1.d0      !< temperature [eV]
REAL(r8) :: a_half = 0.5d0    !< duct half-width (for the core velocity mask)
CHARACTER(LEN=8) :: drive_mode = 'force'  !< 'force' = g_accel body force; 'ramp' = Sdrv field-ramp source
REAL(r8) :: Sdrv  = 0.d0      !< constant induction source dB/dt (z-component), box-1's disruption driver
REAL(r8) :: sdrv_scale = 1.d0 !< calibration factor for the ramp source (sign+magnitude)
REAL(r8) :: sdrv_tramp = 0.d0 !< smooth source turn-on time (0 = instant); kills the turn-on shock
REAL(r8) :: wall_tw = 0.d0    !< conjugate Hartmann-wall thickness (a_half<|x|<a_half+wall_tw -> region 2); 0 = no wall
NAMELIST/hartmann3d_options/order,minlev,B0,drive,n0,t0,a_half,drive_mode,Sdrv,sdrv_scale,sdrv_tramp,wall_tw
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
! Conjugate wall (finite-c): tag Hartmann-wall cells a_half<|x|<a_half+wall_tw as
! region 2. The solver's <xmhd><region id=2 type=2 (solid)> block then makes them a
! solid wall (no flow) with finite resistivity eta_reg -> c=(1/eta_reg)*(tw/a).
! The box must extend to +/-(a_half+wall_tw) in x (set via &cube_options rscale/shift).
!------------------------------------------------------------------------------
IF(wall_tw > 0.d0)THEN
  DO i=1,mg_mesh%mesh%nc
    cx=0.d0
    DO j=1,SIZE(mg_mesh%mesh%lc,1)
      cx=cx+mg_mesh%mesh%r(1,mg_mesh%mesh%lc(j,i))
    END DO
    cx=cx/REAL(SIZE(mg_mesh%mesh%lc,1),8)
    IF(ABS(cx) > a_half .AND. ABS(cx) < a_half+wall_tw+1.d-9) mg_mesh%mesh%reg(i)=2
  END DO
  mg_mesh%mesh%nreg=MAX(mg_mesh%mesh%nreg,2)
  IF(oft_env%head_proc)WRITE(*,'(A,I0,A,I0,A)') 'Conjugate wall: ', &
    COUNT(mg_mesh%mesh%reg==2),' of ',mg_mesh%mesh%nc,' cells -> region 2 (solid wall)'
END IF
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
! Drive: 'force' = uniform body force g_accel(-z); 'ramp' = Sdrv field-ramp source
! (box-1's disruption). Then run with the driver (ramp) + probe.
!------------------------------------------------------------------------------
xmhd_minlev=minlev
xmhd_taxis=2
oft_env%pm=.FALSE.
equil_fields%B=>b
equil_fields%V=>vel
equil_fields%Ne=>den
equil_fields%Ti=>temp
IF(TRIM(drive_mode)=='ramp')THEN
  g_accel=0.d0
  ! Fluid-only disruption source b_S = INT(Sdrv*zhat . phi) over region 1 (fluid);
  ! the conjugate wall (region 2) carries no source, driven only by field diffusion
  ! (box-1's reference). Reduces to a whole-domain source when no wall is tagged.
  sdrv_rfield%sval=Sdrv
  sdrv_rfield%reg=>mg_mesh%mesh%reg
  CALL xmhd_ML_hcurl_grad%vec_create(b_S)
  CALL oft_hcurl_grad_project(xmhd_ML_hcurl_grad%current_level,sdrv_rfield,b_S)
  CALL hcurl_grad_gzerop%apply(b_S)
  driver%b_S=>b_S
  driver%scale=sdrv_scale
  driver%t_ramp=sdrv_tramp
  CALL xmhd_run(equil_fields,driver=driver,probes=probe)
ELSE
  g_accel=drive
  CALL xmhd_run(equil_fields,probes=probe)
END IF
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
