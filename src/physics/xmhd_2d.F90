!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file xmhd_2d.F90
!
!> Solve time-dependent MHD equations in 2D with lagrange basis set
!---------------------------------------------------------------------------
MODULE xmhd_2d
#if !defined(XMHD_RST_LEN)
#define XMHD_RST_LEN 5
#endif
#define XMHD_ITCACHE 10
USE oft_base
USE oft_io, ONLY: hdf5_read, hdf5_write, oft_file_exist, &
  hdf5_field_exist, oft_bin_file, xdmf_plot_file
USE oft_quadrature
USE oft_mesh_type, ONLY: oft_bmesh, cell_is_curved
USE multigrid, ONLY: multigrid_mesh
!
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_local_mat, oft_vector_ptr, &
  vector_extrapolate, oft_graph, oft_graph_ptr
USE oft_solver_utils, ONLY: create_solver_xml, create_diag_pre, create_ilu_pre
USE oft_deriv_matrices, ONLY: oft_noop_matrix, oft_mf_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_nksolver, oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_fem_comp_type
USE fem_utils, ONLY: fem_dirichlet_diag, fem_dirichlet_vec, bfem_map_flag
USE oft_lag_basis, ONLY: oft_lag_setup,oft_scalar_bfem, oft_blag_eval, oft_blag_geval, oft_2D_lagrange_cast
USE oft_blag_operators, ONLY: oft_blag_vproject,oft_blag_project, oft_blag_getmop, oft_lag_bginterp
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: mu0, elec_charge, proton_mass
USE oft_gs, ONLY: gs_epsilon
USE oft_la_utils, ONLY: create_matrix
IMPLICIT NONE
#include "local.h"
#if !defined(TDIFF_RST_LEN)
#define TDIFF_RST_LEN 5
#endif
PRIVATE

!------------------------------------------------------------------------------
!> Object to compute left-hand side of system for nonlinear solves
!------------------------------------------------------------------------------
TYPE, extends(oft_noop_matrix) :: xmhd_2d_nlfun
  REAL(r8) :: dt = -1.d0 !< time step
  REAL(r8) :: diag_vals(5) = 0.d0 !< diagnostic values
  REAL(r8) :: obs_vals(5) = 0.d0 !< LM-MHD observables (Cartesian disruption): (1) Fl_net = int (B0/mu0)*d(by)/dz
                                 !< [Lorentz force density], (2) Fl_abs = int |(B0/mu0)*d(by)/dz|,
                                 !< (3) EJf_rate = int eta*|grad by|^2/mu0 [= eta_res*J^2, since MUG eta = eta_res/mu0],
                                 !< (4) Iw = oint_walls |by|/mu0 dl [thin-wall sheet current K=by_wall/mu0],
                                 !< (5) EJw_rate = oint_walls eta*by^2/(c_wall*a*mu0) dl [thin-wall wall Joule; a=a_half]
  TYPE(oft_xmhd_2d_sim), pointer :: parent_sim => NULL() !< pointer to parent simulation object for access to parameters
  CONTAINS
  !> Apply the matrix
  PROCEDURE :: apply_real => nlfun_apply
END TYPE xmhd_2d_nlfun

!------------------------------------------------------------------------------
!> Object to compute left-hand side of system for linear solves
!------------------------------------------------------------------------------
TYPE, extends(oft_noop_matrix) :: xmhd_2d_mfun
  REAL (r8) :: diag_vals(5) = 0.d0 !< diagnostic values
  TYPE(oft_xmhd_2d_sim), pointer :: parent_sim => NULL() !< pointer to parent simulation object for access to parameters
  CONTAINS
  !> Apply the matrix
  PROCEDURE :: apply_real => mfun_apply
END TYPE xmhd_2d_mfun
!------------------------------------------------------------------------------
!> Main 2D MHD simulation object
!------------------------------------------------------------------------------
TYPE, public :: oft_xmhd_2d_sim
  LOGICAL :: mfnk = .FALSE. !< Use matrix free method?
  LOGICAL :: cyl_flag = .FALSE. !Use cylindrical coordinates
  LOGICAL :: linear = .FALSE. !Run sim. in linear mode
  LOGICAL :: timestep_cn = .FALSE. ! Use Crank-Nicolson timestep?
  INTEGER(i4) :: nsteps = -1 !< # of timesteps
  INTEGER(i4) :: rst_base = 0 !< starting index of restart files
  INTEGER(i4) :: rst_freq = 10 !< frequency to generate restart files
  INTEGER(i4) :: pre_freq = 1 !< frequency to update preconditioner
  INTEGER(i4) :: ittarget = 40 !< Needs docs
  REAL(r8) :: dt = -1.d0 !< timestep
  REAL(r8) :: jac_dt = -1.d0
  REAL(r8) :: t = 0.d0 !< current time
  REAL(r8) :: chi = -1.d0 !< thermal diffusivity
  REAL(r8) :: nu = -1.d0 !< kinematic viscosity [m^2/2]
  REAL(r8) :: gamma = 2.d0/3.d0 !< adiabatic index
  REAL(r8) :: D_diff = -1.d0 !< diffusivity
  REAL(r8) :: k_boltz = elec_charge !< Boltzmann constant (use elec_charge for T in eV)
  REAL(r8) :: den_scale = 1.d19 !< scale factor used to normalize density
  REAL(r8) :: m_i=proton_mass !< proton mass
  REAL(r8) :: lin_tol = 1.d-8 !< absolute tolerance for linear solver
  REAL(r8) :: nl_tol = 1.d-5 !< absolute tolerance for nonlinear solver
  REAL (r8) :: B_0(3) = 0.d0 !< background, static magnetic field
  REAL(r8) :: eta !< electrical resistivity, in units of mu0
  !--- Reduced Hartmann / wall-drag parameters (LM-MHD extension, default off)
  LOGICAL :: use_wall_drag = .FALSE. !< Enable reduced wall/Hartmann drag source term
  REAL(r8) :: drag_coeff = 0.d0 !< Drag coefficient alpha [1/time in problem units]
  REAL(r8) :: drag_bhat(3) = [0.d0,1.d0,0.d0] !< Unit vector along applied-field direction; drag is NOT applied along this direction
  !--- Effective-drag CLOSURE selector. With drag_closure /= 'raw', drag_coeff is
  !    computed from the correct geometry/conductance-aware Hartmann-layer closure
  !    alpha_eff = (nu/a^2) * f(Ha,c) (a = a_half). The naive 1-D Hartmann-CHANNEL
  !    closure f = Ha^2 OVER-BRAKES a rectangular duct by a factor ~Ha (it misses
  !    the Shercliff side-layer current return); use 'duct_*' for ducts/blankets.
  !    'raw'            : use drag_coeff as supplied (backward compatible)
  !    'channel'        : f = Ha^2                 (1-D Hartmann channel; core ~ Ha^-2)
  !    'duct_insulating': f = Ha                   (rectangular duct, insulating walls; core ~ Ha^-1)
  !    'duct_conducting': f = Ha^2 / (1 + 1/c)     (conducting Hunt duct; core ~ (1+1/c)/Ha^2)
  CHARACTER(LEN=16) :: drag_closure = 'raw' !< Hartmann-drag closure: raw|channel|duct_insulating|duct_conducting
  REAL(r8) :: drag_Ha = 0.d0 !< Hartmann number for the drag closure (drag_closure /= 'raw')
  REAL(r8) :: drag_c  = -1.d0 !< wall conductance ratio c for 'duct_conducting' closure
  !--- Uniform body-force (per unit mass) momentum source, default off.
  !    Used to drive fully-developed channel/duct flow (e.g. Hartmann verification).
  !    Adds S_u = body_force to the momentum equation (acceleration units, like the
  !    pressure-gradient and viscous terms). Constant => no Jacobian contribution.
  LOGICAL :: use_body_force = .FALSE. !< Enable uniform body-force momentum source
  REAL(r8) :: body_force(3) = [0.d0,0.d0,0.d0] !< Body force per unit mass [accel units]
  !--- LM-MHD disruption/ELM drive: a UNIFORM time-varying source S(t) = -dB0z/dt in the
  !    by (out-of-plane induced field) equation. Represents a decaying poloidal field ramping
  !    over tau_q. S(t) = (dB/tau_q)*(t<tau_q)*min(t/(tau_q/50),1); default off => no source
  !    (regression-safe). Uniform in space + independent of the solution => residual-only, no
  !    Jacobian contribution. Matches the P1 reduced-model source (docs comsol side).
  LOGICAL :: use_by_source = .FALSE. !< enable the disruption source in the by induction
  REAL(r8) :: by_source_dB   = 0.d0  !< total field change dB of the ramp [field units]
  REAL(r8) :: by_source_tauq = 0.d0  !< ramp duration tau_q [time units]
  !--- Conducting-wall (Robin) EM boundary condition: LM-MHD Phase-4, default off.
  !    Thin-wall condition on the axial induced field:  by + c_wall*d(by)/dn = 0, with
  !    c_wall = sigma_w*t_w/(sigma_f*a) the wall-conductance ratio (wall conductivity*thickness
  !    over fluid conductivity*half-width). c_wall = 0 (default) leaves behaviour unchanged
  !    (walls are then whatever the by_bc mask says: Dirichlet by=0 = insulating, or natural
  !    d(by)/dn=0 = perfectly conducting). When c_wall > 0 the Robin boundary term is applied
  !    on boundary edges whose endpoints are FREE in by_bc (driver masks select which walls);
  !    implemented in nlfun_apply (residual) and build_approx_jacobian ((7,7) block) via the
  !    exact order-1 edge mass matrix. Cartesian + order=1 only (guarded in setup). The
  !    analogous psi Robin term (in-plane conjugate cases) is future work.
  REAL(r8) :: c_wall = 0.d0 !< Wall-conductance ratio for thin-wall Robin by BC (0 = off)
  REAL(r8) :: a_half = 1.d0 !< duct half-width (for thin-wall wall-Joule observable EJw; set by driver)
  LOGICAL :: use_ilu = .FALSE. !< use native ILU(0) preconditioner instead of diagonal (for stiff high-Ha solves)
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: n_bc => NULL() !< n BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: velx_bc => NULL() !< vel BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: vely_bc => NULL() !< vel BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: velz_bc => NULL() !< vel BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: T_bc => NULL() !< T BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: psi_bc => NULL() !< psi BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: by_bc => NULL() !< by BC flag
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:,:) :: jacobian_block_mask => NULL() !< Matrix block mask
  TYPE(oft_fem_comp_type), POINTER :: fe_rep => NULL() !< Finite element representation for solution field
  TYPE(xdmf_plot_file) :: xdmf_plot
  CLASS(oft_vector), POINTER :: u => NULL() !< current solution vector
  CLASS(oft_vector), POINTER :: u0 => NULL() !< equilibrium solution vector (used for linearization if linear = True)
  CLASS(oft_matrix), POINTER :: jacobian => NULL() !< approximate jacobian matrix
  TYPE(oft_mf_matrix), POINTER :: mf_mat => NULL() !< Matrix free operator
  TYPE(xmhd_2d_nlfun), POINTER :: nlfun => NULL() !< Nonlinear function
  TYPE(xmhd_2d_mfun), POINTER :: mfun => NULL() !< mass matrix function
  TYPE(xml_node), POINTER :: xml_root => NULL() !< XML root element
  TYPE(xml_node), POINTER :: xml_pre_def => NULL() !< XML element for preconditioner definition
  contains
  !> Setup
  PROCEDURE :: setup => setup
  !> Set default boundary conditions
  PROCEDURE :: setup_bc => setup_bc
  !> Run simulation
  PROCEDURE :: run_simulation => run_simulation
  !> Run simulation in linear mode
  PROCEDURE :: run_lin_simulation => run_lin_simulation
  !> Save restart file
  PROCEDURE :: rst_save => rst_save
  !> Load restart file
  PROCEDURE :: rst_load => rst_load
END TYPE oft_xmhd_2d_sim
TYPE(oft_xmhd_2d_sim), POINTER :: current_sim => NULL() !Active 2D MHD simulation
!
CLASS(multigrid_mesh), POINTER :: mg_mesh => NULL() !< Multigrid mesh object
CLASS(oft_bmesh), POINTER, PUBLIC :: mesh => NULL() !< Current mesh level
TYPE(oft_ml_fem_type), TARGET, PUBLIC :: ML_oft_blagrange !< Multilevel finite element representation
CLASS(oft_scalar_bfem), POINTER :: oft_blagrange => NULL() !< Lagrange finite element representation
PUBLIC xmhd_2d_plot
CONTAINS

!---------------------------------------------------------------------------
!Main driver for 2D MHD nonlinear time advance
!---------------------------------------------------------------------------
subroutine run_simulation(self)
class(oft_xmhd_2d_sim), target, intent(inout) :: self
type(oft_nksolver) :: nksolver
!---Solver objects
class(oft_vector), pointer :: u,v,up
type(oft_native_gmres_solver), target :: solver
type(oft_timer) :: mytimer
!---History file fields
TYPE(oft_bin_file) :: hist_file
integer(i4) :: hist_i4(3)
real(4) :: hist_r4(6)
!---Timestep control
integer(i4) :: itcount(XMHD_ITCACHE)
real(r8) :: dtin,dthist(XMHD_ITCACHE)
!---
!---Extrapolation fields
integer(i4), parameter :: maxextrap=2
integer(i4) :: nextrap
real(r8), allocatable, dimension(:) :: extrapt
type(oft_vector_ptr), allocatable, dimension(:) :: extrap_fields
!---
character(LEN=TDIFF_RST_LEN) :: rst_char
integer(i4) :: i,j,io_stat,rst_tmp,npre
real(r8) :: n_avg, u_avg(3), T_avg, psi_avg, by_avg,elapsed_time
real(r8), pointer :: plot_vals(:),plot_vec(:,:)
!---LM-MHD disruption moments (per-step observables) — written only when use_by_source is on
integer(i4) :: mom_unit, ip_umax
real(r8) :: Umax
real(r8), pointer :: vely_arr(:)
current_sim=>self

!---------------------------------------------------------------------------
! Set boundary conditions not alreadys set
!---------------------------------------------------------------------------
CALL self%setup_bc()
!---------------------------------------------------------------------------
! Normalize wall-drag field direction (drag projects out the component of
! velocity perpendicular to drag_bhat, so drag_bhat must be a unit vector).
!---------------------------------------------------------------------------
IF(self%use_wall_drag)THEN
  IF(NORM2(self%drag_bhat)<=1.d-12)THEN
    CALL oft_abort('use_wall_drag=.TRUE. but drag_bhat is a zero vector', &
      'run_simulation',__FILE__)
  END IF
  self%drag_bhat=self%drag_bhat/NORM2(self%drag_bhat)
  !--- Resolve drag_coeff from the geometry/conductance-aware Hartmann-layer closure
  !    (see drag_closure above). alpha_eff = (nu/a^2)*f(Ha,c); 'raw' leaves the
  !    user-supplied drag_coeff untouched. The 'channel' closure over-brakes a duct
  !    by ~Ha and is provided only for reference; use 'duct_*' for ducts/blankets.
  IF(TRIM(self%drag_closure)/='raw' .AND. (self%nu<=0.d0 .OR. self%a_half<=0.d0)) &
    CALL oft_abort('drag_closure /= raw needs nu>0 and a_half>0','run_simulation',__FILE__)
  SELECT CASE (TRIM(self%drag_closure))
  CASE ('raw')
    ! keep the user-supplied drag_coeff
  CASE ('channel')
    self%drag_coeff = (self%nu/self%a_half**2)*self%drag_Ha**2
  CASE ('duct_insulating')
    self%drag_coeff = (self%nu/self%a_half**2)*self%drag_Ha
  CASE ('duct_conducting')
    IF(self%drag_c<=0.d0) CALL oft_abort( &
      "drag_closure='duct_conducting' requires drag_c>0",'run_simulation',__FILE__)
    self%drag_coeff = (self%nu/self%a_half**2)*self%drag_Ha**2/(1.d0+1.d0/self%drag_c)
  CASE DEFAULT
    CALL oft_abort("unknown drag_closure (raw|channel|duct_insulating|duct_conducting)", &
      'run_simulation',__FILE__)
  END SELECT
END IF
!---------------------------------------------------------------------------
! Create solver fields
!---------------------------------------------------------------------------
call self%fe_rep%vec_create(u)
call self%fe_rep%vec_create(up)
call self%fe_rep%vec_create(v)
!---

ALLOCATE(extrap_fields(maxextrap),extrapt(maxextrap))
DO i=1,maxextrap
  CALL self%fe_rep%vec_create(extrap_fields(i)%f)
  extrapt(i)=0.d0
END DO
nextrap=0
!---Set initial time and load initial solution vector from driver
self%t=0.d0
CALL u%add(0.d0,1.d0,self%u)
!---Create initial conditions restart file
104 FORMAT (I TDIFF_RST_LEN.TDIFF_RST_LEN)
WRITE(rst_char,104)0
CALL self%rst_save(u, self%t, self%dt, 'xmhd2d_'//rst_char//'.rst', 'U')

!---Fill timestep control cache
itcount=self%ittarget
dthist=self%dt
dtin=self%dt

!---------------------------------------------------------------------------
! Set parameters of nonlinear function
!---------------------------------------------------------------------------
ALLOCATE(self%nlfun)
self%nlfun%parent_sim=>self
!---------------------------------------------------------------------------
! Setup linear solver
!---------------------------------------------------------------------------
self%jac_dt=self%dt
IF(self%timestep_cn)self%jac_dt=self%dt/2.d0
CALL build_approx_jacobian(self,u)

IF(self%mfnk)THEN
  ALLOCATE(self%mf_mat)
  CALL up%set(1.d0)
  CALL self%mf_mat%setup(u,self%nlfun,up)
  self%mf_mat%b0=1.d-5
  solver%A=>self%mf_mat
ELSE
  solver%A=>self%jacobian
END IF
solver%its=self%ittarget*4
solver%atol=self%lin_tol
solver%itplot=1
solver%nrits=40
solver%pm=oft_env%pm
NULLIFY(solver%pre)
IF(ASSOCIATED(self%xml_pre_def))THEN
  CALL create_solver_xml(solver%pre,self%xml_pre_def)
ELSE IF(self%use_ilu)THEN
  CALL create_ilu_pre(solver%pre)   ! stronger than diagonal for stiff high-Ha solves
ELSE
  CALL create_diag_pre(solver%pre)
END IF
solver%pre%A=>self%jacobian

!---------------------------------------------------------------------------
! Setup nonlinear solver
!---------------------------------------------------------------------------
nksolver%A=>self%nlfun
nksolver%J_inv=>solver
nksolver%its=30
nksolver%atol=self%nl_tol
nksolver%rtol=1.d-20 ! Disable relative tolerance
IF(self%mfnk)THEN
  nksolver%J_update=>mfnk_update
  nksolver%up_freq=1
ELSE
  nksolver%J_update=>update_jacobian
  nksolver%up_freq=4
END IF

!---------------------------------------------------------------------------
! Setup history file
!---------------------------------------------------------------------------
IF(oft_env%head_proc)THEN
  CALL hist_file%setup('xmhd_2d.hist', desc="History file for non-linear 2D xMHD run")
  CALL hist_file%add_field('ts',   'i4', desc="Time step index")
  CALL hist_file%add_field('lits', 'i4', desc="Linear iteration count")
  CALL hist_file%add_field('nlits','i4', desc="Non-linear iteration count")
  CALL hist_file%add_field('time', 'r4', desc="Simulation time [s]")
  CALL hist_file%add_field('men',  'r4', desc="Magnetic energy [J*(2*mu0)]")
  CALL hist_file%add_field('ven',  'r4', desc="Kinetic energy [J]")
  CALL hist_file%add_field('ne',   'r4', desc="Average density [m^-3]")
  CALL hist_file%add_field('ti',   'r4', desc="Average temperature [eV]")
  CALL hist_file%add_field('stime','r4', desc="Walltime [s]")
  CALL hist_file%write_header
  CALL hist_file%open ! Open history file
END IF

!---LM-MHD disruption moments: per-step fluid observables (peaks + Joule impulse are reduced
!   offline). Only opened when the disruption source is active, so ordinary runs are untouched.
NULLIFY(vely_arr)
IF(oft_env%head_proc .AND. self%use_by_source)THEN
  OPEN(NEWUNIT=mom_unit, FILE='xmhd_2d.moments')
  WRITE(mom_unit,'(A)') '# t  Fl_net  Fl_abs  EJf_rate  Umax  Iw  EJw_rate   (MUG units; by<->bz; '// &
    'Iw=oint|by|/mu0, EJw=oint eta*by^2/(c_wall*a*mu0))'
END IF

!---------------------------------------------------------------------------
! Begin time stepping
!---------------------------------------------------------------------------
npre=0
DO i=1,self%nsteps
  IF(oft_env%head_proc)CALL mytimer%tick()
  IF(self%timestep_cn)THEN
    self%nlfun%dt=-self%dt/2.0
  ELSE
    self%nlfun%dt=0.d0
  END IF
  CALL self%nlfun%apply(u,v)
  IF(self%timestep_cn)THEN
    self%nlfun%dt=self%dt/2.0
    self%jac_dt=self%dt/2.d0
  ELSE
    self%nlfun%dt=self%dt
    self%jac_dt=self%dt
  END IF
  npre = npre + 1
  ! Update the Jacobian + preconditioner every pre_freq steps (default 1 = every step, unchanged).
  ! For stiff high-Ha ILU runs, pre_freq>1 amortizes the expensive assembly+factorization over a
  ! slowly-evolving transient, cutting per-step cost several-fold.
  IF(MOD(npre,self%pre_freq)==0)THEN
    CALL update_jacobian(u)
    CALL solver%pre%update(.TRUE.)
  END IF
  !
  DO j=maxextrap,2,-1
    CALL extrap_fields(j)%f%add(0.d0,1.d0,extrap_fields(j-1)%f)
    extrapt(j)=extrapt(j-1)
  END DO
  IF(nextrap<maxextrap)nextrap=nextrap+1
  CALL extrap_fields(1)%f%add(0.d0,1.d0,u)
  extrapt(1)=self%t
  IF(i>maxextrap)CALL vector_extrapolate(extrapt,extrap_fields,nextrap,self%t+self%dt,u)
  CALL nksolver%apply(u,v)
  IF(nksolver%cits<0)CALL oft_abort("Nonlinear solve failed","run_simulation",__FILE__)
  !---------------------------------------------------------------------------
  ! Write out initial solution progress
  !---------------------------------------------------------------------------
  n_avg=self%nlfun%diag_vals(3)/self%nlfun%diag_vals(5)
  T_avg = self%nlfun%diag_vals(4)/ self%nlfun%diag_vals(5)
  IF(oft_env%head_proc)THEN
    elapsed_time=mytimer%tock()
    hist_i4=[self%rst_base+i-1,nksolver%lits,nksolver%nlits]
    hist_r4=REAL([self%t,self%nlfun%diag_vals(1), self%nlfun%diag_vals(2), n_avg, T_avg,elapsed_time],4)
103 FORMAT(' Timestep',I8,ES14.6,2X,I4,2X,I4,F12.3,ES12.2)
    WRITE(*,103)self%rst_base+i,self%t+self%dt,nksolver%lits,nksolver%nlits,elapsed_time,self%dt
    IF(oft_debug_print(1))WRITE(*,*)
    CALL hist_file%write(data_i4=hist_i4, data_r4=hist_r4)
  END IF
  !---LM-MHD disruption moments: peak velocity + fluid observables at the new time t+dt
  !   (obs_vals reflect the just-converged residual eval). Offline: peak over t, Joule impulse = int EJf dt.
  IF(self%use_by_source)THEN
    CALL u%get_local(vely_arr,3)
    ! Peak |vely| over the central half in x (|x| < 0.5*a_half), matching box-1's Umax mask:
    ! drops the side-layer/corner spurious maxima, keeps both Hartmann layers (order=1: DOF idx = vertex idx).
    Umax=0.d0
    DO ip_umax=1,MIN(SIZE(vely_arr),mesh%np)
      IF(ABS(mesh%r(1,ip_umax)) < 0.5d0*self%a_half) Umax=MAX(Umax,ABS(vely_arr(ip_umax)))
    END DO
    IF(oft_env%head_proc) WRITE(mom_unit,'(7ES20.10)') self%t+self%dt, &
      self%nlfun%obs_vals(1), self%nlfun%obs_vals(2), self%nlfun%obs_vals(3), Umax, &
      self%nlfun%obs_vals(4), self%nlfun%obs_vals(5)
    DEALLOCATE(vely_arr); NULLIFY(vely_arr)
  END IF
  !---------------------------------------------------------------------------
  ! Update timestep and save solution
  !---------------------------------------------------------------------------
  self%t=self%t+self%dt
  CALL self%u%add(0.d0,1.d0,u)
  IF(MOD(i,self%rst_freq)==0)THEN
    IF(oft_env%head_proc)CALL mytimer%tick
    !---Create restart file
    WRITE(rst_char,104)self%rst_base+i
    READ(rst_char,104,IOSTAT=io_stat)rst_tmp
    IF((io_stat/=0).OR.(rst_tmp/=self%rst_base+i))CALL oft_abort("Step count exceeds format width", "run_simulation", __FILE__)
    CALL self%rst_save(u, self%t, self%dt, 'xmhd2d_'//rst_char//'.rst', 'U')
    IF(oft_env%head_proc)THEN
      elapsed_time=mytimer%tock()
      WRITE(*,'(2X,A,F12.3)')'I/O Time = ',elapsed_time
      CALL hist_file%flush
    END IF
    !---
  END IF
  !---
  itcount(MOD(i,XMHD_ITCACHE)+1)=nksolver%lits
  dthist(MOD(i,XMHD_ITCACHE)+1)=self%dt
  self%dt=self%ittarget*SUM(dthist)/SUM(itcount)
  IF(self%dt>dtin)self%dt=dtin
  IF(self%dt<1.d-3*dtin)CALL oft_abort('Time step dropped too low!','run_simulation',__FILE__)
  IF(ABS(1.d0-(self%dt-dthist(MOD(i,XMHD_ITCACHE)+1))/dthist(MOD(i,XMHD_ITCACHE)+1))>0.1d0)npre=-1
END DO
CALL hist_file%close()
IF(oft_env%head_proc .AND. self%use_by_source)THEN
  CLOSE(mom_unit)
  WRITE(*,'(A)') 'Wrote xmhd_2d.moments (per-step LM-MHD disruption observables)'
END IF
IF(ASSOCIATED(vely_arr))DEALLOCATE(vely_arr)
CALL nksolver%delete()
CALL solver%delete()
CALL u%delete()
CALL up%delete()
CALL v%delete()
DEALLOCATE(u,up,v)
end subroutine run_simulation
!---------------------------------------------------------------------------
!> Main driver for 2D MHD linear time advance
!---------------------------------------------------------------------------
SUBROUTINE run_lin_simulation(self)
CLASS(oft_xmhd_2d_sim), TARGET, INTENT(INOUT) :: self
!---Solver objects
CLASS(oft_vector), POINTER :: u, v, up,un1,un2
TYPE(oft_native_gmres_solver), TARGET :: solver
TYPE(OFT_TIMER) :: mytimer
!---History file fields
TYPE(oft_bin_file) :: hist_file
INTEGER(i4) :: hist_i4(3)
REAL(4) :: hist_r4(6)
!---Extrapolation fields
integer(i4), parameter :: maxextrap=2
integer(i4) :: nextrap
real(r8), allocatable, dimension(:) :: extrapt
type(oft_vector_ptr), allocatable, dimension(:) :: extrap_fields
!---
character(LEN=TDIFF_RST_LEN) :: rst_char
integer(i4) :: i,j,io_stat,rst_tmp,npre
real(r8) :: n_avg, u_avg(3), T_avg, psi_avg, by_avg,elapsed_time,diag_save(5)
real(r8), pointer :: plot_vals(:),plot_vec(:,:)
current_sim=>self
!---------------------------------------------------------------------------
! Set boundary conditions not alreadys set
!---------------------------------------------------------------------------
CALL self%setup_bc()
!---------------------------------------------------------------------------
! Normalize wall-drag field direction (drag projects out the component of
! velocity perpendicular to drag_bhat, so drag_bhat must be a unit vector).
!---------------------------------------------------------------------------
IF(self%use_wall_drag)THEN
  IF(NORM2(self%drag_bhat)<=1.d-12)THEN
    CALL oft_abort('use_wall_drag=.TRUE. but drag_bhat is a zero vector', &
      'run_simulation',__FILE__)
  END IF
  self%drag_bhat=self%drag_bhat/NORM2(self%drag_bhat)
  !--- Resolve drag_coeff from the geometry/conductance-aware Hartmann-layer closure
  !    (see drag_closure above). alpha_eff = (nu/a^2)*f(Ha,c); 'raw' leaves the
  !    user-supplied drag_coeff untouched. The 'channel' closure over-brakes a duct
  !    by ~Ha and is provided only for reference; use 'duct_*' for ducts/blankets.
  IF(TRIM(self%drag_closure)/='raw' .AND. (self%nu<=0.d0 .OR. self%a_half<=0.d0)) &
    CALL oft_abort('drag_closure /= raw needs nu>0 and a_half>0','run_simulation',__FILE__)
  SELECT CASE (TRIM(self%drag_closure))
  CASE ('raw')
    ! keep the user-supplied drag_coeff
  CASE ('channel')
    self%drag_coeff = (self%nu/self%a_half**2)*self%drag_Ha**2
  CASE ('duct_insulating')
    self%drag_coeff = (self%nu/self%a_half**2)*self%drag_Ha
  CASE ('duct_conducting')
    IF(self%drag_c<=0.d0) CALL oft_abort( &
      "drag_closure='duct_conducting' requires drag_c>0",'run_simulation',__FILE__)
    self%drag_coeff = (self%nu/self%a_half**2)*self%drag_Ha**2/(1.d0+1.d0/self%drag_c)
  CASE DEFAULT
    CALL oft_abort("unknown drag_closure (raw|channel|duct_insulating|duct_conducting)", &
      'run_simulation',__FILE__)
  END SELECT
END IF
!---------------------------------------------------------------------------
! Create solver fields
!---------------------------------------------------------------------------
IF(self%cyl_flag)CALL oft_abort( &
        "Cylindrical mode not currently supported for linear solve","run_lin_simulation",__FILE__)


call self%fe_rep%vec_create(u)
call self%fe_rep%vec_create(up)
call self%fe_rep%vec_create(v)
!---

ALLOCATE(extrap_fields(maxextrap),extrapt(maxextrap))
DO i=1,maxextrap
  CALL self%fe_rep%vec_create(extrap_fields(i)%f)
  extrapt(i)=0.d0
END DO
nextrap=0

self%t=0.d0
CALL u%add(0.d0,1.d0,self%u)
!---Create initial conditions restart file
104 FORMAT (I TDIFF_RST_LEN.TDIFF_RST_LEN)
WRITE(rst_char,104)0
CALL self%rst_save(self%u0, self%t, self%dt, 'xmhd2d_eq.rst', 'U')
CALL self%rst_save(u, self%t, self%dt, 'xmhd2d_'//rst_char//'.rst', 'U')

ALLOCATE(self%mfun)
self%mfun%parent_sim=>self

! Construct the linear advance matrix with equilibrium fields
self%jac_dt=self%dt
CALL build_approx_jacobian(self,self%u0)
!---------------------------------------------------------------------------
! Setup linear solver
!---------------------------------------------------------------------------
solver%A=>self%jacobian
solver%its=400
solver%atol=self%lin_tol
solver%itplot=1
solver%nrits=20
solver%pm=oft_env%pm
NULLIFY(solver%pre)
IF(ASSOCIATED(self%xml_pre_def))THEN
  CALL create_solver_xml(solver%pre,self%xml_pre_def)
ELSE
  CALL create_diag_pre(solver%pre)
END IF
solver%pre%A=>self%jacobian

!---------------------------------------------------------------------------
! Setup history file
!---------------------------------------------------------------------------
IF(oft_env%head_proc)THEN
  CALL hist_file%setup('xmhd_2d.hist', desc="History file for linear 2D xMHD run")
  CALL hist_file%add_field('ts',   'i4', desc="Time step index")
  CALL hist_file%add_field('lits', 'i4', desc="Linear iteration count")
  CALL hist_file%add_field('nlits','i4', desc="Non-linear iteration count")
  CALL hist_file%add_field('time', 'r4', desc="Simulation time [s]")
  CALL hist_file%add_field('men',  'r4', desc="Magnetic energy [J*(2*mu0)]")
  CALL hist_file%add_field('ven',  'r4', desc="Kinetic energy [J]")
  CALL hist_file%add_field('ne',   'r4', desc="Average density [m^-3]")
  CALL hist_file%add_field('ti',   'r4', desc="Average temperature [eV]")
  CALL hist_file%add_field('stime','r4', desc="Walltime [s]")
  CALL hist_file%write_header
  CALL hist_file%open ! Open history file
END IF
!------------------------------------------------------------------------------
! Begin time stepping
!------------------------------------------------------------------------------
IF(self%timestep_cn)THEN
  CALL self%fe_rep%vec_create(un1)
  CALL self%fe_rep%vec_create(un2)
  CALL un2%add(0.d0,1.d0,u)
  CALL un1%add(0.d0,1.d0,u)
END IF
DO i=1,self%nsteps
  IF(oft_env%head_proc)CALL mytimer%tick
  DO j=maxextrap,2,-1
    CALL extrap_fields(j)%f%add(0.d0,1.d0,extrap_fields(j-1)%f)
    extrapt(j)=extrapt(j-1)
  END DO
  IF(nextrap<maxextrap)nextrap=nextrap+1
  CALL extrap_fields(1)%f%add(0.d0,1.d0,u)
  extrapt(1)=self%t
  IF(self%timestep_cn)THEN
    CALL un2%add(0.d0,1.d0,un1)
    CALL un1%add(0.d0,1.d0,u)
  END IF
  IF(self%timestep_cn.AND.(i>1))THEN
    IF(i==2)THEN
      self%jac_dt=self%dt*2.d0/3.d0
      CALL build_approx_jacobian(self,self%u0)
      CALL solver%pre%update(.TRUE.)
    END IF
    CALL un2%add(-1.d0,4.d0,un1)
    CALL un2%scale(1.d0/3.d0)
    CALL self%mfun%apply(u,v) ! This double call and diag swap should be replaced with a better approach
    diag_save=self%mfun%diag_vals
    CALL self%mfun%apply(un2,v)
    self%mfun%diag_vals=diag_save
  ELSE
    CALL self%mfun%apply(u,v)
  END IF
  IF(i>maxextrap)CALL vector_extrapolate(extrapt,extrap_fields,nextrap,self%t+self%dt,u)
  CALL solver%apply(u,v)
  !---------------------------------------------------------------------------
  ! Write out initial solution progress
  !---------------------------------------------------------------------------
  n_avg=self%mfun%diag_vals(3)/self%mfun%diag_vals(5)
  T_avg = self%mfun%diag_vals(4)/ self%mfun%diag_vals(5)
  IF(oft_env%head_proc)THEN
    elapsed_time=mytimer%tock()
    hist_i4=(/self%rst_base+i-1,solver%cits,0/)
    hist_r4=REAL([self%t,self%mfun%diag_vals(1), self%mfun%diag_vals(2), n_avg, T_avg,elapsed_time],4)
103 FORMAT(' Timestep',I8,ES14.6,2X,I4,F12.3,ES12.2)
    WRITE(*,103)self%rst_base+i,self%t,solver%cits,elapsed_time,self%dt
    IF(oft_debug_print(1))WRITE(*,*)
    CALL hist_file%write(data_i4=hist_i4, data_r4=hist_r4)
  END IF
  !---------------------------------------------------------------------------
  ! Update timestep and save solution
  !---------------------------------------------------------------------------
  self%t=self%t+self%dt
  CALL self%u%add(0.d0,1.d0,u)
  IF(MOD(i,self%rst_freq)==0)THEN
    IF(oft_env%head_proc)CALL mytimer%tick
    !---Create restart file
    WRITE(rst_char,104)self%rst_base+i
    READ(rst_char,104,IOSTAT=io_stat)rst_tmp
    IF((io_stat/=0).OR.(rst_tmp/=self%rst_base+i))CALL oft_abort("Step count exceeds format width", "run_simulation", __FILE__)
    CALL self%rst_save(u, self%t, self%dt, 'xmhd2d_'//rst_char//'.rst', 'U')
    IF(oft_env%head_proc)THEN
      elapsed_time=mytimer%tock()
      WRITE(*,'(2X,A,F12.3)')'I/O Time = ',elapsed_time
      CALL hist_file%flush
    END IF
    !---
  END IF
END DO
CALL hist_file%close()
CALL solver%delete()
CALL u%delete()
CALL up%delete()
CALL v%delete()
DEALLOCATE(u,up,v)
IF(self%timestep_cn)THEN
  CALL un1%delete()
  CALL un2%delete()
  DEALLOCATE(un1,un2)
END IF
END SUBROUTINE run_lin_simulation
!---------------------------------------------------------------------------
!> Compute the mass matrix for RHS
!---------------------------------------------------------------------------
SUBROUTINE mfun_apply(self,a,b)
class(xmhd_2d_mfun), intent(inout) :: self !< mass function object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
type(oft_quad_type), pointer :: quad 
LOGICAL :: cyl_flag 
INTEGER(i4) :: i
REAL(r8) :: diag_vals(5), B_0(3), gamma
REAL(r8), POINTER, DIMENSION(:) :: n_weights,T_weights,psi_weights,by_weights, T_res, &
                              n_res, psi_res, by_res, vtmp, velx_res, vely_res, velz_res
REAL(r8), POINTER, DIMENSION(:,:) :: vel_weights
quad=>oft_blagrange%quad
NULLIFY(n_weights, vel_weights, T_weights, psi_weights, by_weights, &
n_res, velx_res, vely_res, velz_res, T_res, psi_res, by_res)

!---Get weights from solution vector
ALLOCATE(vel_weights(3,oft_blagrange%ne))
CALL a%get_local(n_weights,1)
vtmp => vel_weights(1, :)
CALL a%get_local(vtmp ,2)
vtmp => vel_weights(2, :)
CALL a%get_local(vtmp ,3)
vtmp => vel_weights(3, :)
CALL a%get_local(vtmp, 4)
CALL a%get_local(T_weights,5)
CALL a%get_local(psi_weights,6)
CALL a%get_local(by_weights,7)

!>--Set constant values
B_0 = self%parent_sim%B_0
cyl_flag = self%parent_sim%cyl_flag
gamma = self%parent_sim%gamma

!---Zero result and get storage array
CALL b%set(0.d0)
CALL b%get_local(n_res, 1)
CALL b%get_local(velx_res, 2)
CALL b%get_local(vely_res, 3)
CALL b%get_local(velz_res, 4)
CALL b%get_local(T_res, 5)
CALL b%get_local(psi_res, 6)
CALL b%get_local(by_res, 7)
diag_vals=0.d0

BLOCK
INTEGER(i4) :: m, jr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs
LOGICAL :: curved
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals,T_weights_loc,n_weights_loc, &
                     psi_weights_loc, by_weights_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: vel_weights_loc, basis_grads,res_loc
REAL(r8) :: jac_det,int_factor, btmp(3), tmp1(3), coords(3)
REAL(r8) :: n, vel(3), T, psi, by, dT(3),dn(3),dpsi(3),dby(3),&
         dvel(3,3),div_vel,jac_mat(3,4)
!$omp parallel private(m,jr,cell_dofs,curved,basis_vals,T_weights_loc,n_weights_loc, &
!$omp psi_weights_loc, by_weights_loc,vel_weights_loc,basis_grads,res_loc, &
!$omp jac_det,int_factor,btmp,tmp1,coords,n,vel,T,psi,by,dT,dn,dpsi,dby,dvel,div_vel,jac_mat) reduction(+:diag_vals)
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(T_weights_loc(oft_blagrange%nce),n_weights_loc(oft_blagrange%nce),&
        psi_weights_loc(oft_blagrange%nce), by_weights_loc(oft_blagrange%nce),&
        vel_weights_loc(3, oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce),res_loc(oft_blagrange%nce,7))
!$omp do ordered
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  res_loc = 0.d0 ! Zero local (cell) contribution to function
  n_weights_loc = n_weights(cell_dofs)
  vel_weights_loc = vel_weights(:, cell_dofs)
  T_weights_loc = T_weights(cell_dofs)
  psi_weights_loc = psi_weights(cell_dofs)
  by_weights_loc = by_weights(cell_dofs)
  !---------------------------------------------------------------------------
  ! Quadrature Loop
  !---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,oft_blagrange%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange,i,jr,quad%pts(:,m),basis_vals(jr))
      CALL oft_blag_geval(oft_blagrange,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    !--Extract spatial coordinates at current point
    coords = mesh%log2phys(i,quad%pts(:,m))
    !---Reconstruct values of solution fields
    n = 0.d0; dn = 0.d0; vel = 0.d0; dvel = 0.d0
    T = 0.d0; dT = 0.d0; psi = 0.d0; dpsi=0.d0
    by = 0.d0; dby = 0.d0
    basis_grads(3, :) = basis_grads(2,:)
    basis_grads(2,:) = 0.d0
    int_factor = jac_det*quad%wts(m)
    DO jr=1,oft_blagrange%nce
      n = n + n_weights_loc(jr)*basis_vals(jr)
      vel = vel + vel_weights_loc(:, jr)*basis_vals(jr)
      T = T + T_weights_loc(jr)*basis_vals(jr)
      psi = psi + psi_weights_loc(jr)*basis_vals(jr)
      by = by + by_weights_loc(jr)*basis_vals(jr)
      dn = dn + n_weights_loc(jr)*basis_grads(:,jr)
      !Note: in cylindrical coordinates, dvel is not the actual gradient of the velocity vector,
      ! but the gradient of each component, in a matrix
      dvel(:, 1) = dvel(:, 1) + vel_weights_loc(:, jr)*basis_grads(1, jr)
      dvel(:, 2) = 0.d0
      dvel(:, 3) = dvel(:, 3) + vel_weights_loc(:, jr)*basis_grads(3, jr)
      dT = dT + T_weights_loc(jr)*basis_grads(:,jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads(:,jr)
      dby = dby + by_weights_loc(jr)*basis_grads(:,jr)
    END DO
    n = n * self%parent_sim%den_scale
    dn = dn * self%parent_sim%den_scale
    div_vel = dvel(1,1) + dvel(3,3)
    btmp = cross_product(dpsi, [0.d0,1.d0,0.d0]) + by*[0.d0,1.d0,0.d0] + B_0

    diag_vals(1) = diag_vals(1) + DOT_PRODUCT(btmp,btmp)*int_factor
    diag_vals(2) = diag_vals(2) + 0.5d0*self%parent_sim%m_i*n*DOT_PRODUCT(vel,vel)*int_factor
    diag_vals(3) = diag_vals(3) + n*int_factor
    diag_vals(4) = diag_vals(4) + T*int_factor
    diag_vals(5) = diag_vals(5) + int_factor !total volume
    !---Compute local function contributions
    DO jr=1,oft_blagrange%nce
      !---Diffusion
      res_loc(jr,1) = res_loc(jr, 1) &
        + basis_vals(jr)*n*int_factor
      !---Momentum
      res_loc(jr, 2:4) = res_loc(jr, 2:4) &
      + basis_vals(jr)*vel*int_factor
      !---Temperature
      res_loc(jr,5) = res_loc(jr, 5) &
        + basis_vals(jr)*T*int_factor/(gamma-1)
      !---Induction
      res_loc(jr, 6) = res_loc(jr, 6) &
        + basis_vals(jr)*psi*int_factor 
      res_loc(jr, 7) = res_loc(jr, 7) &
      + basis_vals(jr)*by*int_factor 
    END DO
  END DO
  !---Add local values to full vector
  !$omp ordered
  DO jr=1,oft_blagrange%nce
    !$omp atomic
    n_res(cell_dofs(jr)) = n_res(cell_dofs(jr)) + res_loc(jr,1)/self%parent_sim%den_scale
    !$omp atomic
    velx_res(cell_dofs(jr)) = velx_res(cell_dofs(jr)) + res_loc(jr,2)
    !$omp atomic
    vely_res(cell_dofs(jr)) = vely_res(cell_dofs(jr)) + res_loc(jr,3)
    !$omp atomic
    velz_res(cell_dofs(jr)) = velz_res(cell_dofs(jr)) + res_loc(jr,4)
    !$omp atomic
    T_res(cell_dofs(jr)) = T_res(cell_dofs(jr)) + res_loc(jr,5)
    !$omp atomic
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + res_loc(jr,6)
    !$omp atomic
    by_res(cell_dofs(jr)) = by_res(cell_dofs(jr)) + res_loc(jr,7)
  END DO
  !$omp end ordered
END DO
!---Cleanup thread-local storage
DEALLOCATE(basis_vals,basis_grads,n_weights_loc,T_weights_loc,&
          vel_weights_loc, psi_weights_loc, by_weights_loc,cell_dofs,res_loc)
!$omp end parallel
END BLOCK
IF(oft_debug_print(2))write(*,'(4X,A)')'Applying BCs'
CALL fem_dirichlet_vec(oft_blagrange,n_weights,n_res,self%parent_sim%n_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(1, :),velx_res,self%parent_sim%velx_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(2, :),vely_res,self%parent_sim%vely_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(3, :),velz_res,self%parent_sim%velz_bc)
CALL fem_dirichlet_vec(oft_blagrange,T_weights,T_res,self%parent_sim%T_bc)
CALL fem_dirichlet_vec(oft_blagrange,psi_weights,psi_res,self%parent_sim%psi_bc)
CALL fem_dirichlet_vec(oft_blagrange,by_weights,by_res,self%parent_sim%by_bc)

!---Put results into full vector
CALL b%restore_local(n_res,1,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(velx_res,2,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(vely_res,3,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(velz_res,4,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(T_res,5,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(psi_res,6,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(by_res,7,add=.TRUE.)
self%diag_vals=oft_mpi_sum(diag_vals,5)
!---Cleanup remaining storage
DEALLOCATE(n_res,velx_res,vely_res, velz_res, T_res, psi_res, by_res, &
        n_weights,vel_weights, T_weights, psi_weights, by_weights)
end subroutine mfun_apply
!---------------------------------------------------------------------------
!> Compute the NL error function, where we are solving F(x) = 0
!!
!! b = F(a)
!---------------------------------------------------------------------------
subroutine nlfun_apply(self,a,b)
class(xmhd_2d_nlfun), intent(inout) :: self !< NL function object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
type(oft_quad_type), pointer :: quad
LOGICAL :: linear,cyl_flag
INTEGER(i4) :: i,l
REAL(r8) :: k_boltz = elec_charge
REAL(r8) :: m_i=proton_mass
REAL(r8) :: chi, eta, nu, D_diff, gamma, diag_vals(5), B_0(3), diag_vec(3)
REAL(r8) :: obs_vals(5) !< LM-MHD observables accumulator (Fl_net, Fl_abs, EJf_rate, Iw, EJw_rate)
LOGICAL :: use_wall_drag
REAL(r8) :: drag_coeff, drag_bhat(3)
LOGICAL :: use_body_force
REAL(r8) :: body_force(3)
LOGICAL :: use_by_source
REAL(r8) :: by_source_dB, by_source_tauq, S_by, t_eval
REAL(r8), POINTER, DIMENSION(:) :: n_weights,T_weights,psi_weights,by_weights, T_res, &
                              n_res, psi_res, by_res, vtmp, velx_res, vely_res, velz_res
REAL(r8), POINTER, DIMENSION(:,:) :: vel_weights
quad=>oft_blagrange%quad
NULLIFY(n_weights, vel_weights, T_weights, psi_weights, by_weights, &
n_res, velx_res, vely_res, velz_res, T_res, psi_res, by_res)

!---Get weights from solution vector
ALLOCATE(vel_weights(3,oft_blagrange%ne))
CALL a%get_local(n_weights,1)
vtmp => vel_weights(1, :)
CALL a%get_local(vtmp ,2)
vtmp => vel_weights(2, :)
CALL a%get_local(vtmp ,3)
vtmp => vel_weights(3, :)
CALL a%get_local(vtmp, 4)
CALL a%get_local(T_weights,5)
CALL a%get_local(psi_weights,6)
CALL a%get_local(by_weights,7)

!--- Set constant values
chi = self%parent_sim%chi
m_i = self%parent_sim%m_i
eta = self%parent_sim%eta
nu = self%parent_sim%nu
gamma = self%parent_sim%gamma
D_diff = self%parent_sim%D_diff
B_0 = self%parent_sim%B_0
cyl_flag = self%parent_sim%cyl_flag
use_wall_drag = self%parent_sim%use_wall_drag
drag_coeff    = self%parent_sim%drag_coeff
drag_bhat     = self%parent_sim%drag_bhat
use_body_force = self%parent_sim%use_body_force
body_force     = self%parent_sim%body_force
!--- Disruption/ELM source S(t) = -dB0z/dt in the by induction (uniform in space; evaluated
!    at the backward-Euler NEW time t+dt). Computed ONCE per residual call (time-only).
use_by_source  = self%parent_sim%use_by_source
by_source_dB   = self%parent_sim%by_source_dB
by_source_tauq = self%parent_sim%by_source_tauq
S_by = 0.d0
IF(use_by_source .AND. by_source_tauq > 0.d0)THEN
  t_eval = self%parent_sim%t + self%dt
  IF(t_eval < by_source_tauq) &
    S_by = (by_source_dB/by_source_tauq)*MIN(t_eval/(by_source_tauq/50.d0), 1.d0)
END IF

!---Zero result and get storage array
CALL b%set(0.d0)
CALL b%get_local(n_res, 1)
CALL b%get_local(velx_res, 2)
CALL b%get_local(vely_res, 3)
CALL b%get_local(velz_res, 4)
CALL b%get_local(T_res, 5)
CALL b%get_local(psi_res, 6)
CALL b%get_local(by_res, 7)
diag_vals=0.d0
obs_vals=0.d0

!!$omp parallel reduction(+:diag_vals)
BLOCK
LOGICAL :: curved
INTEGER(i4) :: k,m,jr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs
REAL(r8) :: n,vel(3),T,psi,by,dT(3),dn(3),dpsi(3),dby(3)
REAL(r8) :: dvel(3,3),div_vel,jac_mat(3,4),jac_det,int_factor,btmp(3),tmp1(3),coords(3)
REAL(r8) :: vdotb_drag !< projection of velocity onto drag_bhat (for wall drag)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals,T_weights_loc,n_weights_loc,psi_weights_loc,by_weights_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: vel_weights_loc,basis_grads,res_loc
!$omp parallel private(k,m,jr,curved,coords,cell_dofs,basis_vals,basis_grads,T_weights_loc, &
!$omp n_weights_loc,psi_weights_loc, by_weights_loc,vel_weights_loc,res_loc,jac_mat, &
!$omp jac_det,int_factor,T,n,psi,by,vel,dT,dn,dpsi,dby,dvel,div_vel,btmp,tmp1,vdotb_drag) reduction(+:diag_vals,obs_vals)
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(T_weights_loc(oft_blagrange%nce),n_weights_loc(oft_blagrange%nce),&
        psi_weights_loc(oft_blagrange%nce), by_weights_loc(oft_blagrange%nce),&
        vel_weights_loc(3, oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce),res_loc(oft_blagrange%nce,7))
!$omp do ordered
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  res_loc = 0.d0 ! Zero local (cell) contribution to function
  n_weights_loc = n_weights(cell_dofs)
  vel_weights_loc = vel_weights(:, cell_dofs)
  T_weights_loc = T_weights(cell_dofs)
  psi_weights_loc = psi_weights(cell_dofs)
  by_weights_loc = by_weights(cell_dofs)
  !---------------------------------------------------------------------------
  ! Quadrature Loop
  !---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,oft_blagrange%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange,i,jr,quad%pts(:,m),basis_vals(jr))
      CALL oft_blag_geval(oft_blagrange,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    !--Extract spatial coordinates at current point
    coords = mesh%log2phys(i,quad%pts(:,m))
    !---Reconstruct values of solution fields
    n = 0.d0; dn = 0.d0; vel = 0.d0; dvel = 0.d0
    T = 0.d0; dT = 0.d0; psi = 0.d0; dpsi=0.d0
    by = 0.d0; dby = 0.d0
    basis_grads(3, :) = basis_grads(2,:)
    basis_grads(2,:) = 0.d0
    int_factor = jac_det*quad%wts(m)
    DO jr=1,oft_blagrange%nce
      n = n + n_weights_loc(jr)*basis_vals(jr)
      vel = vel + vel_weights_loc(:, jr)*basis_vals(jr)
      T = T + T_weights_loc(jr)*basis_vals(jr)
      psi = psi + psi_weights_loc(jr)*basis_vals(jr)
      by = by + by_weights_loc(jr)*basis_vals(jr)
      dn = dn + n_weights_loc(jr)*basis_grads(:,jr)
      !Note: in cylindrical coordinates, dvel is not the actual gradient of the velocity vector,
      ! but the gradient of each component, in a matrix
      dvel(:, 1) = dvel(:, 1) + vel_weights_loc(:, jr)*basis_grads(1, jr)
      dvel(:, 2) = 0.d0
      dvel(:, 3) = dvel(:, 3) + vel_weights_loc(:, jr)*basis_grads(3, jr)
      dT = dT + T_weights_loc(jr)*basis_grads(:,jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads(:,jr)
      dby = dby + by_weights_loc(jr)*basis_grads(:,jr)
    END DO
    n = n * self%parent_sim%den_scale
    dn = dn * self%parent_sim%den_scale
    div_vel = dvel(1,1) + dvel(3,3)
    btmp = cross_product(dpsi, [0.d0,1.d0,0.d0]) + by*[0.d0,1.d0,0.d0] + B_0
    IF (cyl_flag) THEN
      div_vel = dvel(1,1) +vel(1)/(coords(1)+gs_epsilon) + dvel(3,3)
      btmp = cross_product(dpsi/(coords(1)+gs_epsilon), [0.d0,1.d0,0.d0]) + by*[0.d0,1.d0,0.d0]/(coords(1)+gs_epsilon) + B_0
    END IF
    IF (cyl_flag) THEN
      diag_vals(1) = diag_vals(1) + DOT_PRODUCT(btmp,btmp)*int_factor*coords(1)
      diag_vals(2) = diag_vals(2) + 0.5d0*m_i*n*DOT_PRODUCT(vel,vel)*int_factor*coords(1)
      diag_vals(3) = diag_vals(3) + n*int_factor*coords(1)
      diag_vals(4) = diag_vals(4) + T*int_factor*coords(1)
      diag_vals(5) = diag_vals(5) + int_factor*coords(1) !total volume
    ELSE
      diag_vals(1) = diag_vals(1) + DOT_PRODUCT(btmp,btmp)*int_factor
      diag_vals(2) = diag_vals(2) + 0.5d0*m_i*n*DOT_PRODUCT(vel,vel)*int_factor
      diag_vals(3) = diag_vals(3) + n*int_factor
      diag_vals(4) = diag_vals(4) + T*int_factor
      diag_vals(5) = diag_vals(5) + int_factor !total volume
    END IF
    !---LM-MHD fluid observables (Cartesian disruption): MUG by<->box-1 bz, dby(3)=d(by)/dz<->bzy,
    !   B_0(3)=Hartmann B0. Physical Lorentz-force density (J x B)_axial and Joule rate eta_res*J^2
    !   (eta is the diffusivity = eta_res/mu0, so eta_res*J^2 = eta*|grad by|^2/mu0). Cartesian only.
    obs_vals(1) = obs_vals(1) + (B_0(3)/mu0)*dby(3)*int_factor                        ! Fl_net
    obs_vals(2) = obs_vals(2) + ABS((B_0(3)/mu0)*dby(3))*int_factor                   ! Fl_abs
    obs_vals(3) = obs_vals(3) + eta*(dby(1)*dby(1)+dby(3)*dby(3))/mu0*int_factor      ! EJf_rate
    !---Precompute wall-drag projection (independent of test function index jr)
    IF(use_wall_drag) vdotb_drag = DOT_PRODUCT(vel, drag_bhat)
    !---Compute local function contributions
    DO jr=1,oft_blagrange%nce
      !---Diffusion
      IF (cyl_flag) THEN
        res_loc(jr,1) = res_loc(jr, 1) &
          + basis_vals(jr)*n*int_factor*coords(1) &
          +self%dt*basis_vals(jr)*DOT_PRODUCT(dn, vel)*int_factor*coords(1) &
          +self%dt*basis_vals(jr)*n*div_vel*int_factor*coords(1) &
          +self%dt*D_diff*DOT_PRODUCT(dn,basis_grads(:,jr))*int_factor*coords(1)
      ELSE
        res_loc(jr,1) = res_loc(jr, 1) &
          + basis_vals(jr)*n*int_factor &
          +self%dt*basis_vals(jr)*DOT_PRODUCT(dn, vel)*int_factor &
          +self%dt*basis_vals(jr)*n*div_vel*int_factor &
          +self%dt*D_diff*DOT_PRODUCT(dn,basis_grads(:,jr))*int_factor
      END IF
      !---Momentum
      IF (cyl_flag) THEN
        res_loc(jr, 2:4) = res_loc(jr, 2:4) &
          + basis_vals(jr)*vel*int_factor*coords(1)&
          + self%dt*DOT_PRODUCT(btmp,basis_grads(:,jr))*btmp*int_factor*coords(1)/(mu0*m_i*n) & 
          - self%dt*DOT_PRODUCT(btmp,btmp)*basis_grads(:,jr)*int_factor*coords(1)/(2*mu0*m_i*n) & 
          - self%dt*basis_vals(jr)*DOT_PRODUCT(dn,btmp)*btmp*int_factor*coords(1)/(mu0*m_i*n**2) &
          + self%dt*basis_vals(jr)*DOT_PRODUCT(btmp,btmp)*dn*int_factor*coords(1)/(2*mu0*m_i*n**2) &
          + self%dt*basis_vals(jr)*2.d0*k_boltz*dT*int_factor*coords(1)/m_i &
          + self%dt*basis_vals(jr)*2.d0*k_boltz*T*dn*int_factor*coords(1)/(m_i*n)
        DO k=1,3
          res_loc(jr,k+1) = res_loc(jr, k+1) &
            + basis_vals(jr)*self%dt*DOT_PRODUCT(vel,dvel(k,:))*int_factor*coords(1) &
            + nu*self%dt*DOT_PRODUCT(basis_grads(:,jr),dvel(k,:))*int_factor*coords(1) &
            - nu*basis_vals(jr)*self%dt*DOT_PRODUCT(dn,dvel(k,:))*int_factor*coords(1)/n
        END DO 
        res_loc(jr,2) = res_loc(jr,2) &
        + self%dt*basis_vals(jr)*(btmp(2)**2-btmp(1)**2-btmp(3)**2)*int_factor/(2.d0*mu0*m_i*n) &
        + self%dt*basis_vals(jr)*nu*vel(1)*int_factor/(coords(1)+gs_epsilon)  &
        - self%dt*basis_vals(jr)*vel(2)**2*int_factor
        res_loc(jr,3) = res_loc(jr,3) &
        - self%dt*basis_vals(jr)*btmp(1)*btmp(2)*int_factor/(mu0*m_i*n) & !!nonzero 
        + self%dt*basis_vals(jr)*nu*vel(2)*int_factor/(coords(1)+gs_epsilon) &
        + self%dt*basis_vals(jr)*vel(1)*vel(2)*int_factor
      ELSE
        res_loc(jr, 2:4) = res_loc(jr, 2:4) &
          + basis_vals(jr)*vel*int_factor&
          + self%dt*DOT_PRODUCT(btmp,basis_grads(:,jr))*btmp*int_factor/(mu0*m_i*n) & 
          - self%dt*DOT_PRODUCT(btmp,btmp)*basis_grads(:,jr)*int_factor/(2*mu0*m_i*n) & 
          - self%dt*basis_vals(jr)*DOT_PRODUCT(dn,btmp)*btmp*int_factor/(mu0*m_i*n**2) &
          + self%dt*basis_vals(jr)*DOT_PRODUCT(btmp,btmp)*dn*int_factor/(2*mu0*m_i*n**2) &
          + self%dt*basis_vals(jr)*2.d0*k_boltz*dT*int_factor/m_i &
          + self%dt*basis_vals(jr)*2.d0*k_boltz*T*dn*int_factor/(m_i*n)
        DO k=1,3
          res_loc(jr,k+1) = res_loc(jr, k+1) &
            + basis_vals(jr)*self%dt*DOT_PRODUCT(vel,dvel(k,:))*int_factor &
            + nu*self%dt*DOT_PRODUCT(basis_grads(:,jr),dvel(k,:))*int_factor &
            - nu*basis_vals(jr)*self%dt*DOT_PRODUCT(dn,dvel(k,:))*int_factor/n
        END DO
      END IF
      !---Reduced Hartmann/wall-drag source: S_u = -alpha_drag * u_perp
      ! u_perp = u - (u.drag_bhat)*drag_bhat  (component perpendicular to applied field)
      ! Backward-Euler residual convention: F += dt * alpha * phi * u_perp  (positive = dissipation)
      ! Sign: drag is a damping term; it appears on the LHS of F(u)=F_old with positive sign,
      ! equivalent to a negative contribution to the RHS of the momentum PDE.
      IF(use_wall_drag) THEN
        DO k=1,3
          IF (cyl_flag) THEN
            res_loc(jr,k+1) = res_loc(jr,k+1) &
              + drag_coeff*self%dt*basis_vals(jr)*(vel(k)-vdotb_drag*drag_bhat(k))*int_factor*coords(1)
          ELSE
            res_loc(jr,k+1) = res_loc(jr,k+1) &
              + drag_coeff*self%dt*basis_vals(jr)*(vel(k)-vdotb_drag*drag_bhat(k))*int_factor
          END IF
        END DO
      END IF
      !---Uniform body-force momentum source: S_u = body_force (per unit mass).
      ! Same sign convention as the pressure-gradient term: a positive RHS force
      ! appears as -dt*phi*f in the residual F(u)=F_old. Constant => no Jacobian term.
      IF(use_body_force) THEN
        DO k=1,3
          IF (cyl_flag) THEN
            res_loc(jr,k+1) = res_loc(jr,k+1) &
              - self%dt*basis_vals(jr)*body_force(k)*int_factor*coords(1)
          ELSE
            res_loc(jr,k+1) = res_loc(jr,k+1) &
              - self%dt*basis_vals(jr)*body_force(k)*int_factor
          END IF
        END DO
      END IF
      !---Temperature
      IF (cyl_flag) THEN
        res_loc(jr,5) = res_loc(jr, 5) &
          + basis_vals(jr)*T*int_factor*coords(1)/(gamma-1) &
          + self%dt*basis_vals(jr)*DOT_PRODUCT(vel, dT)*int_factor*coords(1)/(gamma-1) &
          + self%dt*basis_vals(jr)*T*div_vel*int_factor*coords(1) & 
          + self%dt*chi*DOT_PRODUCT(dT, basis_grads(:,jr))*int_factor*coords(1) &
          - self%dt*chi*basis_vals(jr)*DOT_PRODUCT(dn, dT)*int_factor*coords(1)/n 
      ELSE
        res_loc(jr,5) = res_loc(jr, 5) &
          + basis_vals(jr)*T*int_factor/(gamma-1) &
          + self%dt*basis_vals(jr)*DOT_PRODUCT(vel, dT)*int_factor/(gamma-1) &
          + self%dt*basis_vals(jr)*T*div_vel*int_factor & 
          + self%dt*chi*DOT_PRODUCT(dT, basis_grads(:,jr))*int_factor &
          - self%dt*chi*basis_vals(jr)*DOT_PRODUCT(dn, dT)*int_factor/n 
      END IF
      !---Psi
      tmp1 = cross_product(B_0,vel)
      IF(cyl_flag) THEN
        res_loc(jr, 6) = res_loc(jr, 6) &
        + basis_vals(jr)*psi*int_factor/(coords(1)+gs_epsilon) &
        + basis_vals(jr)*self%dt*DOT_PRODUCT(vel, dpsi)*int_factor/(coords(1)+gs_epsilon) &
        + basis_vals(jr)*self%dt*tmp1(2)*int_factor/(coords(1)+gs_epsilon) &
        + self%dt*eta*DOT_PRODUCT(basis_grads(:,jr), dpsi)*int_factor/(coords(1)+gs_epsilon)
      ELSE
        res_loc(jr, 6) = res_loc(jr, 6) &
        + basis_vals(jr)*psi*int_factor &
        + basis_vals(jr)*self%dt*DOT_PRODUCT(vel, dpsi)*int_factor &
        + basis_vals(jr)*self%dt*tmp1(2)*int_factor &
        + self%dt*eta*DOT_PRODUCT(basis_grads(:,jr), dpsi)*int_factor
      END IF
      ! --By
      tmp1 = cross_product(dpsi,dvel(2, :))
      IF(cyl_flag) THEN
        res_loc(jr,7) = res_loc(jr,7) &
        + basis_vals(jr)*by*int_factor/(coords(1)+gs_epsilon) &
        - self%dt*basis_vals(jr)*tmp1(2)*int_factor/(coords(1)+gs_epsilon) &
        + self%dt*basis_vals(jr)*dvel(1,1)*by*int_factor/(coords(1)+gs_epsilon) &
        + self%dt*basis_vals(jr)*dvel(3,3)*by*int_factor/(coords(1)+gs_epsilon) &
        + self%dt*basis_vals(jr)*DOT_PRODUCT(vel, dby)*int_factor/(coords(1)+gs_epsilon) &
        - self%dt*basis_vals(jr)*vel(1)*by*int_factor/(coords(1)+gs_epsilon)**2 &
        + self%dt*eta*DOT_PRODUCT(basis_grads(:,jr), dby)*int_factor/(coords(1)+gs_epsilon)
      ELSE
        ! Field-line stretching source (B.grad)v_y with the FULL in-plane field
        ! B_in = grad(psi) x yhat + B_0_in ((yhat x B_0) x yhat = B_0_inplane maps B_0
        ! into the flux-gradient slot). With B_in = +grad(psi) x yhat (the momentum
        ! btmp convention), (B_in.grad)v_y = -tmp1(2), so the backward-Euler residual
        ! term is +dt*phi*tmp1(2). The original -dt*phi*tmp1(2) was a latent sign error
        ! in the (previously unexercised) by<->vely polarization: verified empirically
        ! by test_alfven_by2d (old sign: exponential growth at rate k*v_A; this sign:
        ! oscillation) and by the Shercliff duct (this sign: Hartmann braking, 0.004%
        ! vs analytic). Zero effect for purely out-of-plane (axial) B_0 with by=0
        ! initial data paths that never grow by. Cylindrical branch intentionally
        ! untouched pending a cyl-mode test.
        tmp1 = cross_product(dpsi + cross_product([0.d0,1.d0,0.d0],B_0), dvel(2,:))
        res_loc(jr, 7) = res_loc(jr, 7) &
        + basis_vals(jr)*by*int_factor &
        + basis_vals(jr)*self%dt*tmp1(2)*int_factor &
        + basis_vals(jr)*self%dt*DOT_PRODUCT(vel, dby)*int_factor &
        + basis_vals(jr)*self%dt*by*div_vel*int_factor &
        + self%dt*eta*DOT_PRODUCT(basis_grads(:,jr), dby)*int_factor &
        - self%dt*basis_vals(jr)*S_by*int_factor    ! disruption source S(t)=-dB0z/dt (uniform, time-only)
      END IF
    END DO
  END DO
  !---Add local values to full vector
  !$omp ordered
  DO jr=1,oft_blagrange%nce
    !$omp atomic
    n_res(cell_dofs(jr)) = n_res(cell_dofs(jr)) + res_loc(jr,1)/self%parent_sim%den_scale
    !$omp atomic
    velx_res(cell_dofs(jr)) = velx_res(cell_dofs(jr)) + res_loc(jr,2)
    !$omp atomic
    vely_res(cell_dofs(jr)) = vely_res(cell_dofs(jr)) + res_loc(jr,3)
    !$omp atomic
    velz_res(cell_dofs(jr)) = velz_res(cell_dofs(jr)) + res_loc(jr,4)
    !$omp atomic
    T_res(cell_dofs(jr)) = T_res(cell_dofs(jr)) + res_loc(jr,5)
    !$omp atomic
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + res_loc(jr,6)
    !$omp atomic
    by_res(cell_dofs(jr)) = by_res(cell_dofs(jr)) + res_loc(jr,7)
  END DO
  !$omp end ordered
END DO

!---Cleanup thread-local storage
DEALLOCATE(basis_vals,basis_grads,n_weights_loc,T_weights_loc,&
          vel_weights_loc, psi_weights_loc, by_weights_loc,cell_dofs,res_loc)
!$omp end parallel
END BLOCK
!!$omp end parallel
IF(oft_debug_print(2))write(*,'(4X,A)')'Applying BCs'
CALL fem_dirichlet_vec(oft_blagrange,n_weights,n_res,self%parent_sim%n_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(1, :),velx_res,self%parent_sim%velx_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(2, :),vely_res,self%parent_sim%vely_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(3, :),velz_res,self%parent_sim%velz_bc)
CALL fem_dirichlet_vec(oft_blagrange,T_weights,T_res,self%parent_sim%T_bc)
CALL fem_dirichlet_vec(oft_blagrange,psi_weights,psi_res,self%parent_sim%psi_bc)
CALL fem_dirichlet_vec(oft_blagrange,by_weights,by_res,self%parent_sim%by_bc)
!---Phase-4: thin-wall (finite wall-conductance) Robin BC on by:
!  b + c_wall*db/dn = 0 at the wall. Integrating the by diffusion term by parts
!  leaves the boundary integral -dt*eta*OINT (db/dn)*phi ds = +dt*(eta/c_wall)*
!  OINT b*phi ds. For order-1 elements the edge integral is exact via the 1D
!  linear mass matrix |e|/6*[[2,1],[1,2]] — no quadrature needed. Applied only
!  on boundary edges whose BOTH endpoints are free (not by-Dirichlet): the
!  driver's by_bc mask selects thin-wall vs insulating walls per side.
!  Cartesian, order-1 only (guarded in setup); each boundary edge contributes
!  exactly once (serial loop, after the threaded cell assembly, after the
!  Dirichlet masking so Robin rows are never Dirichlet rows).
IF(self%parent_sim%c_wall > 0.d0)THEN
  BLOCK
  INTEGER(i4) :: ke, ed, ip1, ip2
  REAL(r8) :: elen, rfac
  DO ke=1,mesh%nbe
    ed = mesh%lbe(ke)
    ip1 = mesh%le(1,ed); ip2 = mesh%le(2,ed)
    IF(self%parent_sim%by_bc(ip1).OR.self%parent_sim%by_bc(ip2))CYCLE
    elen = SQRT(SUM((mesh%r(:,ip1)-mesh%r(:,ip2))**2))
    rfac = self%dt*eta*elen/(self%parent_sim%c_wall*6.d0)
    by_res(ip1) = by_res(ip1) + rfac*(2.d0*by_weights(ip1) + by_weights(ip2))
    by_res(ip2) = by_res(ip2) + rfac*(by_weights(ip1) + 2.d0*by_weights(ip2))
  END DO
  END BLOCK
END IF
!---LM-MHD thin-wall wall observables (same Hartmann-wall edge set as the Robin BC above:
!   boundary edges whose endpoints are FREE in by_bc). Order=1 so by is linear on each edge:
!   Iw  = oint |by|/mu0 dl  ~ trapezoid on |by|;  EJw = oint eta*by^2/(c_wall*a*mu0) dl, exact int by^2.
IF(self%parent_sim%use_by_source .AND. self%parent_sim%c_wall > 0.d0)THEN
  BLOCK
  INTEGER(i4) :: ke, ed, ip1, ip2
  REAL(r8) :: elen, b1, b2, ejw_fac
  ejw_fac = eta/(self%parent_sim%c_wall*self%parent_sim%a_half*mu0)
  DO ke=1,mesh%nbe
    ed = mesh%lbe(ke)
    ip1 = mesh%le(1,ed); ip2 = mesh%le(2,ed)
    IF(self%parent_sim%by_bc(ip1).OR.self%parent_sim%by_bc(ip2))CYCLE
    elen = SQRT(SUM((mesh%r(:,ip1)-mesh%r(:,ip2))**2))
    b1 = by_weights(ip1); b2 = by_weights(ip2)
    obs_vals(4) = obs_vals(4) + 0.5d0*(ABS(b1)+ABS(b2))*elen/mu0             ! Iw = oint |by|/mu0 dl
    obs_vals(5) = obs_vals(5) + ejw_fac*(b1*b1 + b1*b2 + b2*b2)*elen/3.d0    ! EJw = oint eta*by^2/(c_wall*a*mu0) dl
  END DO
  END BLOCK
END IF
!---Put results into full vector
CALL b%restore_local(n_res,1,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(velx_res,2,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(vely_res,3,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(velz_res,4,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(T_res,5,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(psi_res,6,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(by_res,7,add=.TRUE.)
self%diag_vals=oft_mpi_sum(diag_vals,5)
self%obs_vals=oft_mpi_sum(obs_vals,5)
!---Cleanup remaining storage
DEALLOCATE(n_res,velx_res,vely_res, velz_res, T_res, psi_res, by_res, &
        n_weights,vel_weights, T_weights, psi_weights, by_weights)
END SUBROUTINE nlfun_apply
!---------------------------------------------------------------------------
!> Compute the approximate Jacobian matrix for the nonlinear function being solved
!---------------------------------------------------------------------------
subroutine build_approx_jacobian(self,a)
class(oft_xmhd_2d_sim), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Solution for computing jacobian
LOGICAL :: cyl_flag, linear
LOGICAL :: use_wall_drag
INTEGER(i4) :: i
REAL(r8) :: k_boltz=elec_charge
REAL(r8) :: m_i = proton_mass
REAL(r8) :: chi, eta, nu, D_diff, gamma, B_0(3), diag_vals(7), dt_fac
REAL(r8) :: drag_coeff, drag_bhat(3)
REAL(r8), POINTER, DIMENSION(:) :: n_weights,T_weights, psi_weights, by_weights, vtmp
REAL(r8), POINTER, DIMENSION(:,:) :: vel_weights
integer(KIND=omp_lock_kind), allocatable, dimension(:) :: tlocks
class(oft_vector), pointer :: tmp
type(oft_quad_type), pointer :: quad
quad=>oft_blagrange%quad
CALL self%jacobian%zero
NULLIFY(n_weights,vel_weights, T_weights, &
         psi_weights, by_weights, vtmp)
!---Get weights from solution vector
CALL a%get_local(n_weights,1)
ALLOCATE(vel_weights(3,oft_blagrange%ne))
vtmp => vel_weights(1, :)
CALL a%get_local(vtmp ,2)
vtmp => vel_weights(2, :)
CALL a%get_local(vtmp ,3)
vtmp => vel_weights(3, :)
CALL a%get_local(vtmp, 4)
CALL a%get_local(T_weights,5)
CALL a%get_local(psi_weights,6)
CALL a%get_local(by_weights,7)
!---
chi = self%chi
m_i = self%m_i
eta = self%eta
nu = self%nu
gamma = self%gamma
D_diff = self%D_diff
B_0 = self%B_0
cyl_flag = self%cyl_flag
linear = self%linear
dt_fac = self%jac_dt
use_wall_drag = self%use_wall_drag
drag_coeff    = self%drag_coeff
drag_bhat     = self%drag_bhat

!--Setup thread locks
ALLOCATE(tlocks(self%fe_rep%nfields))
DO i=1,self%fe_rep%nfields
  call omp_init_lock(tlocks(i))
END DO
!!$omp parallel
BLOCK
LOGICAL :: curved
INTEGER(i4) :: k, l, m, ik, jr, jc
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs
REAL(r8) :: n,vel(3),T,psi,by,dT(3),dn(3),dpsi(3),dby(3),dvel(3,3),div_vel
REAL(r8) :: jac_mat(3,4),jac_det,int_factor,btmp(3),tmp2(3),tmp3(3),coords(3)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals,n_weights_loc,T_weights_loc
REAL(r8), ALLOCATABLE, DIMENSION(:) :: psi_weights_loc,by_weights_loc,res_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: vel_weights_loc,basis_grads
TYPE(oft_1d_int), ALLOCATABLE, DIMENSION(:) :: iloc
type(oft_local_mat), allocatable, dimension(:,:) :: jac_loc
!$omp parallel private(ik, k, l, m,jr,jc,curved,coords,cell_dofs,basis_vals,basis_grads,T_weights_loc, &
!$omp n_weights_loc,psi_weights_loc, by_weights_loc,vel_weights_loc,res_loc,jac_mat, &
!$omp jac_det,int_factor,T,n,psi,by,vel,dT,dn,dpsi,dby,dvel,div_vel,btmp,tmp2,tmp3, iloc, jac_loc) reduction(+:diag_vals)
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(n_weights_loc(oft_blagrange%nce),vel_weights_loc(3, oft_blagrange%nce),&
        T_weights_loc(oft_blagrange%nce), psi_weights_loc(oft_blagrange%nce),&
        by_weights_loc(oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce))
ALLOCATE(jac_loc(self%fe_rep%nfields,self%fe_rep%nfields))
ALLOCATE(iloc(self%fe_rep%nfields))
DO ik=1,self%fe_rep%nfields
   iloc(ik)%v=>cell_dofs
END DO
CALL self%fe_rep%mat_setup_local(jac_loc, self%jacobian_block_mask)
!$omp do ordered
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  CALL self%fe_rep%mat_zero_local(jac_loc) ! Zero local (cell) contribution to matrix
  n_weights_loc = n_weights(cell_dofs)
  vel_weights_loc = vel_weights(:, cell_dofs)
  T_weights_loc = T_weights(cell_dofs)
  psi_weights_loc = psi_weights(cell_dofs)
  by_weights_loc = by_weights(cell_dofs)

!---------------------------------------------------------------------------
! Quadrature Loop
!---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,oft_blagrange%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange,i,jr,quad%pts(:,m),basis_vals(jr))
      CALL oft_blag_geval(oft_blagrange,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    !--Extract spatial coordinates at current point
    coords = mesh%log2phys(i,quad%pts(:,m))
    basis_grads(3, :) = basis_grads(2,:)
    basis_grads(2,:) = 0.d0
    int_factor = jac_det*quad%wts(m)
    !---Reconstruct values of solution fields
    n = 0.d0; dn = 0.d0; vel = 0.d0; dvel = 0.d0
    T = 0.d0; dT = 0.d0; psi = 0.d0; dpsi=0.d0
    by = 0.d0; dby = 0.d0
    DO jr=1,oft_blagrange%nce
      n = n + n_weights_loc(jr)*basis_vals(jr)
      vel = vel + vel_weights_loc(:, jr)*basis_vals(jr)
      T = T + T_weights_loc(jr)*basis_vals(jr)
      psi = psi + psi_weights_loc(jr)*basis_vals(jr)
      by = by + by_weights_loc(jr)*basis_vals(jr)
      dn = dn + n_weights_loc(jr)*basis_grads(:,jr)
      dvel(:, 1) = dvel(:, 1) + vel_weights_loc(:, jr)*basis_grads(1, jr)
      dvel(:, 2) = 0.d0
      dvel(:, 3) = dvel(:, 3) + vel_weights_loc(:, jr)*basis_grads(3, jr)
      dT = dT + T_weights_loc(jr)*basis_grads(:,jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads(:,jr)
      dby = dby + by_weights_loc(jr)*basis_grads(:,jr)
    END DO
    n = n * self%den_scale
    dn = dn * self%den_scale
    div_vel = dvel(1,1) + dvel(3,3)
    btmp = cross_product(dpsi, [0.d0,1.d0,0.d0]) + by*[0.d0,1.d0,0.d0] + B_0
    IF (cyl_flag) THEN
      div_vel = dvel(1,1) +vel(1)/(coords(1)+gs_epsilon) + dvel(3,3)
      btmp = cross_product(dpsi/(coords(1)+gs_epsilon), [0.d0,1.d0,0.d0]) + by*[0.d0,1.d0,0.d0]/(coords(1)+gs_epsilon) + B_0
    END IF
    diag_vals = diag_vals + [n, vel(1), vel(2), vel(3), T, psi, by]*int_factor
    IF (cyl_flag) THEN
      btmp = cross_product(dpsi/coords(1), [0.d0,1.d0,0.d0]) + by*[0.d0,1.d0,0.d0] + B_0
    ELSE
      btmp = cross_product(dpsi, [0.d0,1.d0,0.d0]) + by*[0.d0,1.d0,0.d0] + B_0
    END IF
    !---Compute local matrix contributions
    DO jr=1,oft_blagrange%nce
      DO jc=1,oft_blagrange%nce
        ! Diffusion
        ! --n, n
        IF (cyl_flag) THEN
          jac_loc(1, 1)%m(jr,jc) = jac_loc(1, 1)%m(jr, jc) &
          + basis_vals(jr)*basis_vals(jc)*int_factor*coords(1) &
          + dt_fac*basis_vals(jr)*DOT_PRODUCT(basis_grads(:, jc), vel)*int_factor*coords(1) & !delta_n*div(u)
          + dt_fac*basis_vals(jr)*basis_vals(jc)*div_vel*int_factor*coords(1) & ! u dot grad(delta_n)
          + dt_fac*D_diff*DOT_PRODUCT(basis_grads(:, jr),basis_grads(:, jc))*int_factor*coords(1)
        ELSE
          jac_loc(1, 1)%m(jr,jc) = jac_loc(1, 1)%m(jr, jc) &
          + basis_vals(jr)*basis_vals(jc)*int_factor &
          + dt_fac*basis_vals(jr)*DOT_PRODUCT(basis_grads(:, jc), vel)*int_factor & !delta_n*div(u)
          + dt_fac*basis_vals(jr)*basis_vals(jc)*div_vel*int_factor & ! u dot grad(delta_n)
          + dt_fac*D_diff*DOT_PRODUCT(basis_grads(:, jr),basis_grads(:, jc))*int_factor
        END IF
        ! --n, vel
        IF (cyl_flag) THEN
          DO l=1,3
            jac_loc(1, l+1)%m(jr,jc) = jac_loc(1, l+1)%m(jr, jc) &
            + basis_vals(jr)*dt_fac*n*basis_grads(l, jc)*int_factor*coords(1) & ! n*div(delta_u) = n*SUM(basis_grads)?
            + basis_vals(jr)*dt_fac*basis_vals(jc)*dn(l)*int_factor*coords(1) ! u dot grad(n)
          END DO          
          jac_loc(1,2)%m(jr,jc) = jac_loc(1,2)%m(jr,jc) &
          + basis_vals(jr)*dt_fac*n*basis_vals(jc)*int_factor
        ELSE
          DO l=1,3
            jac_loc(1, l+1)%m(jr,jc) = jac_loc(1, l+1)%m(jr, jc) &
            + basis_vals(jr)*dt_fac*n*basis_grads(l, jc)*int_factor & ! n*div(delta_u) = n*SUM(basis_grads)?
            + basis_vals(jr)*dt_fac*basis_vals(jc)*dn(l)*int_factor ! u dot grad(n)
          END DO
        END IF
        !--Momentum 
        !--vel,n
        IF (cyl_flag) THEN
          DO l=1,3
            jac_loc(l+1,1)%m(jr,jc) = jac_loc(l+1,1)%m(jr,jc) &
            ! grad T terms
            + dt_fac*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*dT(l)*int_factor*coords(1)/(m_i*n) & ! 
            + dt_fac*basis_vals(jr)*2.d0*k_boltz*T*basis_grads(l,jc)*int_factor*coords(1)/(m_i*n) &
            - dt_fac*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*n*dT(l)*int_factor*coords(1)/(m_i*n**2.d0) &
            - dt_fac*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*dn(l)*T*int_factor*coords(1)/(m_i*n**2.d0) &
            ! dyadic b terms
            - dt_fac*basis_vals(jc)*DOT_PRODUCT(basis_grads(:,jr), btmp)*btmp(l)*int_factor*coords(1)/(mu0*m_i*n**2.d0) &
            - dt_fac*basis_vals(jr)*DOT_PRODUCT(basis_grads(:,jc), btmp)*btmp(l)*int_factor*coords(1)/(mu0*m_i*n**2.d0) &
            + dt_fac*basis_vals(jr)*2.d0*basis_vals(jc)*DOT_PRODUCT(dn, btmp)*btmp(l)*int_factor*coords(1)/(mu0*m_i*n**3.d0) &
            + dt_fac*basis_vals(jc)*basis_grads(l,jr)*DOT_PRODUCT(btmp, btmp)*int_factor*coords(1)/(2.d0*mu0*m_i*n**2.d0) &
            + dt_fac*basis_vals(jr)*basis_grads(l,jc)*DOT_PRODUCT(btmp, btmp)*int_factor*coords(1)/(2.d0*mu0*m_i*n**2.d0) &
            - dt_fac*basis_vals(jr)*basis_vals(jc)*dn(l)*DOT_PRODUCT(btmp, btmp)*int_factor*coords(1)/(mu0*m_i*n**3.d0) &
            ! velocity dissipation terms
            - dt_fac*nu*basis_vals(jc)*DOT_PRODUCT(basis_grads(:,jr), dvel(l, :))*int_factor*coords(1)/n & 
            - dt_fac*nu*basis_vals(jr)*DOT_PRODUCT(basis_grads(:,jc), dvel(l, :))*int_factor*coords(1)/n &
            + dt_fac*nu*basis_vals(jr)*2.d0*basis_vals(jc)*DOT_PRODUCT(dn, dvel(l, :))*int_factor*coords(1)/n**2
          END DO       
          jac_loc(2,1)%m(jr,jc) = jac_loc(2,1)%m(jr,jc) &
          - dt_fac*basis_vals(jr)*basis_vals(jc)*(btmp(2)**2-btmp(1)**2-btmp(3)**2)*int_factor/(mu0*m_i*n**2*(coords(1)+gs_epsilon)) &
          - dt_fac*basis_vals(jr)*basis_vals(jc)*nu*vel(1)*int_factor/(n*(coords(1)+gs_epsilon))
          jac_loc(3,1)%m(jr,jc) = jac_loc(3,1)%m(jr,jc) &
          + dt_fac*basis_vals(jr)*basis_vals(jc)*(btmp(1)*btmp(2))*int_factor/(mu0*m_i*n**2*(coords(1)+gs_epsilon)) &
          - dt_fac*basis_vals(jr)*basis_vals(jc)*nu*vel(2)*int_factor/(n*(coords(1)+gs_epsilon))
        ELSE
          IF (linear) THEN
            DO l=1,3
              jac_loc(l+1,1)%m(jr,jc) = jac_loc(l+1,1)%m(jr,jc) &
              + dt_fac*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*dT(l)*int_factor/(m_i*n) & ! 
              + dt_fac*basis_vals(jr)*2.d0*k_boltz*T*basis_grads(l,jc)*int_factor/(m_i*n) &
              + basis_vals(jr)*basis_vals(jc)*dt_fac*DOT_PRODUCT(vel,dvel(l,:))*int_factor/n 
            END DO
          ELSE
            DO l=1,3
              jac_loc(l+1,1)%m(jr,jc) = jac_loc(l+1,1)%m(jr,jc) &
              ! grad T terms
              + dt_fac*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*dT(l)*int_factor/(m_i*n) & ! 
              + dt_fac*basis_vals(jr)*2.d0*k_boltz*T*basis_grads(l,jc)*int_factor/(m_i*n) &
              - dt_fac*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*n*dT(l)*int_factor/(m_i*n**2.d0) &
              - dt_fac*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*dn(l)*T*int_factor/(m_i*n**2.d0) &
              ! dyadic b terms
              - dt_fac*basis_vals(jc)*DOT_PRODUCT(basis_grads(:,jr), btmp)*btmp(l)*int_factor/(mu0*m_i*n**2.d0) &
              - dt_fac*basis_vals(jr)*DOT_PRODUCT(basis_grads(:,jc), btmp)*btmp(l)*int_factor/(mu0*m_i*n**2.d0) &
              + dt_fac*basis_vals(jr)*2.d0*basis_vals(jc)*DOT_PRODUCT(dn, btmp)*btmp(l)*int_factor/(mu0*m_i*n**3.d0) &
              + dt_fac*basis_vals(jc)*basis_grads(l,jr)*DOT_PRODUCT(btmp, btmp)*int_factor/(2.d0*mu0*m_i*n**2.d0) &
              + dt_fac*basis_vals(jr)*basis_grads(l,jc)*DOT_PRODUCT(btmp, btmp)*int_factor/(2.d0*mu0*m_i*n**2.d0) &
              - dt_fac*basis_vals(jr)*basis_vals(jc)*dn(l)*DOT_PRODUCT(btmp, btmp)*int_factor/(mu0*m_i*n**3.d0) &
              ! velocity dissipation terms
              - dt_fac*nu*basis_vals(jc)*DOT_PRODUCT(basis_grads(:,jr), dvel(l, :))*int_factor/n & 
              - dt_fac*nu*basis_vals(jr)*DOT_PRODUCT(basis_grads(:,jc), dvel(l, :))*int_factor/n &
              + dt_fac*nu*basis_vals(jr)*2.d0*basis_vals(jc)*DOT_PRODUCT(dn, dvel(l, :))*int_factor/n**2
            END DO
          END IF
        END IF
        !--vel, vel
        IF (cyl_flag) THEN
          DO k=1,3
            jac_loc(k+1,k+1)%m(jr,jc)= jac_loc(k+1,k+1)%m(jr,jc) &
              + basis_vals(jr)*basis_vals(jc)*int_factor*coords(1) &
              + dt_fac*basis_vals(jr)*DOT_PRODUCT(vel, basis_grads(:,jc))*int_factor*coords(1) &
              + dt_fac*nu*DOT_PRODUCT(basis_grads(:,jr), basis_grads(:,jc))*int_factor*coords(1) &
              - dt_fac*nu*basis_vals(jr)*DOT_PRODUCT(dn, basis_grads(:,jc))*int_factor*coords(1)/n
            DO l=1,3
              jac_loc(k+1,l+1)%m(jr,jc)= jac_loc(k+1,l+1)%m(jr,jc) &
              + dt_fac*basis_vals(jr)*basis_vals(jc)*dvel(k, l)*int_factor*coords(1)
            END DO
          END DO
          jac_loc(2,2)%m(jr,jc)= jac_loc(2,2)%m(jr,jc) &
          + dt_fac*nu*basis_vals(jr)*basis_vals(jc)*int_factor/(coords(1)+gs_epsilon)
          jac_loc(2,3)%m(jr,jc)= jac_loc(2,3)%m(jr,jc) &
          -dt_fac*basis_vals(jr)*2.d0*vel(2)*basis_vals(jc)*int_factor
          jac_loc(3,2)%m(jr,jc)= jac_loc(3,2)%m(jr,jc) &
          + dt_fac*basis_vals(jr)*basis_vals(jc)*vel(2)*int_factor
          jac_loc(3,3)%m(jr,jc)= jac_loc(3,3)%m(jr,jc) &
          + dt_fac*basis_vals(jr)*basis_vals(jc)*vel(1)*int_factor &
          + dt_fac*nu*basis_vals(jr)*basis_vals(jc)*int_factor/(coords(1)+gs_epsilon)
          !---Wall-drag Jacobian (cylindrical): same projection, R-weighted
          IF(use_wall_drag) THEN
            DO k=1,3
              DO l=1,3
                jac_loc(k+1,l+1)%m(jr,jc) = jac_loc(k+1,l+1)%m(jr,jc) &
                  + dt_fac*drag_coeff*basis_vals(jr)*basis_vals(jc) &
                    *(merge(1.d0,0.d0,k==l) - drag_bhat(k)*drag_bhat(l))*int_factor*coords(1)
              END DO
            END DO
          END IF
        ELSE
          DO k=1,3
            jac_loc(k+1,k+1)%m(jr,jc)= jac_loc(k+1,k+1)%m(jr,jc) &
              + basis_vals(jr)*basis_vals(jc)*int_factor &
              + dt_fac*basis_vals(jr)*DOT_PRODUCT(vel, basis_grads(:,jc))*int_factor &
              + dt_fac*nu*DOT_PRODUCT(basis_grads(:,jr), basis_grads(:,jc))*int_factor &
              - dt_fac*nu*basis_vals(jr)*DOT_PRODUCT(dn, basis_grads(:,jc))*int_factor/n
            DO l=1,3
              jac_loc(k+1,l+1)%m(jr,jc)= jac_loc(k+1,l+1)%m(jr,jc) &
              + dt_fac*basis_vals(jr)*basis_vals(jc)*dvel(k, l)*int_factor
            END DO
          END DO
          !---Wall-drag Jacobian: d/d(vel_l)[ alpha * phi_jr * (v_k - (v.bhat)*bhat_k) * phi_jc ]
          !   = alpha * phi_jr * phi_jc * (delta_kl - bhat_k * bhat_l)
          IF(use_wall_drag) THEN
            DO k=1,3
              DO l=1,3
                jac_loc(k+1,l+1)%m(jr,jc) = jac_loc(k+1,l+1)%m(jr,jc) &
                  + dt_fac*drag_coeff*basis_vals(jr)*basis_vals(jc) &
                    *(merge(1.d0,0.d0,k==l) - drag_bhat(k)*drag_bhat(l))*int_factor
              END DO
            END DO
          END IF
        END IF
        ! --vel, T
        IF (cyl_flag) THEN
          DO l=1,3
            jac_loc(l+1,5)%m(jr,jc) = jac_loc(l+1,5)%m(jr,jc) &
            + dt_fac*basis_vals(jr)*2*k_boltz*dn(l)*basis_vals(jc)*int_factor*coords(1)/(m_i*n) &
            + dt_fac*basis_vals(jr)*2*k_boltz*basis_grads(l,jc)*int_factor*coords(1)/(m_i) 
          END DO
        ELSE
          DO l=1,3
            jac_loc(l+1,5)%m(jr,jc) = jac_loc(l+1,5)%m(jr,jc) &
            + dt_fac*basis_vals(jr)*2*k_boltz*dn(l)*basis_vals(jc)*int_factor/(m_i*n) &
            + dt_fac*basis_vals(jr)*2*k_boltz*basis_grads(l,jc)*int_factor/(m_i) 
          END DO
        END IF
        ! --vel, psi
        IF (cyl_flag) THEN
          tmp2 = cross_product(basis_grads(:,jc), [0.d0,1.d0/(coords(1)+gs_epsilon),0.d0]) ! this is 'dB'
          DO l=1,3
            jac_loc(l+1,6)%m(jr,jc) = jac_loc(l+1,6)%m(jr,jc) &
            + dt_fac*DOT_PRODUCT(basis_grads(:,jr), tmp2)*btmp(l)*int_factor*coords(1)/(m_i*n*mu0) &
            + dt_fac*DOT_PRODUCT(basis_grads(:,jr), btmp)*tmp2(l)*int_factor*coords(1)/(m_i*n*mu0) &
            - dt_fac*basis_grads(l,jr)*DOT_PRODUCT(tmp2, btmp)*int_factor*coords(1)/(m_i*n*mu0) &
            - dt_fac*basis_vals(jr)*DOT_PRODUCT(dn, tmp2)*btmp(l)*int_factor*coords(1)/(m_i*n**2*mu0) &
            - dt_fac*basis_vals(jr)*DOT_PRODUCT(dn, btmp)*tmp2(l)*int_factor*coords(1)/(m_i*n**2*mu0) &
            + dt_fac*basis_vals(jr)*dn(l)*DOT_PRODUCT(tmp2, btmp)*int_factor*coords(1)/(m_i*n**2*mu0) 
          END DO
          jac_loc(2,6)%m(jr,jc) = jac_loc(2,6)%m(jr,jc) &
          + dt_fac*basis_vals(jr)*(btmp(2)*tmp2(2)-btmp(1)*tmp2(1)-btmp(3)*tmp2(3))*int_factor/(m_i*n*mu0)
          jac_loc(3,6)%m(jr,jc) = jac_loc(3,6)%m(jr,jc) &
          - dt_fac*basis_vals(jr)*(btmp(1)*tmp2(2)+btmp(2)*tmp2(1))*int_factor/(m_i*n*mu0)
        ELSE
          tmp2 = cross_product(basis_grads(:,jc), [0.d0,1.d0,0.d0])
          DO l=1,3
            jac_loc(l+1,6)%m(jr,jc) = jac_loc(l+1,6)%m(jr,jc) &
            + dt_fac*DOT_PRODUCT(basis_grads(:,jr), tmp2)*btmp(l)*int_factor/(m_i*n*mu0) &
            + dt_fac*DOT_PRODUCT(basis_grads(:,jr), btmp)*tmp2(l)*int_factor/(m_i*n*mu0) &
            - dt_fac*basis_grads(l,jr)*DOT_PRODUCT(tmp2, btmp)*int_factor/(m_i*n*mu0) &
            - dt_fac*basis_vals(jr)*DOT_PRODUCT(dn, tmp2)*btmp(l)*int_factor/(m_i*n**2*mu0) &
            - dt_fac*basis_vals(jr)*DOT_PRODUCT(dn, btmp)*tmp2(l)*int_factor/(m_i*n**2*mu0) &
            + dt_fac*basis_vals(jr)*dn(l)*DOT_PRODUCT(tmp2, btmp)*int_factor/(m_i*n**2*mu0) 
          END DO
        END IF
        ! --vel, by
        IF (cyl_flag) THEN
          tmp2 = [0.d0,basis_vals(jc)/(coords(1)+gs_epsilon),0.d0]! this is 'B_y'
          DO l=1,3
            jac_loc(l+1,7)%m(jr,jc) = jac_loc(l+1,7)%m(jr,jc) &
            + dt_fac*DOT_PRODUCT(basis_grads(:,jr),tmp2)*btmp(l)*int_factor*coords(1)/(m_i*n*mu0) &
            + dt_fac*DOT_PRODUCT(basis_grads(:,jr), btmp)*tmp2(l)*int_factor*coords(1)/(m_i*n*mu0) &
            - dt_fac*basis_grads(l,jr)*DOT_PRODUCT(tmp2, btmp)*int_factor*coords(1)/(m_i*n*mu0) &
            - dt_fac*basis_vals(jr)*DOT_PRODUCT(dn, tmp2)*btmp(l)*int_factor*coords(1)/(m_i*n**2*mu0) &
            - dt_fac*basis_vals(jr)*DOT_PRODUCT(dn, btmp)*tmp2(l)*int_factor*coords(1)/(m_i*n**2*mu0) &
            + dt_fac*basis_vals(jr)*dn(l)*DOT_PRODUCT(tmp2, btmp)*int_factor*coords(1)/(m_i*n**2*mu0)
          END DO
          jac_loc(2,7)%m(jr,jc) = jac_loc(2,7)%m(jr,jc) &
          + dt_fac*basis_vals(jr)*(btmp(2)*tmp2(2)-btmp(1)*tmp2(1)-btmp(3)*tmp2(3))*int_factor/(m_i*n*mu0)
          jac_loc(3,7)%m(jr,jc) = jac_loc(3,7)%m(jr,jc) &
          - dt_fac*basis_vals(jr)*(btmp(1)*tmp2(2)+btmp(2)*tmp2(1))*int_factor/(m_i*n*mu0)
        ELSE
          tmp2 = [0.d0,basis_vals(jc),0.d0] 
          DO l=1,3
            jac_loc(l+1,7)%m(jr,jc) = jac_loc(l+1,7)%m(jr,jc) &
            + dt_fac*DOT_PRODUCT(basis_grads(:,jr),tmp2)*btmp(l)*int_factor/(m_i*n*mu0) &
            + dt_fac*DOT_PRODUCT(basis_grads(:,jr), btmp)*tmp2(l)*int_factor/(m_i*n*mu0) &
            - dt_fac*basis_grads(l,jr)*DOT_PRODUCT(tmp2, btmp)*int_factor/(m_i*n*mu0) &
            - dt_fac*basis_vals(jr)*DOT_PRODUCT(dn, tmp2)*btmp(l)*int_factor/(m_i*n**2*mu0) &
            - dt_fac*basis_vals(jr)*DOT_PRODUCT(dn, btmp)*tmp2(l)*int_factor/(m_i*n**2*mu0) &
            + dt_fac*basis_vals(jr)*dn(l)*DOT_PRODUCT(tmp2, btmp)*int_factor/(m_i*n**2*mu0)
          END DO
        END IF
        ! Temperature
        ! T, n
        IF (cyl_flag) THEN
          jac_loc(5, 1)%m(jr,jc) = jac_loc(5, 1)%m(jr, jc) &  
            - dt_fac*chi*basis_vals(jr)*DOT_PRODUCT(basis_grads(:, jc), dT)*int_factor*coords(1)/n & ! grad(delta_n) 
            + dt_fac*chi*basis_vals(jr)*basis_vals(jc)*DOT_PRODUCT(dn, dT)*int_factor*coords(1)/(n**2) ! delta_n 
        ELSE
          IF (linear) THEN
            jac_loc(5, 1)%m(jr,jc) = jac_loc(5, 1)%m(jr, jc) &  
            + dt_fac*basis_vals(jr)*basis_vals(jc)*DOT_PRODUCT(vel, dT)*int_factor/(n*(gamma-1))&
            + dt_fac*basis_vals(jr)*basis_vals(jc)*k_boltz*T*div_vel*int_factor/n &
            - dt_fac*basis_vals(jc)*chi*DOT_PRODUCT(basis_grads(:, jr), dT)*int_factor/n &
            - dt_fac*basis_vals(jr)*chi*basis_vals(jc)*DOT_PRODUCT(dn,dT)*int_factor/(n**2) 
          ELSE
            jac_loc(5, 1)%m(jr,jc) = jac_loc(5, 1)%m(jr, jc) &  
            - dt_fac*chi*basis_vals(jr)*DOT_PRODUCT(basis_grads(:, jc), dT)*int_factor/n & ! grad(delta_n) 
            + dt_fac*chi*basis_vals(jr)*basis_vals(jc)*DOT_PRODUCT(dn, dT)*int_factor/(n**2) ! delta_n 
          END IF
        END IF
        ! T, vel
        IF (cyl_flag) THEN
          DO l=1,3
            jac_loc(5, l+1)%m(jr,jc) = jac_loc(5, l+1)%m(jr, jc) &  
            + basis_vals(jr) * dt_fac*basis_vals(jc)*dT(l)*int_factor*coords(1)/(gamma-1) & ! delta_u dot Delta_T
            + basis_vals(jr) * dt_fac*T*basis_grads(l, jc)*int_factor*coords(1) ! div(delta u) = SUM(basis_grads)?
          END DO
          jac_loc(5, 2)%m(jr,jc) = jac_loc(5, 2)%m(jr, jc) &
            + basis_vals(jr)*dt_fac*T*basis_vals(jc)*int_factor 
        ELSE
          DO l=1,3
            jac_loc(5, l+1)%m(jr,jc) = jac_loc(5, l+1)%m(jr, jc) &  
            + basis_vals(jr) * dt_fac*basis_vals(jc)*dT(l)*int_factor/(gamma-1) & ! delta_u dot Delta_T
            + basis_vals(jr) * dt_fac*T*basis_grads(l, jc)*int_factor ! div(delta u) = SUM(basis_grads)?
          END DO
        END IF
        ! T, T
        IF(cyl_flag) THEN
          jac_loc(5, 5)%m(jr,jc) = jac_loc(5,5)%m(jr, jc) &  
          + basis_vals(jr) * basis_vals(jc)*int_factor*coords(1)/(gamma-1) &! delta_T
          + basis_vals(jr) * dt_fac*DOT_PRODUCT(vel, basis_grads(:, jc))*int_factor*coords(1)/(gamma-1) & ! nabla(dT)
          + basis_vals(jr) * dt_fac*basis_vals(jc)*div_vel*int_factor*coords(1) & ! dT != nabla(dT)
          + dt_fac * chi * DOT_PRODUCT(basis_grads(:, jc),basis_grads(:, jr))*int_factor*coords(1) & ! dT_Chi
          - dt_fac*basis_vals(jr)*chi*DOT_PRODUCT(dn, basis_grads(:, jc))*int_factor*coords(1)/n 
        ELSE
          jac_loc(5, 5)%m(jr,jc) = jac_loc(5,5)%m(jr, jc) &  
          + basis_vals(jr) * basis_vals(jc)*int_factor/(gamma-1) &! delta_T
          + basis_vals(jr) * dt_fac*DOT_PRODUCT(vel, basis_grads(:, jc))*int_factor/(gamma-1) & ! nabla(dT)
          + basis_vals(jr) * dt_fac*basis_vals(jc)*div_vel*int_factor & ! dT != nabla(dT)
          + dt_fac * chi * DOT_PRODUCT(basis_grads(:, jc),basis_grads(:, jr))*int_factor & ! dT_Chi
          - dt_fac*basis_vals(jr)*chi*DOT_PRODUCT(dn, basis_grads(:, jc))*int_factor/n 
        END IF
        ! Induction
        !--psi, vel
        IF (cyl_flag) THEN
          DO l=1,3
            jac_loc(6,l+1)%m(jr,jc) = jac_loc(6,l+1)%m(jr,jc) &
            + basis_vals(jr)*dt_fac*basis_vals(jc)*dpsi(l)*int_factor/(coords(1)+gs_epsilon)
            tmp2 = [0.d0, 0.d0,0.d0]
            tmp2(l) = 1.d0
            tmp3 = cross_product(B_0, tmp2)
            jac_loc(6,l+1)%m(jr,jc) = jac_loc(6,l+1)%m(jr,jc) &
            + basis_vals(jr)*dt_fac*basis_vals(jc)*tmp3(2)*int_factor/(coords(1)+gs_epsilon)
          END DO
        ELSE
          DO l=1,3
            jac_loc(6,l+1)%m(jr,jc) = jac_loc(6,l+1)%m(jr,jc) &
            + basis_vals(jr)*dt_fac*basis_vals(jc)*dpsi(l)*int_factor
            tmp2 = [0.d0, 0.d0,0.d0]
            tmp2(l) = 1.d0
            tmp3 = cross_product(B_0, tmp2)
            jac_loc(6,l+1)%m(jr,jc) = jac_loc(6,l+1)%m(jr,jc) &
            + basis_vals(jr)*dt_fac*basis_vals(jc)*tmp3(2)*int_factor
          END DO
        END IF
        ! --psi, psi
        IF (cyl_flag) THEN
          jac_loc(6, 6)%m(jr,jc) = jac_loc(6, 6)%m(jr,jc) &
          + basis_vals(jr)*basis_vals(jc)*int_factor/(coords(1)+gs_epsilon) &
          + dt_fac*basis_vals(jr)*DOT_PRODUCT(vel,basis_grads(:,jc))*int_factor/(coords(1)+gs_epsilon) &
          + dt_fac*eta*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*int_factor/(coords(1)+gs_epsilon)
        ELSE
          jac_loc(6, 6)%m(jr,jc) = jac_loc(6, 6)%m(jr,jc) &
          + basis_vals(jr)*basis_vals(jc)*int_factor &
          + dt_fac*basis_vals(jr)*DOT_PRODUCT(vel,basis_grads(:,jc))*int_factor &
          + dt_fac*eta*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*int_factor
        END IF
        !--by, vel
        tmp2 = cross_product(dpsi,basis_grads(:,jc))
        IF(cyl_flag) THEN
          DO l=1,3
            jac_loc(7,l+1)%m(jr,jc) = jac_loc(7, l+1)%m(jr,jc) &
            + basis_vals(jr)*dt_fac*basis_vals(jc)*dby(l)*int_factor/(coords(1)+gs_epsilon) &
            + basis_vals(jr)*dt_fac*by*basis_grads(l,jc)*int_factor/(coords(1)+gs_epsilon)
          END DO
          jac_loc(7,2)%m(jr,jc) = jac_loc(7, 2)%m(jr,jc) &
          + basis_vals(jr)*dt_fac*by*basis_vals(jc)*int_factor/(coords(1)+gs_epsilon)**2
          jac_loc(7,3)%m(jr,jc) = jac_loc(7, 3)%m(jr,jc) &
          -basis_vals(jr)*dt_fac*tmp2(2)*int_factor/(coords(1)+gs_epsilon)
        ELSE
          ! Full in-plane field (psi + background) in the stretch coupling, with the
          ! corrected sign — see nlfun_apply By block.
          tmp2 = cross_product(dpsi + cross_product([0.d0,1.d0,0.d0],B_0), basis_grads(:,jc))
          DO l=1,3
            jac_loc(7,l+1)%m(jr,jc) = jac_loc(7, l+1)%m(jr,jc) &
            + basis_vals(jr)*dt_fac*basis_vals(jc)*dby(l)*int_factor &
            + basis_vals(jr)*dt_fac*by*basis_grads(l,jc)*int_factor
            IF (l==2) THEN
              jac_loc(7,l+1)%m(jr,jc) = jac_loc(7, l+1)%m(jr,jc) &
              + basis_vals(jr)*dt_fac*tmp2(l)*int_factor
            END IF
          END DO
        END IF
        !-- by, psi
        tmp2 = cross_product(basis_grads(:,jc),dvel(2,:))
        IF (cyl_flag) THEN
          jac_loc(7, 6)%m(jr,jc) = jac_loc(7, 6)%m(jr,jc) &
          - basis_vals(jr)*dt_fac*tmp2(2)*int_factor/(coords(1)+gs_epsilon)
        ELSE
          ! Sign corrected with the stretch term (see nlfun_apply By block)
          jac_loc(7, 6)%m(jr,jc) = jac_loc(7, 6)%m(jr,jc) &
          + basis_vals(jr)*dt_fac*tmp2(2)*int_factor
        END IF
        !-- by, by
        IF (cyl_flag) THEN
          jac_loc(7, 7)%m(jr,jc) = jac_loc(7, 7)%m(jr,jc) &
          + basis_vals(jr)*basis_vals(jc)*int_factor/(coords(1)+gs_epsilon) &
          + dt_fac*basis_vals(jr)*dvel(1,1)*basis_vals(jc)*int_factor/(coords(1)+gs_epsilon) &
          + dt_fac*basis_vals(jr)*dvel(3,3)*basis_vals(jc)*int_factor/(coords(1)+gs_epsilon) &
          - dt_fac*basis_vals(jr)*vel(1)*basis_vals(jc)*int_factor/(coords(1)+gs_epsilon)**2 &
          + dt_fac*basis_vals(jr)*DOT_PRODUCT(vel, basis_grads(:,jc))*int_factor/(coords(1)+gs_epsilon) &
          + dt_fac*eta*DOT_PRODUCT(basis_grads(:,jr), basis_grads(:,jc))*int_factor/(coords(1)+gs_epsilon)
        ELSE
          jac_loc(7, 7)%m(jr,jc) = jac_loc(7, 7)%m(jr,jc) &
          + basis_vals(jr)*basis_vals(jc)*int_factor &
          + basis_vals(jr)*dt_fac*DOT_PRODUCT(basis_grads(:,jc),vel)*int_factor & 
          + basis_vals(jr)*dt_fac*basis_vals(jc)*div_vel*int_factor & 
          + dt_fac*eta*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*int_factor
        END IF
      END DO
    END DO
  END DO
  DO jr = 1,7
    jac_loc(1, jr)%m = jac_loc(1, jr)%m / self%den_scale
    jac_loc(jr, 1)%m = jac_loc(jr, 1)%m * self%den_scale
  END DO
  !---Phase-4 thin-wall Robin boundary term on by (see nlfun_apply): per-cell
  !  edge-mass contribution +dt_fac*(eta/c_wall)*|e|/6*[[2,1],[1,2]] on (7,7)
  !  for boundary edges of this cell with both endpoints free. Each boundary
  !  edge belongs to exactly one cell, so this adds each edge once.
  IF(self%c_wall > 0.d0)THEN
    BLOCK
    INTEGER(i4) :: ke, ed, ip1, ip2, jr1, jr2, jq
    REAL(r8) :: elen, rfac
    DO ke=1,3
      ed = ABS(mesh%lce(ke,i))
      IF(.NOT.mesh%be(ed))CYCLE
      ip1 = mesh%le(1,ed); ip2 = mesh%le(2,ed)
      IF(self%by_bc(ip1).OR.self%by_bc(ip2))CYCLE
      jr1 = 0; jr2 = 0
      DO jq=1,oft_blagrange%nce
        IF(cell_dofs(jq)==ip1)jr1 = jq
        IF(cell_dofs(jq)==ip2)jr2 = jq
      END DO
      IF(jr1==0.OR.jr2==0)CYCLE
      elen = SQRT(SUM((mesh%r(:,ip1)-mesh%r(:,ip2))**2))
      rfac = dt_fac*eta*elen/(self%c_wall*6.d0)
      jac_loc(7,7)%m(jr1,jr1) = jac_loc(7,7)%m(jr1,jr1) + 2.d0*rfac
      jac_loc(7,7)%m(jr1,jr2) = jac_loc(7,7)%m(jr1,jr2) + rfac
      jac_loc(7,7)%m(jr2,jr1) = jac_loc(7,7)%m(jr2,jr1) + rfac
      jac_loc(7,7)%m(jr2,jr2) = jac_loc(7,7)%m(jr2,jr2) + 2.d0*rfac
    END DO
    END BLOCK
  END IF
  !---Apply Boundary Conditions
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%n_bc(cell_dofs),1)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%velx_bc(cell_dofs), 2)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%vely_bc(cell_dofs), 3)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%velz_bc(cell_dofs), 4)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%T_bc(cell_dofs),5)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%psi_bc(cell_dofs),6)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%by_bc(cell_dofs), 7)
  !----Add local contributions to matrix
  !$omp ordered
  CALL self%fe_rep%mat_add_local(self%jacobian,jac_loc,iloc,tlocks)
  !$omp end ordered
END DO
!---Cleanup thread-local storage
CALL self%fe_rep%mat_destroy_local(jac_loc)
DEALLOCATE(basis_vals,basis_grads,T_weights_loc,vel_weights_loc, &
          n_weights_loc, psi_weights_loc, by_weights_loc, cell_dofs,jac_loc,iloc)
!$omp end parallel
END BLOCK
!--Destroy thread locks
DO i=1,self%fe_rep%nfields
  CALL omp_destroy_lock(tlocks(i))
END DO
DEALLOCATE(tlocks)
IF(oft_debug_print(2))write(*,'(4X,A)')'Setting BCs'
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%n_bc,1)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%velx_bc,2)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%vely_bc,3)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%velz_bc,4)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%T_bc,5)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%psi_bc,6)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%by_bc,7)
!
call self%fe_rep%vec_create(tmp)
call self%jacobian%assemble(tmp)
call tmp%delete
DEALLOCATE(tmp,n_weights,vel_weights, T_weights, &
          by_weights, psi_weights)
end subroutine build_approx_jacobian

!---------------------------------------------------------------------------
!> Update matrix-free Jacobian on all levels with new solution
!---------------------------------------------------------------------------
subroutine mfnk_update(uin)
class(oft_vector), target, intent(inout) :: uin !< Current field
IF(oft_debug_print(1))write(*,*)'Updating 2D MUG MF-Jacobian'
CALL current_sim%mf_mat%update(uin)
END SUBROUTINE mfnk_update
!---------------------------------------------------------------------------
!> Update Jacobian matrices on all levels with new solution
!---------------------------------------------------------------------------
subroutine update_jacobian(uin)
class(oft_vector), target, intent(inout) :: uin !< Current solution
IF(oft_debug_print(1))write(*,*)'Updating 2D MUG approximate Jacobian'
CALL build_approx_jacobian(current_sim,uin)
END SUBROUTINE update_jacobian
!---------------------------------------------------------------------------
!> Setup composite FE representation and ML environment
!---------------------------------------------------------------------------
subroutine setup(self,mg_mesh_in, order)
class(oft_xmhd_2d_sim), intent(inout) :: self
CLASS(multigrid_mesh), TARGET, intent(in) :: mg_mesh_in
integer(i4), intent(in) :: order
integer(i4) :: i,j, ierr,io_unit, cond_ind, coil_ind, type
LOGICAL, ALLOCATABLE :: vert_flag(:),edge_flag(:), boundary_flag(:)
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs
mg_mesh=>mg_mesh_in
mesh=>mg_mesh%smesh
IF(ASSOCIATED(self%fe_rep))CALL oft_abort("Setup can only be called once","setup",__FILE__)
IF(ASSOCIATED(oft_blagrange))CALL oft_abort("FE space already built","setup",__FILE__)

!---Look for XML defintion elements
#ifdef HAVE_XML
IF(ASSOCIATED(oft_env%xml))THEN
  CALL xml_get_element(oft_env%xml,"xmhd2d",self%xml_root,ierr)
  IF(ierr==0)THEN
    !---Look for pre node
    CALL xml_get_element(self%xml_root,"pre",self%xml_pre_def,ierr)
    IF(ierr/=0)NULLIFY(self%xml_pre_def)
  ELSE
    NULLIFY(self%xml_root)
  END IF
END IF
#endif

!---Setup FE representation
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Building lagrange FE space'
CALL oft_lag_setup(mg_mesh,order,ML_blag_obj=ML_oft_blagrange,minlev=-1)
IF(.NOT.oft_2D_lagrange_cast(oft_blagrange,ML_oft_blagrange%current_level))CALL oft_abort("Invalid lagrange FE object","setup",__FILE__)

!---Build composite FE definition for solution field
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Creating FE type'
ALLOCATE(self%fe_rep)
self%fe_rep%nfields=7
ALLOCATE(self%fe_rep%fields(self%fe_rep%nfields))
ALLOCATE(self%fe_rep%field_tags(self%fe_rep%nfields))
self%fe_rep%fields(1)%fe=>oft_blagrange
self%fe_rep%field_tags(1)='n'
self%fe_rep%fields(2)%fe=>oft_blagrange
self%fe_rep%field_tags(2)='velx'
self%fe_rep%fields(3)%fe=>oft_blagrange
self%fe_rep%field_tags(3)='vely'
self%fe_rep%fields(4)%fe=>oft_blagrange
self%fe_rep%field_tags(4)='velz'
self%fe_rep%fields(5)%fe=>oft_blagrange
self%fe_rep%field_tags(5)='T'
self%fe_rep%fields(6)%fe=>oft_blagrange
self%fe_rep%field_tags(6)='psi'
self%fe_rep%fields(7)%fe=>oft_blagrange
self%fe_rep%field_tags(7)='by'

!---Create solution vector
CALL self%fe_rep%vec_create(self%u)
CALL self%fe_rep%vec_create(self%u0)

!---Create Jacobian matrix
ALLOCATE(self%jacobian_block_mask(self%fe_rep%nfields,self%fe_rep%nfields))
self%jacobian_block_mask=1
CALL self%fe_rep%mat_create(self%jacobian,self%jacobian_block_mask)
end subroutine setup

!---------------------------------------------------------------------------
!> Set any boundary conditions not already set with default values
!---------------------------------------------------------------------------
subroutine setup_bc(self)
class(oft_xmhd_2d_sim), intent(inout) :: self
IF(self%c_wall > 0.d0)THEN
  IF(oft_blagrange%order > 1)CALL oft_abort( &
    'Thin-wall Robin BC (c_wall > 0) requires order = 1', 'setup_bc', __FILE__)
  IF(self%cyl_flag)CALL oft_abort( &
    'Thin-wall Robin BC (c_wall > 0) not implemented for cylindrical mode', 'setup_bc', __FILE__)
END IF
IF(.NOT.ASSOCIATED(self%n_bc))self%n_bc=>oft_blagrange%global%gbe
IF(.NOT.ASSOCIATED(self%velx_bc))self%velx_bc=>oft_blagrange%global%gbe
IF(.NOT.ASSOCIATED(self%vely_bc))self%vely_bc=>oft_blagrange%global%gbe
IF(.NOT.ASSOCIATED(self%velz_bc))self%velz_bc=>oft_blagrange%global%gbe
IF(.NOT.ASSOCIATED(self%T_bc))self%T_bc=>oft_blagrange%global%gbe
IF(.NOT.ASSOCIATED(self%psi_bc))self%psi_bc=>oft_blagrange%global%gbe
IF(.NOT.ASSOCIATED(self%by_bc))self%by_bc=>oft_blagrange%global%gbe
end subroutine setup_bc

!---------------------------------------------------------------------------
!> Save xMHD solution state to a restart file
!---------------------------------------------------------------------------
subroutine rst_save(self,u,t,dt,filename,path)
class(oft_xmhd_2d_sim), intent(inout) :: self
class(oft_vector), pointer, intent(inout) :: u !< Solution to save
real(r8), intent(in) :: t !< Current solution time
real(r8), intent(in) :: dt !< Current timestep
character(LEN=*), intent(in) :: filename !< Name of restart file
character(LEN=*), intent(in) :: path !< Path to store solution vector in file
DEBUG_STACK_PUSH
CALL self%fe_rep%vec_save(u,filename,path)
IF(oft_env%head_proc)THEN
  CALL hdf5_write(t,filename,'t')
  CALL hdf5_write(dt,filename,'dt')
END IF
DEBUG_STACK_POP
end subroutine rst_save

!---------------------------------------------------------------------------
!> Plot the current solution state using XDMF/HDF5 output
!! Runtime options are set in the main input file using the group `xmhd_plot_options`.
!! **Option group:** `xmhd_plot_options`
!! |  Option    |  Description  | Type [dim] |
!! |------------|---------------|------------|
!! | `rst_start=0`      | First restart file to read | int |
!! | `rst_end=2000`     | Last restart file to read | int |
!------------------------------------------------------------------------------
!---------------------------------------------------------------------------
subroutine xmhd_2d_plot(self)
class(oft_xmhd_2d_sim), intent(inout) :: self
class(oft_vector), pointer :: ux,uy,uz,v_lag
type(oft_lag_bginterp) :: grad_psi
class(oft_matrix), pointer :: lmop => NULL()
real(r8), pointer :: plot_vals(:),plot_vec(:,:),plot_u0(:)
!---Solver objects
CLASS(oft_solver), POINTER :: lminv => NULL()
class(oft_vector), pointer :: u,v,up, u0
class(oft_matrix), pointer :: lap_mat => NULL()
INTEGER(i4):: rst_start=0
INTEGER(i4):: rst_end=2000
INTEGER(i4) :: rst_cur, rst_tmp, ierr, io_stat, io_unit
CHARACTER(LEN=OFT_PATH_SLEN) :: file_tmp
LOGICAL :: rst_exist
real(r8) :: t
TYPE(xdmf_plot_file) :: xdmf_plot
character(LEN=XMHD_RST_LEN) :: rst_char
namelist/xmhd_plot_options/rst_start,rst_end
!---------------------------------------------------------------------------
! Read plotting options
!---------------------------------------------------------------------------
open(NEWUNIT=io_unit,FILE=oft_env%ifile)
read(io_unit,xmhd_plot_options,IOSTAT=ierr)
close(io_unit)
!---------------------------------------------------------------------------
! Create solver fields
!---------------------------------------------------------------------------
call self%fe_rep%vec_create(u)
call self%fe_rep%vec_create(up)
call self%fe_rep%vec_create(v)
call self%fe_rep%vec_create(u0)
!---
call oft_blagrange%vec_create(grad_psi%u)
call oft_blagrange%vec_create(ux)
call oft_blagrange%vec_create(uy)
call oft_blagrange%vec_create(uz)
call oft_blagrange%vec_create(v_lag)
NULLIFY(lmop)
call oft_blag_getmop(oft_blagrange,lmop)
CALL create_cg_solver(lminv)
lminv%A=>lmop
lminv%its=-2
CALL create_diag_pre(lminv%pre)
ALLOCATE(plot_vec(3,v_lag%n))
NULLIFY(plot_vals,plot_u0)
CALL grad_psi%setup(oft_blagrange)

CALL xdmf_plot%setup("xmhd_2d")
CALL mesh%setup_io(xdmf_plot,oft_blagrange%order)

!---------------------------------------------------------------------------
! If linear, extract equilibrium fields
!---------------------------------------------------------------------------
IF (self%linear) THEN
  file_tmp='xmhd2d_eq.rst'
  rst_exist=oft_file_exist(TRIM(file_tmp))
  CALL oft_mpi_barrier(ierr)
  IF(.NOT.rst_exist)THEN
    CALL oft_abort('Equilibrium restart file does not exist.','xmhd_plot',__FILE__)
  ELSE
    CALL hdf5_read(t,file_tmp,'t')
    CALL self%rst_load(u0,file_tmp, 'U')
  END IF
END IF

!---------------------------------------------------------------------------
! Loop through rst files
!---------------------------------------------------------------------------
100 FORMAT (I XMHD_RST_LEN.XMHD_RST_LEN)
rst_cur = rst_start
DO
  IF(rst_cur > rst_end) EXIT
  ! Determine file name
  WRITE(rst_char,100)rst_cur
  READ(rst_char,100,IOSTAT=io_stat)rst_tmp
  IF((io_stat/=0).OR.(rst_tmp/=rst_cur))CALL oft_abort("Step count exceeds format width", "xmhd_plot", __FILE__)
  file_tmp='xmhd2d_'//rst_char//'.rst'
  !Read file
  rst_exist=oft_file_exist(TRIM(file_tmp))
  CALL oft_mpi_barrier(ierr)
  IF(.NOT.rst_exist)THEN
    EXIT !CALL oft_abort('Restart file does not exist.','xmhd_plot',__FILE__)
  ELSE
    CALL hdf5_read(t,file_tmp,'t')
    CALL self%rst_load(u,file_tmp, 'U')
  END IF
  IF (self%linear) THEN
    CALL u%add(1.d0, 1.d0, u0)
  END IF
  !Plot data
  CALL xdmf_plot%add_timestep(t)

  !Plot density
  CALL u%get_local(plot_vals,1)
  plot_vals = plot_vals*self%den_scale
  CALL mesh%save_vertex_scalar(plot_vals,xdmf_plot,'n')
  !Plot velocity
  CALL u%get_local(plot_vals,2)
  plot_vec(1,:)=plot_vals
  CALL u%get_local(plot_vals,3)
  plot_vec(3,:)=plot_vals
  CALL u%get_local(plot_vals,4)
  plot_vec(2,:)=plot_vals
  CALL mesh%save_vertex_vector(plot_vec,xdmf_plot,'V')
  !Plot temperature
  CALL u%get_local(plot_vals,5)
  CALL mesh%save_vertex_scalar(plot_vals,xdmf_plot,'T')
  !Plot psi
  CALL u%get_local(plot_vals,6)
  CALL mesh%save_vertex_scalar(plot_vals,xdmf_plot,'psi')
  !Plot B
  CALL grad_psi%u%restore_local(plot_vals)
  CALL grad_psi%setup(oft_blagrange)
  CALL oft_blag_vproject(oft_blagrange,grad_psi,ux,uy,uz)
  CALL v_lag%set(0.d0)
  CALL lminv%apply(v_lag,ux)
  CALL ux%add(0.d0,1.d0,v_lag)
  CALL v_lag%set(0.d0)
  CALL lminv%apply(v_lag,uy)
  CALL uy%add(0.d0,1.d0,v_lag)
  CALL uy%get_local(plot_vals)
  plot_vec(1,:)=-plot_vals
  CALL ux%get_local(plot_vals)
  plot_vec(2,:)=plot_vals
  CALL u%get_local(plot_vals,7)
  plot_vec(3,:)=plot_vals
  IF (self%cyl_flag) THEN
    CALL mesh%save_vertex_vector(plot_vec,xdmf_plot,'B*R')
  ELSE
    CALL mesh%save_vertex_vector(plot_vec,xdmf_plot,'B')
  END IF

  !Move to next file
  rst_cur = rst_cur + self%rst_freq
END DO
CALL u%delete
CALL up%delete
CALL v%delete
CALL u0%delete
CALL ux%delete
CALL uy%delete
CALL uz%delete
DEALLOCATE(u,up,v,u0,ux,uy,uz)
IF(ASSOCIATED(plot_vals))DEALLOCATE(plot_vals)
IF(ASSOCIATED(plot_vec))DEALLOCATE(plot_vec)
IF(ASSOCIATED(plot_u0))DEALLOCATE(plot_u0)
end subroutine xmhd_2d_plot
!---------------------------------------------------------------------------
!> Load xMHD solution state from a restart file
!---------------------------------------------------------------------------
subroutine rst_load(self,u,filename,path,t,dt)
class(oft_xmhd_2d_sim), intent(inout) :: self
class(oft_vector), pointer, intent(inout) :: u !< Solution to load
character(LEN=*), intent(in) :: filename !< Name of restart file
character(LEN=*), intent(in) :: path !< Path to store solution vector in file
real(r8), optional, intent(out) :: t !< Time of loaded solution
real(r8), optional, intent(out) :: dt !< Timestep at loaded time
DEBUG_STACK_PUSH
CALL self%fe_rep%vec_load(u,filename,path)
IF(PRESENT(t))CALL hdf5_read(t,filename,'t')
IF(PRESENT(dt))CALL hdf5_read(dt,filename,'dt')
DEBUG_STACK_POP
end subroutine rst_load

SUBROUTINE r_init(pt,val)
  REAL(r8), INTENT(in) :: pt(3)
  REAL(r8), INTENT(out) :: val
  val = pt(1)
END SUBROUTINE r_init
END MODULE xmhd_2d