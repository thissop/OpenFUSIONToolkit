# Item 3 build — segmented multi-channel + transient inter-channel leakage cross-talk

The one adversary-survivor (novel + un-anchored): transient dB/dt-driven leakage cross-talk
`flow_split_i(t)=Q_i(t)/<Q>` and septum leakage current `I_leak(t)` across N conducting-septa
channels. Steady maldistribution = VALIDATION cross-check (Smolentsev 2014, CFETR-COOL Nucl.Fusion
2024, Mistrangelo-Buhler); Tier-2 anchor = Molokov/JFM-2015 asymptotics; Tier-1 = insulating septa ->
N decoupled Shercliff/Hunt. Template for the whole region-material port = the 3D solver `xmhd.F90`,
which already has it.

## Stage A — region-tagged 2D mesh  [DONE, guarded]
`src/grid/mesh_cube.F90` `smesh_square_load`: added namelist `septum_lo, septum_hi`; cells with
center-x in [lo,hi] -> reg 2 (septum), x>hi -> reg 3, else reg 1. Default lo>hi => single region
(byte-identical). Tags both the tri (mesh_type=1) and quad paths. **Must compile-verify on Ginsburg.**

## Stage B — region-material state in xmhd_2d.F90  [TODO — mirror xmhd.F90 exactly]
xmhd.F90 template line refs: `solid_node`(182), `solid_cell`(326), `eta_reg`(332),
`eta_curr=eta*eta_reg(mesh%reg(i))`(1656), solid-cell branch(1679), solid_node zeroing(708/1262/2561).
Port into oft_xmhd_2d_sim:
1. TYPE (after a_half ~L139): `REAL(r8),CONTIGUOUS,POINTER :: eta_reg(:)=>NULL()` (per-region resistivity
   scale, default 1); `LOGICAL,CONTIGUOUS,POINTER :: solid_reg(:)=>NULL()`; `LOGICAL,ALLOCATABLE ::
   solid_cell(:)`, `solid_node(:)`. Unassociated => scalar-eta path, single-region runs bit-identical.
2. setup (in setup_bc or run_simulation): if eta_reg associated, build solid_cell from solid_reg via
   mesh%reg; build solid_node (a node is solid iff ALL incident cells solid) mirroring xmhd.F90:2561.
3. nlfun_apply quadrature loop: `eta_curr = eta`; if eta_reg associated `eta_curr=eta*eta_reg(mesh%reg(m))`.
   Replace scalar `eta` with `eta_curr` in the psi (L~1150) and by (L~1174-1181) diffusion terms ONLY
   (not the velocity eqs). In solid cells (solid_cell(m)): keep the by/psi diffusion (wall carries
   induced current) but DROP the momentum/advection/buoyancy (no flow in a solid).
4. Velocity Dirichlet in solid: force vel=0 at solid_node (mirror xmhd.F90:708/1262) — add to the
   fem_dirichlet_vec BC application after L1214-1216, or set velx/y/z_bc TRUE at solid_node.
5. build_approx_jacobian: use eta_curr in the psi/by diffusion Jacobian blocks; zero momentum rows at
   solid_node.
6. Guard: all gated on `ASSOCIATED(self%eta_reg)`. Regression gate = Shercliff/Hunt bit-identical.

## Stage C — driver  [TODO]
`test_multichannel_mug.F90` (from test_hunt_mug): namelist adds septum_lo/hi (-> cube_options),
eta_reg/solid_reg per region (region 2 = septum: solid_reg(2)=.TRUE., eta_reg(2)=eta_wall/eta_fluid),
c_wall on the outer Hartmann walls, and by_source disruption ramp (dB, tau_q). Two fluid channels
(reg 1,3) + septum (reg 2). Output per-channel flow rate Q_i(t) (integrate vely over reg==1 and
reg==3 separately) and septum leakage current I_leak(t) = oint J.dl around the septum.

## Stage D — validation  [TODO]
- Tier 1: insulating septum (eta_reg(2) huge) => two decoupled Shercliff/Hunt ducts; flow_split=1 each.
- Tier 2 (steady): finite-c septum, high Ha steady => coupled split vs Molokov/JFM-2015; vs Smolentsev
  2014 / CFETR-COOL COMSOL. This is the VALIDATION (reproduce known steady maldistribution).
- Novel (transient): disruption dB/dt ramp => flow_split_i(t), I_leak(t) waveforms. COMSOL multi-domain
  conjugate run (one steady + one time-dependent) = the independent cross-check.

## Build/verify loop
Sync changed files to Ginsburg repo `/burg/astro/users/tjk2147/spf_work/oft_lm_mhd_papers/repo`;
`cd builds/build_release && make oftgrid` (Stage A) then `make test_hunt_mug test_multichannel_mug`
(Stages B/C). Runtime env: LD_LIBRARY_PATH=/cm/local/apps/gcc/11.2.0/lib64:/burg/opt/anaconda3-2023.09/lib.
Regression gate before trusting any multichannel result: a single-region Shercliff run must be
bit-identical to the pre-change binary.
