# Liquid-Metal MHD in OFT/MUG: Current State, Independent Verification, and Upgrade Roadmap

Prepared by Thomas Kiker (Columbia University) for Yuchen Jiang (ORNL).
Branch `lm-mhd-upgrades` of OpenFUSIONToolkit. Compiled 2026-06-30.

This document is a single self-contained account of the liquid-metal MHD (LM-MHD) work on
the OpenFUSIONToolkit (OFT) two-dimensional MHD solver known as MUG. It folds together the
several scattered working notes in the repository into one narrative, and it states plainly
what has been verified against artifacts on disk and what remains a design proposal. The
framing throughout is a personal verification and engineering-readiness exercise. The work
reproduces classical Hartmann, Shercliff, and Hunt duct-flow results as an independent
cross-check and scopes a set of upgrades; it does not validate MUG or claim new physics
capability, and the reduced wall-drag closure was disabled in every run reported here.

## 1. What the MUG solver (`src/physics/xmhd_2d.F90`) solves today

The module `src/physics/xmhd_2d.F90` is a self-contained two-dimensional, two-and-a-half
dimensional in velocity, compressible visco-resistive MHD solver. It is built on a scalar
Lagrange finite-element basis over a two-dimensional surface mesh and is independent of the
three-dimensional `xmhd.F90` and `xmhd_lag.F90` modules, which use different mesh types,
finite-element spaces, and solvers. It advances seven coupled scalar Lagrange fields: density
`n`, the three velocity components `velx`, `vely`, `velz`, temperature `T`, the in-plane
poloidal flux `psi`, and the out-of-plane field component `by`. Carrying three velocity
components on a two-dimensional mesh is what makes the model two-and-a-half dimensional. The
magnetic field is reconstructed as the in-plane part from the flux, the out-of-plane part from
`by`, and a uniform applied background `B_0`. Both a Cartesian path and a cylindrical R-Z path
(`cyl_flag`) are present, although the cylindrical linear solve aborts at present and only the
nonlinear advance supports R-Z.

The governing system, in Cartesian form, is continuity with a numerical density diffusivity,
a momentum equation carrying the Lorentz force `J x B`, the equation-of-state pressure gradient,
advection, and viscous diffusion, an internal-energy equation for temperature, and the full
induction equations for `psi` and `by`. Pressure is the thermodynamic quantity `p = 2 n k_B T`
(electron plus ion, hence the factor of two) computed from the density and temperature fields;
there is no independent pressure degree of freedom and no incompressibility constraint. The
advective EMF `u x B` is present and is nonlinear in the evolved flux and field fields; the
background-field contribution `B_0 x u` is linear in velocity, so the induction operator is not
uniformly nonlinear and should not be described as such. The transport coefficients `eta`
(resistivity in a magnetic-diffusivity normalization that includes `mu0`), `nu` (kinematic
viscosity), `chi` (thermal diffusivity), and `D_diff` (density diffusivity) are uniform scalars.
There is no anisotropic or tensor transport, no temperature-dependent or Spitzer resistivity, no
Hall term, and no Braginskii closure in this module. Velocity boundary conditions are strictly
Dirichlet, set per boundary node as either homogeneous (no-slip) or natural (zero traction).
There is no periodic, outflow, traction, or partial-slip velocity condition in the module itself.

Two points bear directly on how the Hartmann verification below was constructed. First, the
solver is compressible. It has no incompressible or low-Mach mode. Where the verification needed
incompressible behaviour, that was imposed at the case level by holding `n` and `T` fixed through
Dirichlet boundary conditions on those fields, which removes the pressure-gradient response and
the density evolution, leaving a constant-density momentum balance. This is a property of the
driving case, not of the solver. Second, because the velocity boundary conditions are
Dirichlet-only and there is no periodic option in the module, a fully-developed channel was
realized by meshing a periodic-in-streamwise domain at the mesh construction level and applying
no-slip only on the wall-normal boundaries.

## 2. The Phase-1 reduced wall-drag closure (already-validated background)

A separate and earlier piece of work added an optional reduced wall-drag, or Hartmann-friction,
source term to the same module. It is gated behind a new logical `use_wall_drag`, which defaults
to false, with companion parameters `drag_coeff` (the damping rate alpha, in inverse time) and
`drag_bhat` (a unit vector along the applied-field direction). The term is a volumetric momentum
sink `S_u = -alpha * u_perp`, where `u_perp = u - (u . bhat) bhat` is the velocity component
perpendicular to the applied field, so that the flow along the field is not damped. It is
integrated implicitly by backward Euler, which for a uniform field reduces per degree of freedom
to `v_perp(n+1) = v_perp(n) / (1 + alpha dt)`, the discrete form of exponential decay at rate
alpha. The term is added to both the nonlinear residual and the analytic Jacobian, in the
Cartesian and cylindrical paths, the Jacobian contribution carrying the perpendicular projector
`delta_kl - bhat_k bhat_l`.

This closure is a reduced model. It represents the integrated braking effect of the unresolved
Hartmann boundary layers as a bulk friction; it does not resolve the layer and cannot reproduce
the detailed velocity profile pointwise. In the broader planning notes it is identified with the
Sommeria-Moreau quasi-two-dimensional model, whose Hartmann-wall and side-wall friction
coefficients scale as `alpha_H ~ (sigma B^2 / rho)(1 / Ha)` and `alpha_S ~ alpha_H Ha^(-1/2)`.
Two cautions belong with any future use of those expressions. The dimensional normalization of
`drag_coeff` is not pinned down in the code or its comments, and the planning notes contain two
mutually inconsistent expressions for the Hartmann coefficient that differ by a factor of `Ha`.
Calibrating `drag_coeff` against a resolved reference is therefore a prerequisite to any
quantitative use of the closure.

What was validated about the closure is narrow and honest. The regression test
`test_drag_decay` confirms three things. The perpendicular velocity component decays at exactly
the backward-Euler rate, with a measured relative error of about five times ten to the minus five
across finite-element orders two and three and damping rates two and five (5.17e-5 in the recorded
`drag_decay.results`), set by the convergence of the iterative nonlinear solve rather than by the
drag model and far below the two percent pass threshold. The parallel velocity component is
undamped to a relative error of order ten to the minus thirteen (5.3e-13 recorded), which confirms
the perpendicular projection. And because the feature defaults to off, the existing Alfven-wave and
sound-wave tests are unchanged, which confirms the new path is inert when disabled. What was not
validated is everything spatial: the test uses a uniform field and a uniform initial condition,
so it exercises the decay rate and the projection, not any Hartmann-layer structure, not a
flow-rate reduction, and not the calibration of `drag_coeff` against a real duct flow. The
cylindrical drag path is implemented but is not yet exercised by any test.

The reduced closure is not the subject of this deliverable. It is described here only as
already-validated background, and it was held off in every verification run reported below. The
purpose of the verification was precisely to show that a resolved channel does not need it.

## 3. Independent three-way verification of the insulating Hartmann channel (the headline)

The central result is an independent, three-way cross-check of the classical insulating-wall
Hartmann channel. Three solutions of the same problem were computed and compared against one
another: the closed-form analytic profile as the referee, the resolved MUG solver with the drag
closure off, and the OpenFOAM solver `mhdFoam`. The artifacts behind every number quoted here
are in the accompanying `results/` directory, and all numbers in this section were recomputed
directly from the stored per-profile CSV files during this packaging pass rather than copied from
any earlier note.

The analytic referee is implemented in `analytic_hartmann.py` and was re-run during packaging.
Its self-test passes: the zero-Hartmann limit recovers the Poiseuille ratio `u_mean / u_c = 2/3`
and a parabolic profile to a relative error of order ten to the minus eight, the large-Hartmann
limit gives a plug core with `u_mean / u_c` approaching one, the dimensionless flow rate follows
the classical `Q* ~ 1/Ha` suppression, and the closed-form mean-over-centerline ratio agrees with
direct numerical integration of the profile to a relative error of order ten to the minus eight
across the sweep. The analytic module is therefore trustworthy as the referee.

The OpenFOAM leg used `mhdFoam`, the canonical incompressible MHD solver shipped with OpenFOAM,
on a wall-normal-graded mesh that resolves the Hartmann layer of thickness `a / Ha` with several
cells up to the highest Hartmann number. It was run for `Ha` in the set one, five, ten, twenty,
fifty. Velocity profiles were extracted by conversion to VTK followed by cell-centroid selection,
because the OpenFOAM sampling function object in the installed version fails with an IO error.
Recomputed from the stored profiles, `mhdFoam` agrees with the analytic profile to a normalized
L2 error of order one times ten to the minus three, the largest being 1.045 times ten to the minus
three at Hartmann number one, a normalized maximum error at or below two and a half times ten to
the minus three, and a relative error in the discriminating scalar `u_mean / u_c` of about one
times ten to the minus three across the whole sweep. The discriminating
scalar tracks the analytic value from the near-Poiseuille profile at `Ha = 1` to the plug core at
`Ha = 50`. This is the solid, guaranteed leg of the comparison.

The resolved MUG leg used the same module discussed in Section 1, with the drag closure off,
driven to a fully-developed steady state by a uniform body-force momentum source that was added
to the module as a second default-off option (`use_body_force`, `body_force(3)`). That source
enters the nonlinear residual only; being constant it contributes nothing to the Jacobian, and it
is inert when disabled. The driving case (`src/tests/physics/test_hartmann_mug.F90`) sets a
periodic-in-streamwise channel with no-slip walls, a wall-normal applied field, the insulating
natural boundary condition on the flux, and frozen density and temperature to enforce
incompressible behaviour. The streamwise velocity profile across the channel was extracted and
compared to analytic. The leg is confirmed at `Ha` equal to two, three, and five. At each of
those, the Hartmann number fitted from the computed profile shape matches the Hartmann number
predicted from the input parameters, through `Ha = B_0 a / sqrt(mu0 eta rho nu)`, to within
about one part in a thousand (fitted-over-predicted ratios of 1.0003, 1.0004, and 1.0007). The
normalized profile matches analytic to an L2 error of order one times ten to the minus four, and
the discriminating scalar matches to a relative error of order one times ten to the minus four.
This simultaneously confirms the geometry mapping and the parameter scaling. The flow was
verified to be fully developed, with the streamwise velocity uniform along the channel to about
one part in a million.

A note on two numbers that appeared inconsistent in the earlier write-ups, now reconciled. The
MUG L2 error at `Ha = 5` was quoted in one place as 1.25 times ten to the minus four and in
another as 1.53 times ten to the minus four. Both are correct. The smaller value is the error
against the analytic profile evaluated at the fitted Hartmann number; the larger value is the
error against the analytic profile evaluated at the integer Hartmann number. The rebuilt error
table reports both columns explicitly so the distinction is visible rather than confusing.

The headline three-way error table follows. Profiles are normalized by centerline velocity. The
L2 and maximum errors are on the normalized profile against the analytic referee, computed on each
code's own grid against analytic at the integer Hartmann number. The quantity `u_mean / u_c` is
the discriminating scalar, two-thirds at zero Hartmann number and approaching one in the plug
core. A blank entry means no run was produced for that source at that Hartmann number; no value
is ever synthesized. OpenFOAM was run at `Ha` in one, five, ten, twenty, fifty; the resolved MUG
channel was run at `Ha` in two, three, five; `Ha = 5` is the only value with both codes.

| Ha | u_mean/u_c (analytic) | mhdFoam L2 | mhdFoam Linf | mhdFoam u_mean/u_c | mhdFoam relerr | MUG L2 | MUG Linf | MUG u_mean/u_c | MUG relerr |
|----|----|----|----|----|----|----|----|----|----|
| 1  | 0.6774 | 1.045e-03 | 1.080e-03 | 0.6783 | 1.298e-03 | —         | —         | —      | —         |
| 2  | 0.7055 | —         | —         | —      | —         | 1.328e-04 | 1.843e-04 | 0.7056 | 7.174e-05 |
| 3  | 0.7420 | —         | —         | —      | —         | 2.089e-04 | 3.141e-04 | 0.7421 | 1.229e-04 |
| 5  | 0.8109 | 1.951e-04 | 2.944e-04 | 0.8118 | 1.079e-03 | 1.532e-04 | 4.816e-04 | 0.8108 | 1.491e-04 |
| 10 | 0.9001 | 6.902e-04 | 1.195e-03 | 0.9012 | 1.282e-03 | —         | —         | —      | —         |
| 20 | 0.9500 | 5.793e-04 | 1.464e-03 | 0.9511 | 1.156e-03 | —         | —         | —      | —         |
| 50 | 0.9800 | 5.971e-04 | 2.485e-03 | 0.9810 | 1.067e-03 | —         | —         | —      | —         |

At `Ha = 5`, the one value with all three sources, the discriminating scalar is 0.8109 analytic,
0.8118 from `mhdFoam`, and 0.8108 from MUG. The corresponding overlay, with a residual inset
showing agreement at the level of one part in ten thousand, is `results/report/threeway_overlay_Ha5.png`.

The resolved MUG leg stops at `Ha = 5` in this exercise for a meshing reason, not a solver reason.
The two-dimensional cube mesh used here is uniform in the wall-normal direction, so resolving the
thinning Hartmann layer requires adding cells everywhere. At `Ha = 10` the uniform mesh was pushed
to one hundred sixty cells across the channel, the implicit Newton-Krylov solve stiffened sharply,
and the run did not reach steady state within the available window. It was stopped rather than
reported as converged. The `mhdFoam` leg reached `Ha = 50` cheaply because it used a graded mesh
with a wall cell two orders of magnitude smaller than the channel half-width. Extending the
resolved MUG leg to higher Hartmann number is a wall-normal mesh-grading task, and is listed in
the next steps.

## 4. The Hunt conducting-wall duct (analytic and numerical reference only)

The conducting-wall case is Hunt's flow: a rectangular duct in a transverse field with perfectly
conducting Hartmann walls perpendicular to the field and insulating side walls parallel to it. Its
signature physics, and the feature that discriminates it from the insulating Hartmann channel, is
a pair of high-velocity jets near the insulating side walls over a near-uniform core. This leg was
addressed with a reference solution only; no MUG conducting-wall run was performed.

Rather than transcribe Hunt's closed-form double series, which is error-prone to reproduce
correctly, the reference (`cases/hunt/hunt_duct.py`) is a finite-difference solution of the exact
inductionless MHD duct equations that Hunt solved: a momentum equation and an induction equation
coupled through the Hartmann number, with no-slip on all four walls, the induced field vanishing
on the insulating side walls, and its wall-normal derivative vanishing on the conducting Hartmann
walls. Solved on a grid at `Ha` equal to five, ten, twenty, fifty, it reproduces the side jets and
their growth. The mid-plane peak-to-centerline ratio rises from 1.008 at `Ha = 5`, which is below
the side-jet threshold and therefore essentially jet-free, to 1.793 at `Ha = 10`, 6.225 at
`Ha = 20`, and 9.556 at `Ha = 50`, with the peak migrating toward the side walls as the Hartmann
number rises. This is the textbook Hunt behaviour and is qualitatively distinct from the jet-free,
flat-core insulating channel of Section 3.

The honest limitations are these. The reference is a finite-difference solution on a grid of
moderate resolution, not cross-checked here against the analytic series or a published value, so
its absolute accuracy is unquantified and the side-jet ratios should be read as qualitative
evidence of the correct physics rather than as validated precise numbers. The reported sign of the
peak location is not physically meaningful, since the two side jets are symmetric and the argument
of the maximum simply selects one of them. And no MUG comparison exists, because Hunt's duct
requires a different MUG geometry from the Section 3 channel, namely an out-of-plane axial flow
over the two-dimensional cross-section with walls on all four sides and a conducting Hartmann-wall
electromagnetic boundary condition, which the in-plane channel mapping does not provide. Building
that case is listed as the top next step.

## 5. Upgrade roadmap, ranked by effort to impact

The repository contains four planning notes. Two of them, a technical briefing and a LaTeX
roadmap, analyze the public `main` branch at an earlier commit and predate this branch; their
statements that no liquid-metal term exists anywhere describe `main`, not the Phase-1 work already
on this branch. The other two, an upgrade plan and a validation plan, are on-branch documents
that lay out the staged roadmap and the effort ranking summarized here and that record Phase 1 as
implemented and tested. The effort estimates below follow the upgrade plan's risk assessment; they
vary somewhat between the notes, for instance regional material support is rated more heavily in
the LaTeX roadmap, so they should be read as planning estimates rather than settled costs. With
that caveat, the staged plan and its effort ranking are as follows, ordered here by the ratio of
payoff to cost. Only Phase 1 is implemented. Everything in Phases 2 through 5 is a design proposal,
with no module or solver infrastructure yet in place.

The single highest-payoff upgrade is a low magnetic-Reynolds-number, inductionless,
electric-potential mode. In place of evolving the flux and field directly, the solver would solve
a scalar Poisson problem for the electric potential, set the current from Ohm's law, and apply the
resulting Lorentz force. For engineering liquid metals, where the magnetic Reynolds number is far
below one, this is both more appropriate and dramatically cheaper, since it removes the Alfvenic
timestep constraint and allows time steps larger by two to three orders of magnitude. The notes
identify this as the main computational unlock and as the change that would make MUG comparable in
cost to the established engineering codes. Its cost is correspondingly significant, at the
medium-high to high end of the effort scale, because it introduces a new field structure and a new
solve.

The conducting-wall boundary condition is the next most valuable, because it is the gateway to the
conducting-wall benchmarks. It would be implemented as a Robin condition on the flux,
`psi + c_wall d(psi)/dn = 0`, parameterized by a wall-conductance ratio `c_wall` equal to the wall
conductivity times thickness over the fluid conductivity times the duct half-width, whose limiting
values recover insulating and perfectly-conducting walls. A typed, default-off entry point for this
parameter, `c_wall = 0`, has been placed in the solver type as an inert, unused stub so the
conducting-wall rung has a home; the precise mapping of its limits to the two wall types still has
to be fixed against the flux sign convention, since the insulating Hartmann case on this branch
uses the natural zero-normal-derivative wall condition rather than a Dirichlet one. This is medium
effort, but it depends on the inductionless or incompressible substrate being in place first, and
it is what would let the Hunt and Shercliff side-layer cases be run in MUG and compared to the
reference of Section 4.

Regional and anisotropic material properties are a self-contained, medium-effort step that
generalizes the present uniform scalar coefficients to per-region values and a conductivity tensor,
reusing the regional-coefficient pattern already present in the three-dimensional solver. It is
validated naturally by an anisotropic manufactured solution and is a prerequisite for realistic
multi-material blanket geometries.

An incompressible velocity-pressure formulation is the foundational but heavier step. The relevant
liquid-metal regime is incompressible, which is architecturally incompatible with the present
compressible solver, so the recommended path is a new module with a mixed velocity-pressure
finite-element pair and a saddle-point solve, rather than retrofitting a pressure degree of freedom
into the existing module. This is medium-high effort and is the substrate that a high-Hartmann-number
resolved Hartmann benchmark would ultimately need.

The coupling to TokaMaker is the highest-effort and the most strategically distinctive item, and is
discussed on its own below. At present only a one-way export exists, from a TokaMaker equilibrium
into a MUG background field through an HDF5 snapshot. There is no reverse coupling of blanket
currents back onto the plasma equilibrium, no coupled timestep loop, and no Python interface to MUG
at all, the solver being Fortran-only. Building the bidirectional loop is high effort and gated on a
MUG scripting interface.

### The full-induction question and OFT's distinctive niche

A word on positioning, which is the author's strategic read rather than a claim drawn from the
repository notes. The standalone liquid-metal blanket problem in its transient, higher
magnetic-Reynolds-number regime, of order one to ten, where induced fields are not negligible and
where disruption and edge-localized-mode transients matter, is increasingly being taken up with
full-induction MHD in dedicated codes. MUG can run in that full-induction regime today, since it is
a full magnetic-Reynolds-number compressible solver, and that is genuinely its differentiator
against the low magnetic-Reynolds-number, incompressible engineering codes such as the HIMAG and
PTOLEMY families. The argument of this roadmap, however, is that competing on the standalone duct is
not where OFT is most distinctive. OFT's distinctive niche is the coupled plasma-blanket problem.
The background field that drives the blanket flow is itself produced by an evolving plasma
equilibrium, and OFT already carries that equilibrium in TokaMaker on a shared mesh and shared
input-output infrastructure. A bidirectional TokaMaker-to-MUG-and-back loop would let the blanket
respond to a time-dependent plasma, and the plasma respond to the induced blanket currents, which
the standalone full-induction duct codes have no mechanism to address. The full-magnetic-Reynolds-number
capability is most valuable kept for the plasma side of that coupled problem, while a cheaper
inductionless mode covers the engineering duct interior. This positioning is the reason the
low-magnetic-Reynolds-number mode and the TokaMaker coupling sit at the two ends of the roadmap: one
makes the duct interior cheap, the other realizes the coupled problem that is OFT's reason to be in
this space at all.

## 6. Honest scope and next steps

This is a personal verification and engineering-readiness exercise. It reproduces classical
Hartmann, and qualitatively Hunt, duct-flow results as an independent cross-check against an
analytic referee and against OpenFOAM, in order to learn the workflow and to scope upgrades. It
does not validate MUG, does not advance its physics, and does not claim new capability. The reduced
wall-drag closure was disabled in every run reported here; the resolved channel does not need it,
which was the point of resolving it. What is solidly established is the analytic referee, the
OpenFOAM-to-analytic agreement across the full Hartmann sweep at the one-part-in-a-thousand level,
and the resolved-MUG-to-analytic agreement at Hartmann numbers two, three, and five with the
parameter scaling confirmed to one part in a thousand. What is not established, and is labeled as
such, is the resolved MUG leg above Hartmann number five, any MUG conducting-wall result, and the
absolute accuracy of the Hunt reference.

The next steps, in order, are to add wall-normal mesh grading to the two-dimensional MUG mesh so
the resolved Hartmann leg reaches Hartmann numbers of twenty to fifty without a uniform-mesh
blow-up in cost; to build the conducting-wall MUG geometry, an out-of-plane axial-flow duct with a
conducting Hartmann-wall boundary condition, and compare it to the Hunt side-jet reference; to
implement the low-magnetic-Reynolds-number inductionless mode that makes those runs cheap; and to
calibrate the reduced wall-drag coefficient against these resolved results, with its normalization
ambiguity resolved first. The strategic fork that most shapes all of this, and the question I would
most like Yuchen's view on, is in the accompanying email note.

## Appendix: artifact map

All numbers in Section 3 were recomputed during packaging directly from the stored profile CSV
files by `cases/hartmann/recompute_headline.py`, and the figures were rebuilt from the same CSV
files by `cases/hartmann/plots_from_csv.py`, with no dependency on re-running OpenFOAM. The
evidence is collected under `results/`:

- `results/error_table.md`: the complete three-way error table, including the integer-Hartmann and
  fitted-Hartmann MUG columns.
- `results/report/threeway_overlay_Ha5.png`: the three-way overlay at `Ha = 5` with the residual inset.
- `results/report/umean_uc_vs_Ha.png`: the discriminating scalar against Hartmann number, all sources.
- `results/report/error_vs_Ha.png`: profile L2 error against analytic for both codes.
- `results/report/flowrate_vs_Ha.png` and `results/flowrate_vs_Ha.png`: the analytic flow-rate suppression.
- `results/overlay_Ha*.png`: per-Hartmann-number profile overlays, including the MUG-only cases at two and three.
- `results/profile_Ha*.csv`: the per-Hartmann-number normalized profiles (analytic, mhdFoam, MUG where present).
- `results/mug_profiles/Ha_{2,3,5}_profile.csv`: the resolved MUG channel profiles.
- `results/hunt/`: the Hunt side-jet table and the mid-plane and contour figures.
