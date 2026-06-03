# OpenFOAM `mhdFoam` Cross-Validation Plan for the OFT Reduced Wall-Drag Model

**Branch:** `lm-mhd-upgrades`
**Date:** 2026-06-02
**Status:** Plan only — OFT Phase 1 validation is complete (see `lm_mhd_validation_plan.md`);
OpenFOAM setup is the next step and is documented here for execution.

---

## 0. What this validates — and what it does NOT

The OFT Phase 1 feature is a **reduced volumetric momentum sink**:
```
S_u = -α · u_⊥ ,   u_⊥ = u - (u·b̂) b̂
```
It does **not** resolve the thin Hartmann boundary layer on the mesh. It models the
*bulk* effect of that layer as a linear friction. Therefore:

- ✅ **Fair comparison:** bulk velocity damping / decay rate, flow-rate reduction vs an
  effective friction, and the scaling of that friction with field strength.
- ❌ **NOT a fair comparison:** the detailed `cosh`-shaped Hartmann velocity profile across
  the channel. The reduced model has no boundary layer to compare pointwise.

`mhdFoam` (a *resolved*, low-Rm-capable incompressible MHD solver) provides the "truth"
Hartmann flow. From it we extract an **effective bulk friction** and check that the OFT
reduced model, run with the matching `drag_coeff`, reproduces the same bulk behaviour.

This is the bridge between "the drag term decays velocity at the right numerical rate"
(already shown by `test_drag_decay`) and "the drag coefficient means something physical."

---

## 1. Analytic Hartmann channel flow (the reference both codes target)

1-D fully-developed flow between insulating plates at `y = ±a`, transverse field `B₀ ŷ`,
streamwise pressure gradient `G = -dp/dx`, kinematic viscosity `ν`, conductivity `σ`,
density `ρ`:

```
Ha = B₀ a sqrt(σ / (ρ ν))                         (Hartmann number)
u_x(y) = (G/(ρ ν)) (a²/Ha) [ 1 - cosh(Ha y/a)/cosh(Ha) ] · (1/ (1 - tanh(Ha)/Ha))  (insulating)
```

Useful bulk quantities:
```
u_mean / u_mean(Ha=0)  = 3 (Ha - tanh Ha) / Ha³      (flow-rate reduction vs Poiseuille)
Hartmann layer thickness  δ ≈ a / Ha
```

For **Ha ≫ 1** the core is nearly flat and the friction is dominated by the Hartmann
layers. The effective bulk friction (momentum sink per unit mass balancing the pressure
gradient at steady state) is:
```
α_eff ≈ σ B₀² / ρ · (1/Ha) = ν Ha / a²       (leading-order, high-Ha)
```
This is the mapping the OFT `drag_coeff` should match. Document the exact constant from
the steady balance `G/ρ = α_eff · u_mean` once `mhdFoam` gives `u_mean(G, Ha)`.

---

## 2. Install OpenFOAM (Ubuntu, apt)

OpenFOAM is available in the distro repos on this machine (confirmed via `apt-cache`):

```bash
sudo apt-get update
sudo apt-get install -y openfoam openfoam-examples
# Locate the environment script (path depends on package version):
dpkg -L openfoam | grep -E 'bashrc|etc/bashrc'
source /usr/share/openfoam/etc/bashrc     # adjust to the path found above
which mhdFoam blockMesh                    # confirm tools on PATH
```

`mhdFoam` is the canonical incompressible MHD solver shipped with OpenFOAM and is exactly
the Hartmann-problem solver (the OpenFOAM tutorial `tutorials/.../mhdFoam/hartmann` is the
reference case).

If the apt `openfoam` package does not ship `mhdFoam`, fall back to the openfoam.org or
ESI (`openfoam.com`) distribution, both of which include it; do **not** switch to FreeMHD
for this first cross-check.

---

## 3. Reference case: `mhdFoam/hartmann`

```bash
source /usr/share/openfoam/etc/bashrc
cp -r $FOAM_TUTORIALS/incompressible/mhdFoam/hartmann ~/hartmann_case   # path may vary
cd ~/hartmann_case
```

Key inputs to set (in `constant/` and `0/`):
- `constant/physicalProperties` (or `transportProperties`): `nu`, `rho`
- `constant/physicalProperties`: `DB` (magnetic diffusivity) / `B0` field direction & magnitude
- `0/U`: no-slip at the Hartmann walls (`y = ±a`), cyclic in x
- `0/B`: imposed transverse `B₀ ŷ`, with `dB/dn = 0` at insulating walls
- `0/p`, `0/pB`: pressure / magnetic pressure BCs
- `system/blockMeshDict`: 2-D channel, `2a` wall-normal extent, **graded mesh** so the
  Hartmann layer `δ ≈ a/Ha` is resolved with ≥ 8–10 cells (critical for high Ha)
- `system/controlDict`: run to steady state; write `u_x(y)` profile

Run:
```bash
blockMesh
mhdFoam            # or `foamRun -solver incompressibleMHD` on newer versions
postProcess -func 'graphCell' -latestTime    # extract u_x(y) along the wall-normal line
```

Sweep `Ha ∈ {1, 5, 20, 50}` by varying `B₀` (keep `a, ν, σ, ρ` fixed). Record:
- `u_x(y)` profile (for the analytic `cosh` check)
- `u_mean` (for the flow-rate-reduction check)

### Validation step A (OpenFOAM vs analytic)
Confirm `mhdFoam` reproduces §1: profile RMS error < few % and `u_mean/u_Poiseuille`
matching `3(Ha - tanh Ha)/Ha³`. This establishes the trusted reference.

---

## 4. OFT reduced-drag comparison

The OFT reduced model cannot reproduce the `cosh` profile, so compare **bulk** quantities.

### Step B1 — calibrate `drag_coeff` from the resolved solution
From the `mhdFoam` steady state at each `Ha`, compute the effective bulk friction:
```
α_eff(Ha) = (G/ρ) / u_mean(Ha)
```
This is the physically-correct `drag_coeff` the reduced model should use.

### Step B2 — run OFT reduced model with that α
Two equivalent checks, both using the existing Phase 1 machinery:

**(i) Decay form (already wired):** an unforced uniform flow under drag decays as
`u(t) = u₀ exp(-α t)` (backward-Euler discrete `u_{n+1}=u_n/(1+α dt)`). `test_drag_decay`
already confirms the OFT solver integrates this exactly for the calibrated `α_eff`.
→ This shows the *time-constant* the reduced model imposes equals the Hartmann-layer
friction time-constant `1/α_eff` extracted from `mhdFoam`.

**(ii) Forced steady form (new small driver, optional):** add a uniform body force `f_x`
to a periodic OFT box with `use_wall_drag=.TRUE.`, `drag_bhat=ŷ`, `drag_coeff=α_eff`. At
steady state the OFT momentum balance gives `u_x = f_x/α_eff`. Check this equals the
`mhdFoam` `u_mean` for the same `f_x = G/ρ`. This is the most direct bulk cross-check.
(Requires a tiny driver analogous to `test_drag_decay.F90` with a constant momentum
source instead of an initial condition — not yet written; see §6.)

### Step B3 — scaling check
Plot `α_eff(Ha)` from `mhdFoam` against the analytic high-Ha law `α_eff ≈ ν Ha / a²`
(= `σ B₀²/(ρ Ha)`). Confirm the OFT reduced model, fed `α_eff(Ha)`, tracks the
`mhdFoam` flow-rate reduction across the sweep. This demonstrates the reduced model is a
defensible closure in the high-Ha regime and quantifies where it breaks down (low Ha,
where the core is not flat and a single bulk friction is a poor description).

---

## 5. Deliverable artifacts

- `mhdFoam` profiles `u_x(y; Ha)` and `u_mean(Ha)` for `Ha ∈ {1,5,20,50}`
- A table: `Ha | u_mean(mhdFoam) | α_eff | u_x(OFT)=f_x/α_eff | rel. diff`
- A plot: flow-rate reduction vs `Ha` — analytic, `mhdFoam`, OFT-reduced overlaid
- A short note on the validity envelope (high-Ha good; low-Ha bulk-friction breakdown)

---

## 6. Exact next steps (in order)

1. `sudo apt-get install -y openfoam openfoam-examples`; `source .../etc/bashrc`;
   `which mhdFoam`.
2. Copy and run the `mhdFoam/hartmann` tutorial; confirm it converges on this machine.
3. Add wall-normal mesh grading; sweep `Ha`; export `u_x(y)` and `u_mean`.
4. Validation step A: OpenFOAM vs analytic (§1). Gate: profiles within a few %.
5. Compute `α_eff(Ha)` (§4 B1).
6. Reuse `test_drag_decay` to confirm the OFT decay time-constant = `1/α_eff` (§4 B2-i).
7. *(Optional)* write `test_drag_forced.F90` — a periodic box with a constant `f_x`
   momentum source + `use_wall_drag` — and check `u_x = f_x/α_eff` vs `mhdFoam` `u_mean`
   (§4 B2-ii). This is the cleanest bulk cross-check and the recommended new test.
8. Produce the table/plot (§5) and write the validity-envelope note.

**Do not** proceed to FreeMHD until this OpenFOAM cross-check is complete and documented.
**Do not** start incompressibility implementation until this confirms the reduced model's
calibration is understood (the incompressible solver is the proper tool for the regimes
where the bulk-friction closure breaks down).
