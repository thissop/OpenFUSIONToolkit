You are Codex inside my Lima Ubuntu VM (Ubuntu 25.10, aarch64). This is an UNATTENDED
overnight run — I am asleep. Never wait for my input. If you hit a decision point, make the
most reasonable conservative choice, log it, and keep going. Your goal is that I wake up to
the maximum amount of VERIFIED, HONESTLY-LABELED results plus a clear report of anything that
could not be completed and exactly why.

HARD RULES
- Do NOT start, stop, or modify any Lima VM. I already have one running. Work only in the
  current shell.
- You HAVE full install latitude: apt installs and source builds are both allowed.
- Be honest. Never report agreement you did not measure. If something fails, say so plainly
  and move to the next task. A truthful "this rung failed because X" is worth more than a
  faked pass.
- Keep a running log file at docs/dev/lm_mhd_verification_overnight_log.md — append a
  timestamped line at every major step (start, build done, case set up, run done, fail, fix,
  fallback taken). I want to reconstruct the night from this file.

REPO / CONTEXT
- Repo: /Users/tkiker/Documents/GitHub/OpenFUSIONToolkit (Mac mount via virtiofs; build on Linux).
- Branch: lm-mhd-upgrades. Run `git status` and `git branch --show-current` first; log them.
- The standalone 2D MUG solver is src/physics/xmhd_2d.F90: compressible visco-resistive MHD,
  fields n, velx, vely, velz, T, psi, by; scalar coefficients eta, nu, chi, D_diff. It CAN
  resolve Hartmann layers directly — that is the entire premise of this exercise.
- A previous session already built OFT successfully on this VM and passed test_alfven_2d,
  test_sound_2d, test_drag_decay (21 passed / 36 skipped). The Phase 1 reduced-drag term
  (use_wall_drag, drag_coeff, drag_bhat) exists and is DEFAULT-OFF.

SCOPE OF THIS RUN — this is a personal VERIFICATION / learning exercise, NOT a capability
claim about MUG. The team has already done Hartmann/Shercliff benchmarking; I am reproducing
known-good results independently and cross-checking against OpenFOAM to learn the workflow.
Frame all output that way. Do NOT claim to have "validated MUG" or advanced its physics.

THE DRAG TERM STAYS OFF. use_wall_drag = .FALSE. everywhere in this exercise. The whole point
is that a resolved channel does NOT need the reduced-drag closure. Do not enable it.

================================================================================
TASK A (CORE — must finish): Insulating-wall Hartmann channel, three-way verification
================================================================================
Verify three independent solutions of the SAME problem against each other:
  (1) the analytic Hartmann profile,
  (2) resolved MUG (xmhd_2d.F90, drag OFF),
  (3) OpenFOAM mhdFoam.

PROBLEM DEFINITION (insulating walls — this is the classic Hartmann 1937 case):
  Fully-developed, pressure-driven, incompressible, steady flow between two parallel plates at
  y = -a and y = +a. Streamwise velocity u(y) in x. Uniform transverse applied field B0 in y.
  Insulating walls (no wall current return; net streamwise current = 0).
  Hartmann number Ha = B0 * a * sqrt(sigma / (rho * nu)), where sigma = electrical
  conductivity, rho = density, nu = kinematic viscosity.

ANALYTIC SOLUTION (verify everything against this — it is the referee):
  u(y) / u_mean is fixed by Ha. The closed form for the velocity profile is:

     u(y) = u_c * [ cosh(Ha) - cosh(Ha * y/a) ] / [ cosh(Ha) - 1 ]      ... (insulating)

  where u_c is the centerline velocity. Equivalently, normalized so that the profile is driven
  by a constant pressure gradient G = -dp/dx:

     u(y) = (G / (sigma * B0^2)) * [ 1 - cosh(Ha * y/a) / cosh(Ha) ]

  Key verifiable scalars (these are what you actually compare — they are unambiguous):
   - Profile shape: flat core + boundary layers of thickness ~ a/Ha. As Ha grows the profile
     flattens (plug-like core) and the layers thin.
   - Mean/centerline ratio:  u_mean / u_c  ->  1 - tanh(Ha)/Ha.
   - Flow-rate-vs-Ha scaling: for fixed G, the dimensionless flow rate Q* satisfies
        Q* = (1/Ha) * [ 1 - tanh(Ha)/Ha ]   (up to the constant prefactor G/(sigma B0^2)),
     i.e. Q ~ 1/Ha at large Ha (the classic Hartmann flow-rate suppression).
  Implement these analytic expressions in a Python reference module
  (cases/hartmann/analytic_hartmann.py) with a function u_profile(y, a, Ha) and
  flow_rate(Ha), with docstrings citing the formulas above. Unit-test the limits:
   - Ha -> 0 must recover the parabolic Poiseuille profile (check u_mean/u_c -> 2/3).
   - Ha large must give u_mean/u_c -> 1/Ha behavior (plug core).
  If these limits do not come out right, your analytic implementation is wrong — fix it before
  trusting any comparison.

DO A SWEEP: Ha in {1, 5, 10, 20, 50}. For each Ha, compute analytic profile + flow rate.

----- (2) Resolved MUG side -----
First DETERMINE whether standalone MUG can be driven into a steady Hartmann channel. It needs:
pressure-gradient / body-force forcing in x, no-slip at y=+/-a, periodicity or
developed-flow in x, a uniform applied B in y, insulating EM BCs, and a way to march to
steady state and extract u(y) across the channel.
  - Inspect existing tests/examples: src/tests/physics/*, any examples/ dirs, test_alfven_2d
    and test_sound_2d harnesses, and the input/config parsing in xmhd_2d.F90. Determine what
    forcing and BC options actually exist.
  - IF a driven-duct setup is feasible from existing machinery: build the case
    (cases/hartmann/mug/), run the Ha sweep to steady state, extract u(y), and compare to
    analytic (L2 and Linf error on the normalized profile; relative error on u_mean/u_c and on
    flow rate). Target agreement: profile within a few percent for resolved meshes; refine mesh
    near walls so the Hartmann layer (~a/Ha) is resolved by >= ~8 cells at the highest Ha.
  - FALLBACK (you have my permission — chosen explicitly): if MUG has no clean driven-channel
    path, ADAPT THE TEST HARNESS as needed — e.g. add a minimal body-force/pressure-gradient
    forcing option or a dedicated cases/hartmann test driver, MIRRORING existing code style,
    DEFAULT-OFF and non-invasive (must not change existing test results). If you add a forcing
    option to xmhd_2d.F90, gate it behind a new default-off flag, document it, and confirm
    test_alfven_2d + test_sound_2d still pass unchanged afterward. Log exactly what you added.
  - If after reasonable effort MUG still cannot be driven cleanly, STOP the MUG rung, document
    precisely why (what BC/forcing is missing, what you tried), and proceed — do NOT spin on it
    all night. Cap total MUG-setup effort at a sensible bound and move on so Task A's OpenFOAM
    and analytic rungs definitely complete.

----- (3) OpenFOAM mhdFoam side -----
  - Install OpenFOAM (apt: the `openfoam` package; if the apt route misbehaves, source build
    is permitted — you have full latitude — but prefer apt for speed). Log the install path and
    version. Source the OpenFOAM environment before running.
  - Build a Hartmann channel case (cases/hartmann/openfoam/): mhdFoam solves induction-form
    incompressible MHD (evolves U, p, and B with the Lorentz force). Set up a 2D channel
    (parallel plates), uniform applied B in y, pressure-driven flow in x, no-slip walls, and
    INSULATING EM BC on B at the walls (this is mhdFoam's natural case — it has NO
    wall-conductance model, which is fine for the insulating verification). Use a graded mesh
    resolving the Hartmann layer. Run each Ha in the sweep to steady state; extract centerline
    u(y) profile.
  - Compare mhdFoam u(y) to analytic (same error metrics). Note mhdFoam's limitations honestly
    in the writeup (minimal solver, insulating only).

----- Task A output -----
  - cases/hartmann/results/: per-Ha CSVs of y, u_analytic, u_mug (if available), u_mhdfoam.
  - cases/hartmann/results/overlay_Ha*.png: overlay plots, analytic vs MUG vs mhdFoam, one per Ha.
  - cases/hartmann/results/flowrate_vs_Ha.png: Q vs Ha, all three sources, with the 1/Ha
    reference line.
  - cases/hartmann/results/error_table.md: L2/Linf profile error and flow-rate error for each
    source at each Ha. This is the headline deliverable — the three-way agreement table.

================================================================================
TASK B (CHERRY — only if Task A is complete and time remains): Conducting-wall duct vs Hunt
================================================================================
This is MUG-vs-analytic ONLY. NO OpenFOAM. mhdFoam has no wall-conductance model and the
conducting-wall OpenFOAM solvers (epotFoam / FreeMHD lineage) are an install rabbit hole —
explicitly OUT OF SCOPE tonight. Do not attempt them.

PROBLEM: Hunt's flow — fully-developed flow in a rectangular duct with a transverse applied
field, with two PERFECTLY CONDUCTING Hartmann walls (perpendicular to B) and two INSULATING
side walls (parallel to B). Hunt (1965) gives the analytic solution. The signature physics:
high-velocity side jets near the insulating side walls (the "Hunt jets") and a near-uniform
core — qualitatively different from the insulating Hartmann case, and the discriminating
feature to verify.
  - Implement Hunt's analytic solution in cases/hunt/analytic_hunt.py (series solution; cite
    Hunt 1965, J. Fluid Mech. 21, 577). Verify it produces the characteristic side-jet profile
    at moderate-to-high Ha. If the series implementation is shaky, it is acceptable to verify
    against the well-known qualitative features (side jets present, conducting-wall flow-rate
    scaling Q ~ 1/Ha differing from insulating) rather than fake precision — log which.
  - IF the MUG side of Task A worked (driven duct available): set the Hartmann-wall EM BC to
    conducting and re-run at one or two Ha values; compare to Hunt (look for the side jets).
  - IF MUG could not be driven in Task A: do analytic Hunt ONLY, and document that the MUG
    conducting run is blocked for the same reason as Task A. Do not fake a MUG result.
  - Output: cases/hunt/results/ overlay plot(s) MUG vs Hunt (or analytic-only), and a short note.

================================================================================
WRAP-UP (always do this, even on partial completion)
================================================================================
Write docs/dev/lm_mhd_verification_results.md containing:
  1. One-paragraph honest summary of what was verified and what was not.
  2. Task A three-way error table (analytic / MUG / mhdFoam) — the headline.
  3. Which rungs completed, which were skipped/blocked, and the EXACT reason for each block.
  4. Any code added to MUG (the forcing option, if you added one) and confirmation that
     existing tests (test_alfven_2d, test_sound_2d, test_drag_decay) still pass unchanged.
  5. Honest scope statement: this reproduces known Hartmann/Shercliff benchmarking as a
     personal verification exercise; it does NOT validate MUG or claim new capability; the
     reduced-drag term was OFF throughout (resolved-layer regime).
  6. Clear next steps.
Also keep docs/dev/lm_mhd_verification_overnight_log.md (the timestamped running log).
Do NOT git commit anything — leave changes uncommitted on the branch for me to review.

Priority order if time is short: Task A analytic + mhdFoam (guaranteed completable) >
Task A MUG > Task B analytic Hunt > Task B MUG. Always leave me a coherent results doc.
