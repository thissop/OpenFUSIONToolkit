# Coupled-workflow validation plan: Plasma → ThinCurr → MUG

*How to earn confidence in the **composed** transient loop (TokaMaker plasma → ThinCurr solid
eddy currents → MUG fluid-blanket MHD), not just each code alone. Three individually-validated
codes do NOT make a validated coupled system — coupling has its own failure modes. This plan is
the methods backbone for the "reach" paper (self-consistent fluid-in-the-loop, candidate #3).*

---

## 0. Why this is a separate problem
Unit validation (each code vs analytic/benchmark) cannot see coupling-specific failures:
**interface transfer errors** (field/current mis-mapped between meshes/codes), **operator-splitting /
time-lag errors** (how often codes exchange data), **feedback instabilities** (a two-way loop can be
numerically unstable even when each piece is stable), and **conservation violations** across code
boundaries. So the composition needs its own validation layer.

## 1. The coupling graph and its edges
Validate each **edge** in isolation before closing the loop:

| Edge | Physics | Status |
|---|---|---|
| Plasma → ThinCurr | equilibrium field/current → solid eddy currents | already coupled+validated in OFT (TokaMaker↔ThinCurr passive-structure papers) |
| **(Plasma/ThinCurr field) → MUG** | applied B delivered to the fluid solver (`save_mug()` / interpolator) | **new edge to validate** |
| **MUG → field (feedback)** | fluid current distribution → induced B (Biot–Savart) | machinery validated to machine precision (LM_MHD module); re-check *at the loop interface* |

The genuinely new graph edge is the **current-carrying *fluid* blanket as a loop participant** —
ThinCurr handles *solid* conductors, MUG handles the *fluid* whose flow-driven currents ThinCurr
cannot represent. Everything else is present or validated.

## 2. Staged validation ladder (unit → edge → one-way → weak two-way → tight two-way)
- **Stage A — Unit (done/doing):** MUG (Shercliff/Hunt/V2 + finite-c QS map), ThinCurr (eddy benchmarks),
  TokaMaker (equilibria). *Paper #1 lives here.*
- **Stage B — Each edge, known transfer:**
  - *Field→MUG:* push a **known analytic field** through the interface; MUG must receive it to round-off
    (tests transfer, decoupled from physics).
  - *MUG→field:* verify the fluid-current→Biot–Savart export at the exact interface the loop uses.
- **Stage C — One-way loop:** full chain (plasma field → ThinCurr → MUG) on a case with a *known* delivered
  field; MUG must match its **stand-alone** result under that field. Mismatch ⇒ the pipeline corrupts the input.
- **Stage D — Weak two-way (steady, O(Rm≈0.1)):** add feedback as a converged sub-iteration. Require
  (i) sub-iteration residual → 0 each step, and (ii) feedback *effect* is O(Rm) (~10%). Wrong magnitude ⇒ artifact.
- **Stage E — Tight two-way (transient/disruption, O(1) feedback):** the reach-paper regime; validate against
  an independent coupled reference (§4) + all §3 diagnostics.

## 3. Coupling-artifact detection (the diagnostics unit tests can't provide)
1. **Turn-off recovery (cheapest, strongest):** set feedback coupling = 0 ⇒ must **exactly recover** the
   validated stand-alone MUG result. Any deviation = coupling is contaminating the base solution.
2. **Interface conservation:** total current, magnetic energy, momentum continuous across the code boundary.
   A **jump at the interface** = transfer/interpolation error.
3. **Coupling-scheme (exchange-frequency) convergence:** refine how often codes exchange → answer must
   converge. Dependence on the exchange interval = operator-splitting (time-lag) error.
4. **One-way vs two-way delta:** the difference must equal the physically-expected O(Rm) steady / O(1)
   transient. Wrong *magnitude* flags an artifact.
5. **Symmetry/parity:** the coupled solution must respect the problem's symmetries (cf. the even-b_y/odd-v
   parity we used for the disruption source). A broken symmetry is a coupling-bug fingerprint.
6. **Coupled residual monotonicity:** the outer coupled-iteration residual should decrease monotonically;
   oscillation/growth = feedback instability, fix with under-relaxation or a tighter coupling scheme.

## 4. Confidence multipliers (extend what made paper #1 credible)
- **Manufactured / analytic coupled benchmark FIRST.** Before trusting either code on the real problem,
  verify the *composition* against ground truth: a reduced geometry with a semi-analytic coupled answer —
  e.g., Zikanov's standing-Alfvén conjugate-wall duct **driven by a prescribed current loop standing in for
  the plasma**, where the coupled response is tractable. MMS (add source terms to a chosen coupled solution)
  is the rigorous version and tests the coupling machinery independent of physics.
- **Independent *coupled* reference (the two-code move, applied to the loop).** COMSOL is genuinely
  multiphysics: it can host eddy-structure + conjugate-wall fluid + a prescribed plasma-current source in one
  model. Box 1 building a **COMSOL coupled model** on a reduced geometry gives the independent check on the
  OFT (TokaMaker+ThinCurr+MUG) loop that the 0.1–1.2% loads agreement gave paper #1. Two independently-coupled
  codes agreeing on the *loop* is the credibility standard for the reach paper.

## 5. Mapping to the papers
- **Paper #1 (finite-c QS map + loads):** Stage A only — publishable now, no coupling.
- **Reach paper #3 (self-consistent fluid-in-the-loop):** Stages B→E + §3 diagnostics + §4 references **are
  its methods section.** The validation plan is not overhead; the coupled-validated demonstration IS the paper.

## 6. Concrete near-term first steps (cheap, high-information)
1. **Edge B (field→MUG):** impose a uniform/known analytic field through the interface, confirm MUG reads it
   to round-off. (Days; no new physics.)
2. **Turn-off recovery harness:** a one-line diagnostic that runs the loop with feedback=0 and diffs against
   the stand-alone MUG result. Build it once; reuse at every later stage.
3. **Pick the reduced analytic coupled benchmark** (standing-Alfvén conjugate duct + prescribed current loop)
   and derive/scope its semi-analytic coupled answer — this becomes the composition's ground truth.

*Companion to `coupled_loop_scoping.md` (physics/ordering) and `next_paper_scoping.md` (paper ranking).*
