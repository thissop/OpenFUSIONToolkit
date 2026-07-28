# Next-paper scoping: ranked small-paper candidates (2026-07-27)

*PI tasking: find the robust, genuinely-novel result within our current capacity — the kind the
field is pointing at needing — for a small Smolentsev / Y. Jiang / Zikanov-audience LM-MHD paper.
Grounded in a 2024–26 literature scan with verbatim open-questions, mapped to our validated stack.*

---

## 1. What the field itself says is missing (verbatim, 2025–26 frontier)

Transient/disruption LM-MHD with **full induction** is the live frontier, and the leading groups
state their own gaps explicitly:

- **Vertex-CFD** — *A Full-Induction MHD Solver for LM Fusion Blankets* (Endeve, Slattery, **Smolentsev**
  et al., ORNL, arXiv:2511.15549, 2025): only *"initial verification,"* validated *against "recently
  published quasi-2D simulations"* (one code, quasi-2D), **finite-conductance wall BCs "left unspecified
  for now,"** disruption transients *"mentioned but not directly simulated,"* future = *"open avenues for
  future extensions… computational foundation for future simulations of transient MHD phenomena."*
- **Smolyanov & Zikanov** — *Exploratory study of LM response to rapid variation of applied field*
  (arXiv:2502.02699, 2025): **prescribed/imposed field** (not self-consistent), idealized walls, two-stage
  response, and flags *"the assumption of liquid metal incompressibility is found to be questionable."*
- **Oscillatory LM under rapidly-decaying field** (arXiv:2606.25167, 2026): standing Alfvén waves as the
  controlling mechanism — imposed field, channel (not conjugate-wall duct). This is *our validated V2
  physics*, published from the fluid side without the conjugate wall or the code cross-check.
- **Free-surface MHD beyond inductionless** (arXiv:2606.18745, 2026): the "beyond-inductionless" push is
  spreading to free-surface too — confirms the whole field is moving off the QS assumption.

**Net:** every frontier paper (a) imposes the field one-way, (b) uses idealized or unspecified walls, and
(c) has NO independent second-code cross-validation of transient loads. The **finite-c wall**, the
**quantified QS-validity boundary**, and **two-code robustness** are unclaimed — and we have all three.

**Novelty check on box 1's momentum-dip kernel:** a targeted search for *non-monotonic* QS/inductionless
validity in wall conductance returned **nothing** — the literature treats wall conductance monotonically.
A non-monotonic-in-`c` failure of the QS surrogate appears **unreported**.

---

## 2. Our capacity (validated, weeks-not-years)

Validated 2.5D finite-c COMSOL **and** MUG, cross-checked on transient disruption loads to **0.1–1.2%**
(blind pre-registered); the full **(Ha, c, τ_q) QS-validity map**; **Walker/inductionless
complementary-failure** result (both surrogates fail, in different corners ⇒ full model required);
Tier5 S&Z `c→0` machinery; conjugate-wall Robin BC; ILU+pre_freq solver upgrades; TokaMaker/ThinCurr/MUG
stack for the coupling reach. Known limit: pointwise near-wall *velocity* has an odd-even artifact in
`xmhd_2d` (loads unaffected; the 3D upwinded solver is the fix).

---

## 3. Ranked candidates (novelty × feasibility × Smolentsev/Jiang/Zikanov fit)

### ★ #1 — "When and why the quasi-static surrogate fails at disruption rates: a finite-conductance-wall validity map" — **RECOMMENDED, FASTEST**
- **The gap it fills (verbatim):** Vertex-CFD leaves finite-c walls "unspecified" and disruptions "not
  simulated"; Zikanov imposes the field with idealized walls. Nobody has mapped **QS-error vs (Ha, c, τ_q)**
  with a **conjugate wall**, and nobody reports the **non-monotonic-in-c** minimum.
- **Novelty:** HIGH — the c-dependence of QS validity, the non-monotonic dip (box 1's `c_dip ~ Ha⁻¹`
  mechanism if it closes), and a *two independent full-induction codes* certificate of the loads.
- **Feasibility:** HIGH — **data already in hand** (this project). 2.5D only; **needs no 3D**. Add box 1's
  mechanism closure + write-up. Weeks.
- **Audience fit:** direct hit — this is the "when is my cheap model safe?" question Smolentsev/Zikanov
  keep raising, answered quantitatively with the wall-conductance axis they've all omitted.

### #2 — "Two-code (COMSOL + MUG) cross-validated transient LM-MHD disruption benchmark with conjugate walls"
- **Gap:** Vertex-CFD's only cross-check is a quasi-2D code; there is no independent full-induction
  two-code transient benchmark, and none with finite-c walls. A blind two-code benchmark to ~1% is a
  V&V contribution the community (which runs on V&V) values.
- **Novelty:** MEDIUM-HIGH (V&V/benchmark, not new physics). **Feasibility:** HIGH (in hand). Overlaps #1;
  could be the methods half or a companion.

### #3 — "First self-consistent demonstration of a current-carrying *fluid* blanket in the plasma+solid eddy-current loop" — **THE REACH**
- **Gap:** every frontier paper imposes the field. The self-consistent 3-way loop (plasma equilibrium ↔
  solid eddy currents ↔ *fluid*-blanket induced currents) is **unclaimed** — Smolyanov&Zikanov 2025 is
  still one-way. Our scoping identifies the exact missing graph edge: MUG (fluid) as a participant
  alongside TokaMaker+ThinCurr (already coupled for solids).
- **Novelty:** HIGHEST. **Feasibility:** LOWEST — needs the coupling build (and 3D for a real geometry);
  weeks-to-months, and gated on the 3D velocity-gate robustness. The **3D Stage-1→3 plan feeds THIS**.
- Frame as a *reduced-problem first demonstration* to keep it small-paper scope.

### #4 — "Where thin-wall boundary conditions become unsafe for disruption transients" (from our thin-wall vs resolved-wall boundary map)
- **Gap:** thin-wall BCs are standard in blanket practice; nobody has quantified where they break under
  fast ramps. We measured it (MUG-thin vs COMSOL-resolved Iw: 0.7→5.2→13.6% across the δ_w/tw≈1 crossover).
- **Novelty:** MEDIUM (modeling-fidelity guideline). **Feasibility:** HIGH (in hand). A crisp letter, or a
  section of #1.

### #5 — "Incompressibility limits of transient LM disruption response" (probing Zikanov's own flagged limitation)
- **Gap:** Zikanov 2025 flags incompressibility as "questionable" (pressure waves) but doesn't resolve it.
  MUG's `true_pressure` compressible mode could quantify it. **Novelty:** MEDIUM. **Feasibility:** MEDIUM
  (needs the compressible-mode validation we parked). Lower priority.

---

## 4. Recommendation & implication for the 3D decision

**Write #1 first.** It is the highest novelty-per-week: the data exists, it needs no 3D, it answers the
exact question (finite-c QS validity + the non-monotonic-c mechanism) the frontier papers leave open, and
it lands squarely in the Smolentsev/Zikanov audience. #2 and #4 fold in as companion/sections.

**This reprioritizes 3D honestly:** the fast, high-novelty paper (#1) is 2.5D with data in hand — so the
3D build is **not** the critical path to a near-term paper; it is the enabling work for the **reach** paper
(#3, the self-consistent fluid-in-the-loop), which is the genuinely biggest but slowest result. So:
**#1 now (weeks) → 3D Stage-1→3 in parallel/after, targeting #3.** That sequences a guaranteed small paper
ahead of the ambitious one, which is the safer research bet.

Open item for box 1: whether its `c_dip ~ Ha⁻¹` (τ_w=T_A matched-damping) mechanism closes — if it does,
#1 gains a named mechanism, not just a mapped phenomenon, which is what elevates it from "useful map" to
"the result the field was pointing at."
