# Scoping: a coupled plasma ↔ liquid-metal-blanket loop for OFT/MUG

*Purpose: decide, before building, what physically couples to what, at what order of
magnitude, and therefore whether the loop should be one-way, weakly two-way, or fully
bidirectional — and how it maps onto OFT's existing pieces. Grounded in the LM-blanket
MHD literature (sources at the end) and this project's validated MUG capabilities.*

---

## 1. The headline answer

**The coupling is strongly asymmetric and multi-channel, so the right build is *staged and
hierarchical*, not "fully bidirectional" from day one.**

- **Plasma → blanket is the dominant, order-1, essentially *imposed* direction.** The
  plasma's confining field sets the entire blanket MHD (Hartmann number ~10⁴). Almost all
  of the useful physics is captured one-way.
- **Blanket → plasma is a *weak* feedback in steady operation** — order `Rm ≈ 0.1`, the
  standard "inductionless" small parameter — and is routinely neglected in the blanket-MHD
  literature.
- **The feedback only becomes order-1 (genuinely two-way) in *fast transients*** — field
  ramps and disruptions on millisecond timescales — which is exactly the regime MUG's
  *full-induction* capability (the thing we just validated) exists for.

So: build **one-way first** (it is most of the physics and most of the existing
infrastructure), add a **loose two-way** feedback as a converged sub-iteration (cheap
because the feedback is ~10%), and reserve **tight bidirectional coupling** for the
transient/disruption problem where it is actually required.

---

## 2. The three coupling channels (they are not one thing)

"Plasma–blanket coupling" is really three physically distinct channels with different
directions, magnitudes, and timescales. Conflating them is the main way to get the
architecture wrong.

| # | Channel | Plasma → blanket | Blanket → plasma | Timescale |
|---|---|---|---|---|
| **A** | **Electromagnetic (field/current)** | confining field `B(R,Z)` sets the MHD — **O(1), dominant** | induced field `δB/B ~ Rm ~ 0.1` (steady); O(1) in disruptions | slow (equilibrium) / ms (transient) |
| **B** | **Plasma-facing surface loads** | heat, particle, current, momentum flux to the first wall / free surface — **O(1) for thermal/surface** | evaporated-metal impurity influx (chemistry, separate) | ~ms–s (thermal) |
| **C** | **Transient EM / disruption** | plasma-current decay drives eddy currents in blanket+structure — **O(1)** | eddy-current forces (MN-scale) + structure motion intercepting flux — **two-way** | ms (disruption) |

**Channel A** is the LM-MHD coupling proper — the field the fluid moves through, and the
field the fluid's currents make. **Channel B** is what `plasma_interface.py` scaffolds
(sheath heat/particle/current/momentum loads, and the liquid-surface thermal/Marangoni
/capillary response). **Channel C** is the disruption electromagnetics problem, which the
fusion community already solves with plasma-equilibrium → structure-eddy-current tools.

A subtlety worth stating: in the blanket-MHD literature the word "coupling" most often
means **channel-to-channel** electromagnetic coupling *inside* the blanket (leakage currents
across shared conducting walls between adjacent ducts), not plasma↔blanket. MUG already has
the ingredient for that (the validated Phase-4 conjugate-wall Robin BC); it is a *within*-blanket
effect, not part of the plasma loop, and should not be confused with it.

---

## 3. The magnitudes — "what order is the effect?"

Representative ITER/DEMO PbLi outboard-blanket numbers (`B ≈ 4 T`, `L ≈ 0.1 m`,
`u₀ ≈ 0.1 m/s`, PbLi properties):

- **Hartmann number `Ha ≈ 9000`** (≈8800 measured for the EU-WCLL frontal region). The
  field utterly dominates viscosity; boundary layers are ~`1/Ha ≈ 10⁻⁴` of the duct.
- **Interaction parameter `N = Ha²/Re ≈ 0.1–100`.** Lorentz forces dominate inertia.
- **Magnetic Reynolds number `Rm ≈ 0.1`** (steady). This is *the* number that sets the
  blanket→plasma feedback: the induced field is ~`Rm` of the applied field.

Consequences for the loop:

- **Steady:** `Rm ≈ 0.1` ⇒ induced field ≲10% of applied ⇒ **inductionless / quasi-static is
  the standard, and one-way is physically justified.** The blanket is a (near-)passive
  recipient of the plasma's field.
- **Transient (ms field variation, disruptions):** the induction term can no longer be
  dropped — the field "cannot keep up" and dynamically couples to the flow. This is the
  explicit motivation for the new full-induction blanket solvers (e.g. Vertex-CFD, 2025):
  *"traditional blanket modeling during normal steady operation commonly employs the
  inductionless approximation … [full induction is needed] when the plasma-confining
  magnetic field varies on millisecond time scales."* **This is MUG's niche.**

**So the ordering is: plasma→blanket `O(1)` always; blanket→plasma `O(Rm)≈0.1` steady,
`O(1)` in ms transients.**

---

## 4. Should it be fully bidirectional? — a regime map

| Regime | Induction | Coupling needed | Why |
|---|---|---|---|
| **Steady / slow operation** | inductionless (`Rm≪1`) | **One-way** (plasma field → blanket) | feedback is a ~10% perturbation, conventionally neglected |
| **Precision equilibrium / error-field** | inductionless + feedback | **Loose two-way** (converged sub-iteration) | the ~10% `δB` shifts the equilibrium at the percent level |
| **Field ramps / pulsed operation / disruptions** | **full induction** | **Tight two-way** | eddy/induced currents are O(1); structure & plasma co-evolve on ms |

**Verdict: do not commit to full bidirectional coupling globally.** It is *unnecessary* for
the bulk of blanket design (steady) and *essential* only for the transient regime. A staged
loop matches the physics and is far cheaper and more robust.

---

## 5. Map onto OFT — what already exists vs. what to build

OFT is unusually well-positioned because the *plasma ↔ conducting-structure* half of this
problem is already solved in-toolkit:

- **TokaMaker** — time-dependent Grad–Shafranov (axisymmetric plasma equilibrium). Already
  couples to ThinCurr; hybrid 2D–3D reconstructions exist. **Already has a one-way MUG field
  export (`save_mug()`)** + point-wise field interpolator (this project's roadmap).
- **ThinCurr** — 3D thin-wall eddy-current code for *solid* conducting structures; already
  inductively coupled to coils/plasma/changing plasma current.
- **MUG (`xmhd_2d`)** — the *fluid* blanket solver, now validated (this project): correct
  full induction, conjugate-wall Robin BC, stable to physical `Pm=10⁻⁶`, correct transient
  Alfvén dynamics.
- **`LM_MHD` coupling module** — the blanket→plasma field machinery: Biot–Savart from a
  current distribution (validated here to machine precision), regime metrics, sheath/surface
  loads.

**The missing piece is exactly one edge of the graph: the current-carrying *fluid* blanket
(MUG) as a participant in the plasma-equilibrium / eddy-current loop.** ThinCurr handles
*solid* conductors; MUG handles the *fluid* whose flow-driven currents ThinCurr cannot
represent. Everything else is either present or validated.

```
        TokaMaker  ──(field, exists: save_mug)──►  MUG (blanket fluid)     ◄── Stage 1
        (plasma eq)  ◄─(δB feedback, Biot-Savart, ─  currents             ◄── Stage 2
                        module validated)              
             ▲                                     
             └───── ThinCurr (solid-structure eddy currents, exists) ──────◄── Stage 3
                    ◄── MUG fluid currents into the same eddy-current bus ──►
```

---

## 6. Recommended staged build

**Stage 1 — one-way plasma → blanket (the 80/20; partly exists).**
- TokaMaker equilibrium field → MUG background `B₀(R,Z)` (extend the existing `save_mug`
  export + interpolator to feed MUG's `B_0` field over the blanket domain).
- Requires: an **axisymmetric (R–Z) MUG blanket case** (MUG already has cylindrical R–Z
  support) driven by a spatially-varying imposed field rather than a uniform `B₀`.
- Captures: MHD pressure drop, flow redistribution, Hartmann/side layers in the real field —
  i.e. essentially all of steady blanket design.
- **Validation:** against a known duct in a uniform field (recover the ladder), then against
  a published graded-field blanket case (Smolentsev/Bühler benchmarks).

**Stage 2 — loose two-way feedback (blanket → plasma, weak).**
- Compute the blanket's induced field `δB` from MUG's currents via the *validated* coupling
  Biot–Savart, project onto the plasma boundary, and re-solve the TokaMaker equilibrium;
  Picard-iterate to convergence.
- Cheap and robust: because `δB/B ~ Rm ~ 0.1`, the fixed-point converges in ~1–2 iterations.
- **Validation:** self-consistency (converged `δB`), and a manufactured case where the
  feedback magnitude is known; check the equilibrium shift scales as `Rm`.

**Stage 3 — tight two-way transient (the hard, novel regime).**
- Co-evolve MUG (fluid + induced currents, *full induction*) with ThinCurr (solid structure
  eddy currents) and TokaMaker (evolving plasma) on the ms timescale of a field ramp or
  disruption. This is where MUG's validated full-induction transient physics (the V2
  capstone regime) is *load-bearing*.
- Requires a shared eddy-current "bus" so the fluid and solid currents see each other's
  fields; time-stepping reconciled across the plasma/structure/fluid timescales.
- **Validation:** a rapid-field-decay channel problem (there are recent analytic/DNS
  references, e.g. "oscillatory LM flow under rapidly decaying applied field," 2026) — which
  is a clean, MUG-shaped transient two-way benchmark.

**Channel B (surface loads) is an independent one-way path** that can be built in parallel:
`plasma_interface.py` already provides the formula-level sheath and surface-response loads
(validated here); wiring them as boundary conditions on a free-surface/first-wall MUG run is
separate from the electromagnetic loop and does not need the plasma solver in the inner loop.

---

## 6a. Novelty of the Stage-3 target (plasma + solid + liquid-metal, self-consistent, transient)

Assessed against the current literature, honestly:

- **The liquid-metal *fluid* response to disruption-rate field changes is an active research
  frontier, not unprecedented.** A handful of 2025–2026 works (arXiv:2502.02699,
  arXiv:2606.25167; the SciDAC-5 plasma–LM effort; Maxwell–Navier–Stokes liquid-wall
  modeling arXiv:2411.11865) treat the LM as a *fluid* with disruption-driven eddy currents,
  Lorentz forces (`~10⁷–10⁸ N/m³`, accelerations `10³–10⁴ m/s²`) and flow response. So MUG
  doing transient LM-MHD is **timely and competitive**, not first-of-kind.
- **But that existing work almost universally *prescribes the field one-way*** — a disruption
  scenario supplies an imposed, time-varying field, and the LM responds to it. The field is
  an input, not a self-consistent participant.
- **The genuine novelty is closing the loop: a self-consistent, three-way
  plasma-equilibrium ↔ solid-structure-eddy-current ↔ liquid-metal-*fluid* coupling** — where
  the moving liquid metal's own currents feed back into the field seen by the plasma and the
  structure — in a single integrated open-source toolkit. The pieces exist *separately*
  (TokaMaker↔ThinCurr is the plasma↔solid loop; MUG is the validated fluid); **no one appears
  to have closed the full three-way loop with the blanket as a live fluid rather than a
  static conductor.** That integration is the contribution.
- **MUG walks in with the load-bearing physics already validated.** The literature identifies
  *standing Alfvén waves* as the key mechanism controlling the LM disruption response — and
  the conjugate (conducting-wall) coupling sets its damping. That is exactly the physics our
  V2 capstone validated to a few percent (transient Alfvén overshoot/timing, `Pm=1`
  full-induction, thin-wall Robin BC, with the wall-diffusion-time caveat). So MUG is
  arguably the **best-validated tool for the specific physics that governs this regime** — a
  strong, defensible position for the Stage-3 paper.

**Verdict:** Stage 3 is a genuine, publishable research contribution. Its novelty is
concentrated in the *self-consistent integration* (plasma + solid + fluid, closed loop),
which is unclaimed, rather than in the fluid-response physics alone (which is contested
frontier). Reaching it through Stages 1–2 is the right way to earn it — each stage is
independently useful and de-risks the next.

## 7. Recommendation and the decisions that gate it

**Recommendation: build Stage 1 first and treat it as the deliverable milestone.** It is the
dominant physics, extends infrastructure that already exists (`save_mug`), needs no new
plasma-side coupling, and is directly validatable against the blanket-MHD literature. Stage 2
is a small, high-value add (converged feedback) once Stage 1 is trustworthy. Stage 3 is a
research contribution in its own right and should be scoped separately once 1–2 are solid.

**Open decisions to make before coding Stage 1:**
1. **Geometry:** axisymmetric R–Z blanket sector (matches TokaMaker/ThinCurr) vs. a
   Cartesian duct module. The loop wants R–Z; MUG supports it but our validated cases are
   Cartesian — so the first task is an R–Z MUG blanket verification against a known case.
2. **Field ingestion — CONFIRMED, this is the first code change.** MUG's `B_0` is today a
   *constant static 3-vector* (`xmhd_2d.F90:91, REAL(r8) :: B_0(3)=0`), read once per cell
   (`B_0 = self%parent_sim%B_0`) and used in the flux and stretch terms. Stage 1 needs it
   promoted to a **spatial field `B₀(x,z)`** sampled at the integration points; Stage 3 needs
   it **time-varying `B₀(x,z,t)`**. This is a well-scoped, localized change (the sites that
   read `B_0` are the residual/Jacobian flux and stretch blocks — lines ~602/683, 801/886,
   1046, 1170/1256, 1544/1586) and is the concrete task-one for the whole build.
3. **Steady vs transient target first:** Stage 1 is steady (inductionless-equivalent, but
   MUG runs full induction at low `Rm` — we proved it is stable there). Confirm we want the
   steady milestone before the transient one.
4. **Feedback tolerance (Stage 2):** what `δB/B` accuracy matters for the intended use
   (equilibrium design vs. error-field studies) — sets whether Stage 2 is even needed.

---

## Sources
- Numerical investigation of moderate-`Rm` fusion LM-MHD flows — escholarship.org/uc/item/6373r5xw
- Full-induction MHD solver for LM fusion blankets (Vertex-CFD), arXiv:2511.15549 — motivation for full induction at ms field variation
- Electro-magnetic flow coupling for LM blanket applications — ScienceDirect S0920379615303811 (channel-to-channel coupling)
- Coupling of electromagnetics and structural/fluid dynamics — dual-coolant blanket under plasma disruptions — INIS yxpcw-hvb86
- Smolentsev et al., MHD code V&V for fusion blankets; HIMAG / characterization of MHD pressure drop — UCLA/ORNL
- TokaMaker (time-dependent Grad–Shafranov), ScienceDirect S0010465524000341; ThinCurr (3D thin-wall eddy currents), arXiv:2412.14962; Design of passive/structural conductors, IOP 10.1088/1741-4326/ad0bcf
- ITER/DEMO PbLi blanket parameters (Ha≈9000, N=0.1–100, Rm≈0.1) — MDPI Energies 14/19/6192 & 14/20/6640; UCLA Smolentsev 2008 FED
- Disruption EM loads / eddy currents on in-vessel components — ITER-FEAT VV (S092037960100494X), K-DEMO (S0920379621003689), KSTAR (S0920379623001357)
- Oscillatory LM flow under rapidly decaying applied field, arXiv:2606.25167 / 2502.02699 (transient two-way benchmark candidates)
