# LM-MHD references

Curated literature for the MUG liquid-metal MHD validation (Paper 2) and the
coupled plasma-blanket scoping. PDFs are in this directory; online-only sources
are listed with URLs. Organized by theme, with a one-line note on why each matters
to *our* work.

---

## 1. Steady LM-blanket MHD and the inductionless / electric-potential formulation

This is how the field community actually solves spatially-varying (fringing)
fields: **do not evolve the field** — solve for the electric potential with the
imposed `B(x)` as a coefficient. Directly relevant because it is the alternative
formulation to MUG's full-induction approach, and it uses the *same* conjugate-wall
BC we implemented.

- **`AlbetsChico2011_fringing-field-consistency.pdf`** — Albets-Chico, Votyakov,
  Radhakrishnan & Kassinos, *Fusion Eng. Des.* **86**, 5-14 (2011). DNS of LM flow
  in a fringing field. Governing set is the **inductionless** one:
  `∂u/∂t + u·∇u = −∇p + (1/Re)∇²u + N(−∇φ + u×B)×B`, `∇·u = 0`,
  `∇²φ = ∇·(u×B)`. Wall BC is the **thin-wall Walker/Robin condition**
  `∂φ/∂n|_wall = −∇_τ·(c ∇_τ φ)`, `c = σ_w t_w/(σ L)` — the *same* conjugate-wall
  conductance ratio as our Phase-4 BC. Key result: using a **Maxwell-consistent
  (div-free AND curl-free) imposed field matters** — the historically neglected
  consistency explains the under-predicted peak transverse pressure gradient.
  *Why it matters to us:* confirms (a) the community sidesteps our full-induction
  `∇B_0`-terms problem by not evolving the field; (b) an imposed graded field must
  be div/curl-free; (c) our conjugate-wall BC is the standard one.
- **`ORNL_MHD-survey.pdf`** — ORNL/UCLA (Smolentsev et al.) survey of LM-blanket
  MHD modeling (past/present/future). The landscape: HIMAG, potential-based solvers,
  V&V benchmarks. *Why:* situates MUG among the established steady-state tools and
  their inductionless assumption.

## 2. Full-induction and transient LM-MHD — MUG's niche

The regime where the inductionless approximation fails and full induction is
required: **millisecond field transients / disruptions**. This is MUG's
comparative advantage (validated V2 transient conjugate capstone).

- **`VertexCFD_full-induction-LM-blanket.pdf`** — A full-induction MHD solver for
  LM fusion blankets (Vertex-CFD, 2025). The closest analog to MUG. Explicit
  motivation: *"traditional blanket modeling during normal steady operation
  commonly employs the inductionless approximation … [full induction is needed]
  when the plasma-confining magnetic field varies on millisecond time scales."*
  *Why:* independent confirmation that full-induction LM-MHD is a live, current
  research target — and a candidate cross-code comparison for MUG.
- **`Oscillatory-LM-rapidly-decaying-field.pdf`** — Oscillatory LM channel flow
  under a rapidly decaying applied field (2026). Treats the LM as a *fluid* with
  disruption-driven eddy currents and Lorentz forces; identifies **standing Alfvén
  waves** as the controlling mechanism. *Why:* this is (a) a clean, MUG-shaped
  transient benchmark and (b) direct evidence that the physics governing the LM
  disruption response is exactly what our V2 capstone validated (transient Alfvén +
  conjugate wall).

## 3. Adjacent — blanket / FWBS design context (not LM-MHD)

- **`Moreno2024_Frontiers-NuclEng.pdf`** — ParaStell: parametric FWBS modeling and
  neutronics for stellarator power plants (Moreno, Bader, Wilson, 2024). Stellarator
  first-wall/blanket/shield + OpenMC neutronics. *Why:* adjacent design-context, not
  LM-MHD; useful for the blanket-integration picture, not the MHD physics.

---

## Online sources (no local PDF)

**OFT coupling machinery (the plasma↔solid loop already exists):**
- TokaMaker (time-dependent Grad-Shafranov): ScienceDirect `S0010465524000341`
- ThinCurr (3D thin-wall eddy currents): arXiv `2412.14962`
- Passive/structural conductor design (TokaMaker↔ThinCurr): IOP `10.1088/1741-4326/ad0bcf`

**Spatially-varying / fringing field studies (inductionless):**
- 3-D MHD in a spatially varying solenoidal field: ScienceDirect `S0920379620302076`
- Heterogeneous-field duct flow: arXiv `0705.0633`
- Access-duct fringing-field pressure drop (Smolentsev): ScienceDirect `S0920379624001157`

**Disruption electromagnetics (the transient regime, solid-structure side):**
- ITER-FEAT VV/in-vessel EM loads: ScienceDirect `S092037960100494X`
- K-DEMO divertor EM loads (eddy + halo): ScienceDirect `S0920379621003689`

**Cross-code validation candidate:**
- OpenFOAM LM-MHD solver (unstructured, validated): ScienceDirect `S0920379615300880`
  — a possible *third* reference for Paper 2 alongside COMSOL + analytic series.
