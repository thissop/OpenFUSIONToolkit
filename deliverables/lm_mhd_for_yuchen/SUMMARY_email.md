Subject: OFT/MUG liquid-metal MHD: independent Hartmann verification and an upgrade-direction question

Hi Yuchen,

This is a short package on the liquid-metal MHD work on OpenFUSIONToolkit's two-dimensional MHD
solver, MUG. It is a personal verification and engineering-readiness exercise: I independently
reproduced the classical Hartmann channel three ways, against the analytic profile, against the
resolved MUG solver, and against OpenFOAM's mhdFoam, and I scoped the upgrades that would take MUG
toward blanket modeling. Nothing here is a new capability claim, and the reduced wall-drag closure
was off in every run, since the whole point was to resolve the Hartmann layer directly rather than
model it.

What was verified, one line per rung:

- Analytic Hartmann referee: re-run and passes its self-test (Poiseuille two-thirds limit, plug
  core, and the one-over-Hartmann flow-rate suppression).
- OpenFOAM mhdFoam against analytic: agreement at the one-part-in-a-thousand level across Hartmann
  numbers one, five, ten, twenty, fifty.
- Resolved MUG (drag off) against analytic: confirmed at Hartmann numbers two, three, five, with the
  fitted Hartmann number matching the predicted Hartmann number to one part in a thousand; higher
  Hartmann numbers are blocked only by the uniform mesh, a grading task, not by the solver.
- Hunt conducting-wall duct: a finite-difference reference reproduces the side jets and their growth
  with Hartmann number; this is a reference solution only, no MUG conducting-wall run was done yet.

The headline three-way error table and the figures are in the attached folder; the full narrative,
including what the solver does and does not do today and the ranked upgrade roadmap, is in
OFT_LM_MHD_state_and_roadmap.md.

The one question that would most shape the next step: should OFT aim MUG at the standalone
full-induction blanket duct, or stay focused on the coupled plasma-blanket problem through
TokaMaker, which is where I think OFT is most distinctive since standalone duct codes have no
mechanism for the evolving plasma background. Two concrete sub-questions follow from whichever way
you lean: which liquid-metal normalization and material parameters does the ORNL side work in, lead
lithium or lithium; and is the conducting-wall benchmark you care about Hunt, Shercliff, or a duct
with a specified wall-conductance ratio. Your answer determines which roadmap branch I build next.

Thanks,
Thomas Kiker (Columbia University, tjk2147@columbia.edu)
