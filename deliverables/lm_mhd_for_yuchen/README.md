# OFT/MUG Liquid-Metal MHD: verification and roadmap package

This folder is a self-contained package on the OFT/MUG liquid-metal MHD work, prepared by
Thomas Kiker (Columbia University, tjk2147@columbia.edu) for Yuchen Jiang (ORNL).

It is a personal verification and engineering-readiness exercise on OpenFUSIONToolkit's
two-dimensional MHD solver MUG (`src/physics/xmhd_2d.F90`). It independently reproduces the
classical Hartmann and Hunt duct-flow benchmarks against an analytic referee and against
OpenFOAM, and it scopes a set of upgrades toward liquid-metal blanket modeling. It does not
validate MUG or claim new physics capability, and the reduced wall-drag closure was disabled
throughout.

## What is in this folder

- `OFT_LM_MHD_state_and_roadmap.md`: the single self-contained report. It covers what the MUG
  solver actually solves today, the Phase-1 reduced wall-drag closure and what was and was not
  validated about it, the independent three-way Hartmann verification (the headline), the Hunt
  conducting-wall reference, the upgrade roadmap ranked by effort to impact, an honest scope
  statement, and next steps.

- `SUMMARY_email.md`: a short note suitable for pasting into an email to Yuchen, ending with the
  one question whose answer most shapes the next step.

- `results/`: the evidence behind the report.
  - `error_table.md`: the complete three-way Hartmann error table.
  - `report/`: presentation-quality figures (the three-way overlay with residual inset, the
    discriminating scalar against Hartmann number, the profile error against Hartmann number, and
    the analytic flow-rate suppression).
  - `overlay_Ha*.png`: per-Hartmann-number profile overlays.
  - `profile_Ha*.csv`: per-Hartmann-number normalized profiles (analytic, mhdFoam, and MUG where
    a MUG run exists).
  - `mug_profiles/`: the resolved MUG channel profiles at Hartmann numbers two, three, and five.
  - `hunt/`: the Hunt side-jet table and figures.

## How the headline numbers were produced

Every number in the report's verification section was recomputed during packaging directly from
the stored profile CSV files, not copied from any earlier note. The analytic referee
(`cases/hartmann/analytic_hartmann.py`) was re-run and passes its self-test. The error table and
figures were rebuilt from the CSV files by `cases/hartmann/recompute_headline.py` and
`cases/hartmann/plots_from_csv.py`, with no dependency on re-running OpenFOAM. Where a source
produced no run at a given Hartmann number, the table entry is blank; no value was synthesized.

## One-line orientation to the underlying repository

The verification cases live in `cases/hartmann/` and `cases/hunt/` on branch `lm-mhd-upgrades`.
The solver under test is `src/physics/xmhd_2d.F90`. The MUG driving case for the resolved
Hartmann channel is `src/tests/physics/test_hartmann_mug.F90`. The reduced wall-drag closure is
covered by `src/tests/physics/test_drag_decay.F90`.
