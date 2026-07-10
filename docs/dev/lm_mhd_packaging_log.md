# LM-MHD Packaging Session — Running Log

Timestamped running log for the packaging / honest-verification session that bundles the
OFT/MUG liquid-metal MHD work into a single deliverable for Yuchen Jiang (ORNL).
Append-only. All times UTC.

## Step 0 — Provenance and inventory

- 2026-06-30T16:03:50Z  SESSION START. Mac-side mount (darwin); build tree is Linux (Lima VM, not driven this session).
- 2026-06-30T16:03:50Z  Provenance: branch=`lm-mhd-upgrades`  HEAD=`3e01e1358dd4b0f47ccedf71742894644a14a5b9`.
  `git status`: clean tracked tree; untracked = `.claude/`, `CLAUDE.MD`, `OFT_MUG_LiquidMetal_Briefing.md`, `mug_liquid_metal_wall_drag_plan.tex`.
  Recent commits: 3e01e13 (smplotlib report plots + MUG Ha={2,3,5}), 09ee98f, 23e57dc, 3a15548, 65ed685 (the overnight three-way Hartmann run).
- 2026-06-30T16:03:50Z  Inventory confirmed on disk: docs/dev/{lm_mhd_verification_results.md, lm_mhd_verification_overnight_log.md, lm_mhd_lima_handoff.md, lm_mhd_upgrade_plan.md, lm_mhd_validation_plan.md, lm_mhd_openfoam_validation.md}; root {OFT_MUG_Phase1_Run_Summary.md, OFT_MUG_LiquidMetal_Briefing.md, mug_liquid_metal_wall_drag_plan.tex}; cases/hartmann/{analytic_hartmann.py, mug/, openfoam/, results/} and cases/hunt/{hunt_duct.py, results/}.

### State table (Step 0.3) — verified against artifacts, not labels

| Rung | Status | Evidence path |
|------|--------|---------------|
| A1 analytic Hartmann | COMPLETED | `cases/hartmann/analytic_hartmann.py` — self-test re-run PASS this session |
| A2 resolved MUG (drag OFF) | COMPLETED at Ha∈{2,3,5} | `cases/hartmann/mug/Ha_{2,3,5}_profile.csv`, driver `src/tests/physics/test_hartmann_mug.F90`. Ha≥7 not finalized (uniform-mesh stiffness). |
| A3 OpenFOAM mhdFoam | COMPLETED at Ha∈{1,5,10,20,50} | `cases/hartmann/results/profile_Ha*.csv`, `cases/hartmann/openfoam/Ha_*/` |
| B analytic Hunt | COMPLETED (finite-difference referee, not closed-form series) | `cases/hunt/hunt_duct.py`, `cases/hunt/results/hunt_sidejet_table.md` |
| B MUG conducting | NOT ATTEMPTED | needs out-of-plane axial-flow geometry + conducting Hartmann-wall EM BC; flagged as next step in results doc §7 |

## Step 1 — Verify the verification (artifacts, not labels)

- 2026-06-30T16:03:50Z  Step 1.1 analytic referee: re-ran `python3 cases/hartmann/analytic_hartmann.py` → RESULT: PASS.
  Confirmed Ha→0 u_mean/u_c = 0.666667 (Poiseuille), profile-vs-parabola relL2 = 1.8e-8, Ha=200 u_mean/u_c = 0.995 (plug core),
  Q* ~ 1/Ha (rel 2e-2 at Ha=50, 5e-3 at Ha=200), eq.(4) vs numeric rel ≤ 1e-8 across sweep. Analytic module trusted.
- 2026-06-30T16:07:39Z  Step 1.2 headline regenerated DIRECTLY FROM CSVs via new `cases/hartmann/recompute_headline.py`
  (no meshio/VTK; mhdFoam+analytic from `results/profile_Ha*.csv`, MUG from `mug/Ha_{2,3,5}_profile.csv`). Findings:
  - mhdFoam recomputed L2/Linf/relerr reproduce the committed `error_table.md` to ~1e-7 absolute → published mhdFoam headline VERIFIED.
  - MUG: reconciled the two prior "conflicting" Ha=5 L2 values — `1.53e-4` is vs analytic at **integer** Ha=5; `1.25e-4`
    is vs analytic at the **fitted** Ha=4.998. Both correct; the rebuilt table now reports both columns explicitly.
  - Ha_fit/Ha_pred = 1.0003 (Ha2), 1.0004 (Ha3), 1.0007 (Ha5); u_mean/u_c = 0.7056/0.7421/0.8108 match the results doc §3.
  - Rebuilt `cases/hartmann/results/error_table.md` complete and honest: real MUG points at Ha∈{2,3,5}, blanks where no run,
    integer-Ha and fit-Ha L2 both shown. The old table silently omitted the Ha=2,3 MUG points; now included.
- 2026-06-30T16:07:39Z  Step 1.2 plots rebuilt from CSVs via new `cases/hartmann/plots_from_csv.py`:
  results/overlay_Ha{1,2,3,5,10,20,50}.png (added Ha=2,3 MUG-only overlays), results/flowrate_vs_Ha.png,
  and results/report/{threeway_overlay_Ha5, umean_uc_vs_Ha, error_vs_Ha, flowrate_vs_Ha}.png (smplotlib).
  Visual check of threeway_overlay_Ha5.png: analytic/mhdFoam/MUG coincide; residual inset ~1e-4. Good.
- 2026-06-30T16:07:39Z  Step 1.3 regression integrity: build tree `builds/build_release/` holds Linux ARM aarch64 ELF
  binaries (`file` check) — NOT runnable on this macOS shell, so the three tests cannot be re-executed here.
  Honest position: the overnight log records `test_alfven_2d`/`test_sound_2d`/`test_drag_decay` → 21 passed / 36 skipped,
  unchanged after the body-force addition; and the source change is inert-by-construction (use_wall_drag and use_body_force
  both default `.FALSE.`; body force enters the residual only, lines 967-977, with no Jacobian term ⇒ zero effect when off).
  Re-running on the VM is a one-command follow-up; not faked here.
- 2026-06-30T16:07:39Z  Confirmed `use_wall_drag` stays OFF: MUG driver `test_hartmann_mug.F90` sets `use_body_force=.TRUE.`,
  `body_force=[fx,0,0]`, insulating `psi` natural BC, no-slip walls, n/T frozen; explicit comment `! drag stays OFF`.
- 2026-06-30T16:08:00Z  Fact-extraction workflow (5 parallel grounded readers over xmhd_2d.F90, the 4 planning docs,
  the briefing, the .tex, and the Hunt/Phase-1 files). Key honesty flags captured: solver is fully COMPRESSIBLE (no
  incompressible/low-Mach mode; incompressibility imposed at case level via frozen n,T); the upgrade/validation plans are
  on-branch while the briefing/.tex analyze `main` at 5faae9a; the reduced-drag SM82 closure has an undefined `drag_coeff`
  normalization and a factor-Ha internal inconsistency in the .tex; Hunt was an FD referee (not the analytic series),
  grid-sensitive, MUG-Hunt not run.

## Step 2 — Deliverable produced

- 2026-06-30T16:16:00Z  Wrote `deliverables/lm_mhd_for_yuchen/`: `README.md` (orientation), `OFT_LM_MHD_state_and_roadmap.md`
  (the single self-contained report: what xmhd_2d.F90 solves today, the Phase-1 drag closure as already-validated
  background, the three-way Hartmann headline, the Hunt reference, the effort-ranked roadmap incl. the full-induction /
  coupled-plasma-blanket positioning, scope + next steps), `SUMMARY_email.md` (paste-ready note ending with the one
  question for Yuchen), and `results/` (verified error_table.md, overlay + report PNGs, per-Ha CSVs, MUG profiles, Hunt
  evidence). Style rules enforced: zero prose em dashes, no run-in bold headings, formal tone.
- 2026-06-30T16:24:00Z  Adversarial verification workflow (4 skeptic agents, one per report section, tasked to REFUTE
  claims against on-disk artifacts). Section 4 (Hunt) sound. Confirmed fixes applied to the report:
  (a) parallel-component drag error stated as ~5.3e-13 (was wrongly "machine precision ~1e-16"); cited recorded
  `drag_decay.results` value 5.17e-5 for the perpendicular component and dropped the unsupported "bounded by nl_tol";
  (b) mhdFoam L2 stated as "of order 1e-3, largest 1.045e-3 at Ha=1" (the flat "<=1e-3" bound was contradicted by the
  Ha=1 row); (c) developed-flow uniformity "about one part in a million" (was "better than", actual ~1.1e-6);
  (d) analytic eq.(4)-vs-numeric "of order 1e-8" (was "below 1e-8"; 1.06e-8 at Ha=50); (e) Section-5 provenance corrected
  to distinguish the on-branch upgrade/validation plans from the main-analyzing briefing/.tex, and flagged that effort
  estimates vary between notes (regional materials rated higher in the .tex). The verifier's one "missing path" flag
  (`results/mug_profiles/`) was a false positive: that dir exists inside the deliverable and the path is correct relative
  to the package.

## Step 3 — Optional item taken: conducting-wall Robin BC stub (default-off, inert)

- 2026-06-30T16:27:00Z  Added an inert, default-off typed entry point `REAL(r8) :: c_wall = 0.d0` to `oft_xmhd_2d_sim`
  in `src/physics/xmhd_2d.F90` (11-line diff, +stub +doc comment), mirroring the existing `drag_coeff` declaration. It is
  an UNUSED derived-type component (not read by setup_bc, the residual, or the Jacobian), so it is inert by Fortran
  semantics and cannot change any test result; type components do not trigger unused-variable warnings even under -Werror.
  HONEST CAVEAT: this Mac shell cannot build or run the Linux aarch64 test binaries, so I did NOT compile-check or re-run
  the tests after this edit. A VM rebuild + the three-test regression re-run is the required confirmation before relying on
  the stub. The accompanying design note (report Section 5) flags that the c_wall limit-to-wall-type mapping must be settled
  against the flux sign convention before wiring, since the insulating Hartmann case on this branch uses natural d(psi)/dn=0.

## Step 4 — Wrap-up

- 2026-06-30T16:27:25Z  HONEST SUMMARY. The deliverable `deliverables/lm_mhd_for_yuchen/` (+ `.zip`, 38 files) packages an
  independent, artifact-grounded verification of the insulating Hartmann channel plus an effort-ranked LM-MHD upgrade
  roadmap. VERIFIED FROM ARTIFACTS THIS SESSION: the analytic referee (self-test re-run PASS); mhdFoam-vs-analytic across
  Ha={1,5,10,20,50} at ~1e-3, recomputed directly from the stored CSVs and matching the prior committed table to ~1e-7;
  resolved-MUG-vs-analytic (drag OFF) at Ha={2,3,5} with Ha_fit/Ha_pred=1.0003/1.0004/1.0007 and profile L2 ~1e-4; the
  earlier "1.25e-4 vs 1.53e-4" MUG discrepancy reconciled as fit-Ha vs integer-Ha references. NOT ESTABLISHED (and labeled
  so): resolved MUG above Ha=5 (uniform-mesh stiffness, a grading task); any MUG conducting-wall run; the absolute accuracy
  of the Hunt FD referee (qualitative side-jet evidence only). The reduced wall-drag closure was OFF in every verification
  run.
- 2026-06-30T16:27:25Z  Task-A three-way error table as it now stands (from `cases/hartmann/results/error_table.md`):

  | Ha | u_mean/u_c (ana) | mhdFoam L2 | mhdFoam relerr | MUG L2 | MUG relerr |
  |----|----|----|----|----|----|
  | 1  | 0.6774 | 1.045e-03 | 1.298e-03 | —         | —         |
  | 2  | 0.7055 | —         | —         | 1.328e-04 | 7.174e-05 |
  | 3  | 0.7420 | —         | —         | 2.089e-04 | 1.229e-04 |
  | 5  | 0.8109 | 1.951e-04 | 1.079e-03 | 1.532e-04 | 1.491e-04 |
  | 10 | 0.9001 | 6.902e-04 | 1.283e-03 | —         | —         |
  | 20 | 0.9500 | 5.793e-04 | 1.156e-03 | —         | —         |
  | 50 | 0.9800 | 5.971e-04 | 1.067e-03 | —         | —         |

  (Ha=5 is the only three-way row: u_mean/u_c = 0.8109 analytic / 0.8118 mhdFoam / 0.8108 MUG.)
- 2026-06-30T16:27:25Z  DRAG STAYED OFF: `use_wall_drag=.TRUE.` appears only in `test_drag_decay.F90` (its own regression
  test), never in any Hartmann/Hunt verification case; the MUG channel is driven by the default-off body-force term instead.
  REGRESSION: the three tests were confirmed passing on the VM (overnight log 21 passed / 36 skipped; on-disk artifacts
  `builds/build_release/tests/physics/{alfven_2d,sound_2d,drag_decay}.results`, drag_decay.results = 5.17e-5 / 5.33e-13,
  binaries dated after the body-force addition). They were NOT re-run this session: the build tree is Linux aarch64 ELF and
  this is a macOS shell. The new c_wall stub post-dates those artifacts and is inert by construction; a VM re-run is the
  pending confirmation. NOTHING WAS COMMITTED: all changes are uncommitted on `lm-mhd-upgrades` (M: error_table.md + 10
  regenerated PNGs + xmhd_2d.F90 stub; ?? deliverables/, the two new analysis scripts, overlay_Ha{2,3}.png, this log).
- 2026-06-30T16:27:25Z  THE ONE QUESTION FOR YUCHEN (most unblocks the next step): should OFT aim MUG at the standalone
  full-induction blanket duct, or stay focused on the coupled plasma-blanket problem via TokaMaker (OFT's distinctive niche,
  since standalone duct codes have no mechanism for the evolving plasma background)? The answer selects the roadmap branch.
  Two concrete sub-questions follow: which LM normalization / material parameters (PbLi vs Li) does the ORNL side use; and
  is the conducting-wall benchmark they care about Hunt, Shercliff, or a duct with a specified wall-conductance ratio.

### Files created or modified this session

Created:
- `deliverables/lm_mhd_for_yuchen.zip`  (one-attachment package; listed first per the brief)
- `deliverables/lm_mhd_for_yuchen/README.md`
- `deliverables/lm_mhd_for_yuchen/OFT_LM_MHD_state_and_roadmap.md`
- `deliverables/lm_mhd_for_yuchen/SUMMARY_email.md`
- `deliverables/lm_mhd_for_yuchen/results/` (error_table.md, flowrate_vs_Ha.png, overlay_Ha{1,2,3,5,10,20,50}.png,
  profile_Ha{1,5,10,20,50}.csv, report/{threeway_overlay_Ha5,umean_uc_vs_Ha,error_vs_Ha,flowrate_vs_Ha}.png,
  mug_profiles/Ha_{2,3,5}_profile.csv, hunt/{hunt_sidejet_table.md, hunt_midplane_Ha*.png, hunt_contour_Ha*.png})
- `cases/hartmann/recompute_headline.py`  (integrity recompute of the headline from stored CSVs)
- `cases/hartmann/plots_from_csv.py`  (rebuild all plots from CSVs, no meshio/OpenFOAM dependency)
- `cases/hartmann/results/overlay_Ha2.png`, `cases/hartmann/results/overlay_Ha3.png`  (new MUG-only overlays)
- `docs/dev/lm_mhd_packaging_log.md`  (this log)

Modified:
- `src/physics/xmhd_2d.F90`  (+11 lines: inert default-off `c_wall` stub + doc comment; nothing else)
- `cases/hartmann/results/error_table.md`  (rebuilt: complete, honest blanks, MUG Ha={2,3,5}, integer+fit L2)
- `cases/hartmann/results/{flowrate_vs_Ha,overlay_Ha{1,5,10,20,50}}.png` and `results/report/*.png`  (regenerated from CSVs)
