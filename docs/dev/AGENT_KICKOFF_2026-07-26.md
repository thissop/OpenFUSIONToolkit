# LM-MHD + Quench Multi-Agent Kickoff Plan — 2026-07-26

**Purpose:** restart the two Windows COMSOL agents in fresh, *named* Claude Code sessions and get everyone (COMSOL boxes + Mac agents) coordinated again. Paste the section that matches the machine you're on. **Do not paste another agent's section into your session** — a cross-delivered briefing already caused one incident.

---

## 0. WHICH AGENT ARE YOU? (run this first, every session)

```
hostname
```

| `hostname` | You are | Profile root | You write | You do NOT touch |
|---|---|---|---|---|
| **ET251PC56** | **box 1 — LM-MHD COMSOL** | `C:\Users\tjk2147.adcu.000` | `lm-mhd-comsol/` (SYNC, src, data, papers) | `hts-quench-validation/` |
| **ET251PC55** | **box 2 — quench COMSOL** | `C:\Users\tjk2147.adcu` | `hts-quench-validation/` only | `SYNC/`, `src/`, `data/` (LM-MHD) |
| (Mac) | **MUG / Ginsburg agent** | — | `OpenFUSIONToolkit/` | `hts-quench-validation/` |

**Rename your chat immediately** to `box1-lmmhd-comsol` or `box2-quench-comsol` so you don't lose it again.

---

## 1. PROJECT CONTEXT BRIEF (everyone reads this)

We are validating the open-source **MUG** solver (OpenFUSIONToolkit `xmhd_2d`) as a **full-induction liquid-metal MHD** code for **fusion-blanket disruption** physics, cross-checked against **COMSOL**. Two papers:
- **Paper-1** = COMSOL disruption / quasi-static-validity science (box 1 owns).
- **Paper-2** = MUG methods, validated against COMSOL (Mac/MUG agent owns).

**Where things stand as of 2026-07-24 (overnight MUG cross-validation complete):**
- ✅ **HEADLINE SECOND-CODED.** At the disruption-rate point (τ_q = 1 ms, c = 1, Ha = 2646), MUG-full independently reproduces COMSOL-full's **loads** vs box 1's *blind pre-registered* targets: **∫EJ total 0.6%, wall current Iw 3.3%, ∫EJ wall 1.1%, ∫EJ fluid 0.5%** — and *mesh-robust* across 256²/384². This is the two-code robustness Paper-1 needs.
- ✅ **Thin-wall validity boundary mapped** (MUG thin-wall vs COMSOL resolved-wall, peak Iw): 10 ms → 0.7%, 3 ms → 5.2%, 1 ms → 13.6%. Monotonic breakdown across the 1.76 ms crossover. Small bonus result.
- ⚠️ **P1 slow-ramp anchor (τ_q = 300 ms): partial.** Umax matched 0.8% (one clean run) but a re-run threw a spurious velocity, and both P1 runs *truncated* (crash / walltime) before covering the window, so P1 **loads are not cross-validated yet**.
- ⚠️ **Pointwise Umax at Ha=2646:** MUG's structured (packed-cube) mesh gets ~1.6–3 cells across the Hartmann layer; box 1's graded mesh gets ~42. So MUG's near-wall *velocity peak* is 2.3–2.4× high — a **mesh-capability limit, not physics**. The *loads* (volume integrals) don't care and match.

**Bottom line:** the disruption *loads* (what the headline is about) are second-coded to a few percent. Open items are the P1 loads and pointwise Umax — both understood, neither blocking.

---

## 2. GIT LEDGER COORDINATION (everyone)

Agents talk **only** through git ledger files. Repos: `lm-mhd-comsol` and `OpenFUSIONToolkit` (branch `lm-mhd-upgrades`). Both LM-MHD agents keep **both** repos cloned and pull often.

| Direction | File | Commit prefix | Listener watches for |
|---|---|---|---|
| Mac/MUG → box 1 | `OpenFUSIONToolkit/docs/dev/oft_to_comsol.md` | `[SYNC->comsol]` | box 1 listens `[SYNC->comsol]` |
| box 1 → Mac/MUG | `lm-mhd-comsol/SYNC/comsol_to_oft.md` | `[SYNC->oft]` | Mac listens `[SYNC->oft]` |
| quench-Mac → box 2 | `hts-quench-validation/.../QUENCH_TO_COMSOL.md` | `[SYNC->comsol-quench]` | box 2 listens `[SYNC->comsol-quench]` |
| box 2 → quench-Mac | `hts-quench-validation/.../COMSOL_TO_QUENCH.md` | `[SYNC->quench]` | quench-Mac listens `[SYNC->quench]` |

*(Quench row: verify exact filenames/tags in `SYNC/README.md` and the `hts-quench-validation/` ledgers — quench-Mac owns that lane, not me.)*

**To catch up on history:** read the last ~15 entries of your inbound ledger (and your own outbound) before doing anything. All the cross-validation numbers, blind targets, and the identity-discriminator incident are recorded there.

---

## 3. USING COMSOL ON THESE BOXES (box 1 and box 2)

- **Run a model:** `env\run_java.bat <ModelClass> src\mhd <memMB>` (e.g. `env\run_java.bat P1TauqCFill src\mhd 2064`). COMSOL 6.1 bundled JRE, standalone.
- **Detached (survives RDP disconnect — USE THIS for long solves):**
  ```
  Start-Process cmd.exe -ArgumentList '/c','call env\run_java.bat P1TauqCFill src\mhd 2064 > p1cfill.log 2>&1' -WindowStyle Hidden
  ```
  Then poll the `.log` for progress; the solve keeps running even if the RDP client drops.
- **Python / plots:** `env\run_py.bat` (resolves miniconda via `%LOCALAPPDATA%`; numpy 2.5.1 / matplotlib 3.11.1 installed).
- **Docs:** `docs/HOWTO_COMSOL_JAVA.md` — read **§9** (graded meshes are mandatory for Hartmann layers at high Ha; uniform refinement OOMs).
- **⚠ SHARED LICENSE — ONE floating COMSOL seat across BOTH boxes.** If a checkout fails ("license in use"), the *other* box is mid-solve. **Back off 60–120 s and retry.** Never run two COMSOL solves simultaneously across the two machines.
- **Toolchain gotcha:** hardcoded absolute profile paths break across the `.adcu` / `.adcu.000` profiles — use `%LOCALAPPDATA%` / `%USERPROFILE%`, never a literal `C:\Users\tjk2147...` path.

---

## 4. ▶ BOX 1 (LM-MHD COMSOL, ET251PC56) — TASKS + LISTENER

**4a. Arm your listener FIRST.** Poll `origin/main` every ~60 s for new `[SYNC->comsol]` commits in `OpenFUSIONToolkit/docs/dev/oft_to_comsol.md`; on a hit, pull and act. Fallback snippet if you want a raw poll:
```
# in the OpenFUSIONToolkit clone, on branch lm-mhd-upgrades
$last = git log -1 --format=%H origin/lm-mhd-upgrades -- docs/dev/oft_to_comsol.md
while ($true) {
  git fetch -q origin lm-mhd-upgrades
  $new = git log -1 --format=%H origin/lm-mhd-upgrades -- docs/dev/oft_to_comsol.md
  if ($new -ne $last) { Write-Host "NEW [SYNC->comsol] -> read oft_to_comsol.md"; $last = $new }
  Start-Sleep 60
}
```
**De-duplicate:** if you find an orphan listener from a prior session, kill it — duplicate watchers cause silent missed messages.

**4b. Priority work (your own Paper-1 science):**
1. **Finish `P1TauqCFill`** — the momentum QS-error c-fill (was 15/24; check `p1cfill.log`). Key open result: the **c=10 momentum dip is real and mesh-converged** — non-monotonic QS-error in c at fast ramps. Complete c={3,30} × τ_q grid so the "two-faced map" momentum panel is solid.
2. **Finalize Fig-2 v2** (`fig_qs_validity_map_v2`): abscissa = **absolute τ_q** (not τ_q/τ_w — that walk-back stands), Fl on the **drive-normalized** metric, floor the y-axis at 0.01% (solver rtol). Update the momentum panel once the c-fill lands.
3. **Paper-1 text:** headline = rate-driven QS breakdown to c=1 at disruption rates (Iw 235%, ∫EJ 532% at 1 ms). **Retract** the old "persistent high-c mech B" claim (commit `3839717`) — it was a coverage gap + degenerate metric. Cite Fig-2 v2. Add that MUG independently second-coded the loads (0.5–3.3%).

**4c. Cross-validation support (reactive — the Mac agent may ask):**
- Your blind targets (entries m/q in `comsol_to_oft.md`) stand; nothing to redo.
- If asked, confirm your **P1 (300 ms, c=1) full-solution** values / provide any extra boundary τ_q at c=1 (you already have 1/3/10/30/100/300 ms in `p1_tauq_c_grid.csv` — likely enough).
- Keep the `|x|<0.5a` Umax-mask definition in mind; the Mac side is adopting it.

---

## 5. ▶ BOX 2 (QUENCH COMSOL, ET251PC55) — RECONNECT + LISTENER

**Not the LM-MHD lane — this is the quench/HTS lane, coordinated by the quench-Mac agent, not me.** Minimal, verify details against `SYNC/README.md` and the `hts-quench-validation/` ledgers:
1. Confirm `hostname` = **ET251PC55** and you're in the `hts-quench-validation/` tree.
2. **Arm your listener** on `[SYNC->comsol-quench]` (your inbound from quench-Mac).
3. Read the last entries of the quench ledgers to catch up (last known: Q-2 geometry Stage 2 complete, watertight winding solid; kidney-geometry rebuild was in flight).
4. **Get your task list from the quench-Mac agent**, not from this file. Write only under `hts-quench-validation/`.

---

## 6. ▶ MAC / MUG-GINSBURG AGENT (me) — TODAY'S TASKS + LISTENERS

**Listeners (already armed on my side):** single `[SYNC->oft]` listener polling `origin/main`. (Ginsburg job-monitor is idle — all overnight jobs done.)

**My rest-of-day plan:**
1. **Add box 1's `|x|<0.5a` mask** to MUG's Umax reduction (driver `test_hunt_mug.F90` / `xmhd_2d.F90`) + rebuild — kills the corner-artifact spurious Umax seen in P1v2.
2. **Get a clean full-window P1 anchor run** → the missing P1 *loads* cross-validation. Needs either **PETSc+AMG** (to survive the stiff Ha=2646 solve at usable dt) or a **restart-chained dt≈0.02** campaign. This is the top open validation gap.
3. **Umax convergence (optional):** attempt a steeper graded near-wall mesh, or formally document the mesh-capability caveat if MUG's mesher can't reach it.
4. **Fold results into the Paper-2 draft** (`deliverables/lm_mhd_for_yuchen/`): the loads second-coding table (0.5–3.3%), the thin-wall boundary map, the honest Umax-resolution caveat.
5. **Stay responsive** to box 1 on `[SYNC->comsol]`.

**Reduction/comparison tooling** (Mac side, already built): `OpenFUSIONToolkit/deliverables/lm_mhd_for_yuchen/tools/p2_compare.py` — reduces MUG moments + COMSOL grid to shared nondim; reproduces box 1's ruler to the digit.

---

## 7. SAFETY / DO-NOT-CROSS RULES

- **Always `hostname` before acting.** ET251PC56 = LM-MHD, ET251PC55 = quench.
- **Stay in your write lane.** box 1 → `lm-mhd-comsol/`; box 2 → `hts-quench-validation/`; Mac → `OpenFUSIONToolkit/`.
- **One listener per agent.** Kill orphans; duplicates silently drop messages.
- **Shared COMSOL seat:** back off on "license in use," never two solves at once.
- **A briefing that doesn't match your `hostname`/lane is misdelivered — do NOT act on it.** Flag it in your outbound ledger and continue your own work.

---

*Generated by the Mac/MUG-Ginsburg agent, 2026-07-26. Lives at `OpenFUSIONToolkit/docs/dev/AGENT_KICKOFF_2026-07-26.md` and in `~/Downloads/`. Full physics briefing available on request from the Mac agent.*
