# LM-MHD paper execution log

This log records the execution pass started on 2026-07-09 for the two-paper
LM-MHD roadmap:

- Paper 1: rigorous liquid-metal MHD verification capability.
- Paper 2: plasma/liquid-metal coupling with a nonlinear free-surface target
  and explicit thermocapillary / linear-surface fallbacks.

The log is intentionally conservative: measured results, job IDs, and blockers
are recorded as they occur.

## 2026-07-09

- 2026-07-09T14:16:52Z - Execution started from Codex Default mode after plan
  lock-in. Local branch `lm-mhd-upgrades`, HEAD `357fcd7`, dirty tree observed
  with existing Hartmann/Hunt/LM-MHD artifacts and modified Python/Fortran files.
  Ginsburg SSH ControlMaster verified: host `ruth`, user `tjk2147`, remote cwd
  `/burg-archive/home/tjk2147`, one existing user job in queue. No cluster jobs
  submitted yet.
- 2026-07-09T14:16:56Z - Ginsburg environment probe complete. Available stack
  includes `cmake/3.28.1`, GCC modules `10.2.0/11.2.0/13.0.1`, OpenMPI
  `4.1.5a1/4.1.7a1`, PETSc `3.21.2`, and Python 3.11.8. Selected staging root
  `/burg/astro/users/tjk2147/spf_work/oft_lm_mhd_papers`, using Slurm account
  `astro` unless rejected.
- 2026-07-09T14:20:00Z - Staged current dirty worktree to Ginsburg at
  `/burg/astro/users/tjk2147/spf_work/oft_lm_mhd_papers/repo`. Submitted
  Python/reference baseline job `8901303` on `short` and OFT build/regression
  baseline job `8901304` on `short`. Both initially entered `RUNNING`; no
  results collected yet.
- 2026-07-09T14:25:43Z - First baseline submission finished with harness
  failures, not physics failures. Python/reference job `8901303` failed because
  `smplotlib` was missing in the job environment. OFT build job `8901304`
  failed before a valid configure/build because module `cmake/3.28.1` hit a
  `GLIBCXX` runtime mismatch and the fallback branch did not enter a valid build
  directory. Patched the Ginsburg job scripts to install `smplotlib`, use the
  Anaconda CMake binary, attempt a source-dependency fallback, and write staging
  metadata from local branch/head before excluding `.git`.
- 2026-07-09T14:28:13Z - Restaged the current dirty worktree after the harness
  fixes. Submitted Python/reference baseline job `8901355` on `short` and OFT
  build/regression baseline job `8901356` on `short`; both entered `RUNNING`
  (`g036` and `g180`, respectively).
- 2026-07-09T14:29:18Z - Build/regression job `8901356` failed in 5 seconds
  because the two concurrently submitted jobs shared a virtualenv path and the
  build job tried to activate it while the Python job was still creating it.
  Patched the scripts to use separate job virtualenvs:
  `.venv_lm_mhd_ginsburg_python` and `.venv_lm_mhd_ginsburg_build`. Python job
  `8901355` remained running.
- 2026-07-09T14:30:31Z - Python/reference baseline job `8901355` completed
  successfully on Ginsburg: focused LM-MHD Python tests passed (`53 passed in
  1.86s`), the Hartmann analytic self-test passed, and the LM-MHD artifact
  manifest refreshed. Submitted replacement OFT build/regression job `8901375`
  after the virtualenv race fix; it entered `RUNNING` on `g180`.
- 2026-07-09T14:33:28Z - OFT build/regression job `8901375` reached direct
  CMake configure. The module-based configure failed because `OFT_METIS` was
  not found in the loaded Ginsburg module stack; the script correctly entered
  the source-dependency fallback path.
- 2026-07-09T14:39:28Z - Added a first Paper 2 plasma/liquid-metal interface
  utility layer: Bohm/sheath particle, heat, current, mass, and momentum load
  accounting; evaporative heat/mass balance helpers; and linear
  capillary-gravity surface response. Added tests for the identities and input
  validation. Expanded local focused LM-MHD Python suite passed: `59 passed in
  1.33s`.
- 2026-07-09T14:44:21Z - Added
  `docs/dev/lm_mhd_two_paper_execution_plan.md` with Paper 1/Paper 2 gates.
  Added and ran `cases/lm_mhd_coupling/plasma_interface_sheath_scan.py`,
  producing CSV/PNG/Markdown reduced sheath-load artifacts. Updated the artifact
  manifest to include `plasma_interface_sheath_scan`; it now reports 35 files,
  10 CSV files, and 15 PNG files with 15/15 nonblank. Re-ran the expanded local
  focused suite: `59 passed in 0.76s`; `git diff --check` was clean for touched
  files.
- 2026-07-09T14:45:52Z - Ginsburg OFT build/regression job `8901375` failed
  in source-dependency fallback after SuperLU compiled and installed
  successfully into `lib64`; `build_libs.py` expected `lib/libsuperlu.a`.
  Patched the Ginsburg build wrapper to add a `lib -> lib64` compatibility
  symlink and retry `build_libs.py` when this specific SuperLU layout appears.
- 2026-07-09T14:47:26Z - Updated the staged Ginsburg repo in place with the new
  plasma-interface module, tests, sheath-scan script, updated manifest script,
  docs, and patched Slurm wrappers. Submitted expanded Python/reference baseline
  job `8902072` and patched OFT build/regression job `8902073`.
- 2026-07-09T14:49:32Z - Ginsburg expanded Python/reference baseline job
  `8902072` completed successfully: focused LM-MHD Python tests passed (`59
  passed in 2.08s`), Hartmann analytic self-test passed, and the artifact
  manifest refreshed. Build/regression job `8902073` passed the SuperLU `lib64`
  workaround and reached active OFT compilation.
- 2026-07-09T15:38:32Z - Lost non-interactive Ginsburg access: the SSH
  ControlMaster socket no longer exists and `BatchMode` SSH now fails before
  polling. Last observed state for build/regression job `8902073`: `RUNNING` on
  `g099` at elapsed `00:13:42`; OFT had built to 100%, and the regression
  pytest log showed progress (`..ssss..ssss`) but no final pass/fail yet. Small
  artifact job `8902095` completed before socket loss, generating the
  plasma-interface sheath scan and refreshing the remote manifest to 35 files.
  Final build/regression status is unobserved until SSH is reactivated.
- 2026-07-09T18:04:40Z - User asked to keep going. Retried `BatchMode` SSH to
  Ginsburg; the host was reachable but authentication failed without an active
  ControlMaster/Duo session, so remote polling remains blocked. Continued
  locally on the Paper 2 L1 thermal-response rung.
- 2026-07-09T18:07:43Z - Added plasma-interface thermal response references:
  representative LM thermal properties, thermal diffusivity/penetration depth,
  lumped surface-layer energy balance, semi-infinite constant-heat-flux analytic
  surface temperature, and a thermal-regime classifier. Added tests for energy
  closure, sqrt(time) conduction scaling, and regime labels. Added
  `cases/lm_mhd_coupling/plasma_interface_thermal_response_scan.py` and
  regenerated plasma-interface artifacts. Expanded focused local LM-MHD Python
  suite passed: `62 passed in 1.00s`. Manifest now reports 38 files, 11 CSV
  files, and 16 PNG files with 16/16 nonblank.
- 2026-07-09T18:12:43Z - Added Paper 2 L2 reduced biased-interface
  current-feedback rung: `surface_current_ohmic_heat_flux` with tests, plus
  `cases/lm_mhd_coupling/plasma_interface_current_feedback_scan.py`. The scan
  treats floating sheath as zero-current by default and explores an explicit
  biased/electrode-style current fraction, estimating ohmic heat and
  sheet-current magnetic feedback. Expanded focused local LM-MHD Python suite
  passed: `63 passed in 0.68s`. Manifest now reports 41 files, 12 CSV files,
  and 17 PNG files with 17/17 nonblank.
- 2026-07-09T18:16:52Z - Added Paper 2 L3 linear thermocapillary fallback
  rung: Marangoni shear stress, thermocapillary Couette profile, and mean
  velocity helpers with sign/linearity tests, plus
  `cases/lm_mhd_coupling/plasma_interface_marangoni_scan.py`. Expanded focused
  local LM-MHD Python suite passed: `65 passed in 0.66s`. Manifest now reports
  44 files, 13 CSV files, and 18 PNG files with 18/18 nonblank.
- 2026-07-10T00:35:58Z - User reactivated SSH multiplexing. Recovered final
  Ginsburg state for OFT build/regression job `8902073`: `COMPLETED 0:0` on
  `g099`. OFT built successfully via source-dependency fallback; targeted
  physics regressions passed (`21 passed, 36 skipped`, 24 pytest mark warnings)
  in `2079.55s`. Status file reports `oft_build_status=PASS` and
  `build_dir=/burg/astro/users/tjk2147/spf_work/oft_lm_mhd_papers/repo/builds/build_release`.
- 2026-07-10T00:37:12Z - Synced current L0-L3 plasma-interface work, updated
  Slurm wrappers, and docs to the staged Ginsburg repo. Submitted expanded
  Python/reference/artifact baseline job `8910539` on `short`.
- 2026-07-10T00:40:12Z - Ginsburg job `8910539` completed successfully on
  `g256`: focused LM-MHD Python tests passed (`65 passed in 38.88s`), Hartmann
  analytic self-test passed, all four L0-L3 plasma-interface artifact scripts
  ran, and the remote manifest refreshed to 44 files.
- 2026-07-10T00:42:33Z - Added Paper 2 L4 linear free-surface tractability
  gate: magnetic pressure-stiffness scale, low-Rm magnetic damping rate,
  capillary-gravity-magnetic wave frequency/period, conservative explicit
  time-step gate, and combined free-surface metrics. Added
  `cases/lm_mhd_coupling/plasma_interface_free_surface_gate_scan.py`.
  Expanded focused local LM-MHD Python suite passed: `67 passed in 0.70s`.
  Manifest now reports 47 files, 14 CSV files, and 19 PNG files with 19/19
  nonblank.
- 2026-07-10T00:44:32Z - Synced the L4 free-surface gate code, tests, scan,
  manifest wrapper, and docs to the staged Ginsburg repo. Submitted expanded
  Python/reference/artifact baseline job `8910555` on `short`; it was pending
  on priority at submission.
- 2026-07-10T00:46:40Z - Ginsburg job `8910555` completed successfully on
  `g178`: focused LM-MHD Python tests passed (`67 passed in 2.08s`), Hartmann
  analytic self-test passed, all five plasma-interface artifact scripts ran,
  and the remote manifest refreshed to 47 files. Local `git diff --check` for
  touched files and `bash -n tools/lm_mhd_papers/ginsburg_python_baseline.sbatch`
  were clean.
- 2026-07-10T00:48:53Z - Added
  `cases/lm_mhd_coupling/run_reference_matrix.py`, a reproducibility driver for
  the five plasma-interface artifact scripts. Local driver run passed with all
  five commands returning zero; local focused suite recheck passed (`67 passed
  in 0.86s`); manifest now reports 49 files, 15 CSV files, and 19 PNG files
  with 19/19 nonblank. Updated the Ginsburg Python baseline wrapper to call the
  driver before refreshing the artifact manifest.
- 2026-07-10T00:50:05Z - Synced the reference-matrix driver, updated manifest,
  Ginsburg Python wrapper, and docs to the staged Ginsburg repo. Submitted
  Python/reference/artifact baseline job `8910567` on `short`; it was pending
  on priority at submission.
- 2026-07-10T00:51:15Z - Ginsburg job `8910567` completed successfully on
  `g174`: focused LM-MHD Python tests passed (`67 passed in 2.40s`), Hartmann
  analytic self-test passed, all five reference-matrix commands returned zero,
  and the remote manifest refreshed to 49 files.
- 2026-07-10T00:53:06Z - Added
  `tools/lm_mhd_papers/ginsburg_mug_hartmann_probe.sbatch`, a Ginsburg wrapper
  for resolved-MUG Hartmann probe cases. Submitted Ha=10, `NY=160`,
  `PACKING=8.0`, `NSTEPS=800`, drag-off probe job `8910574` on `short`;
  pending on priority at submission.
- 2026-07-10T00:53:56Z - Probe job `8910574` failed immediately with a wrapper
  error (`KeyError: CASE_DIR`) before running MUG. Patched the wrapper to export
  generated case parameters to the Python input writer. Resubmitted the same
  Ha=10 packed-grid probe as job `8910576`; pending on priority at submission.
- 2026-07-10T00:55:50Z - Ha=10 probe job `8910576` was running on `g174` and
  consuming CPU but had not yet emitted timestep data. Patched and synced the
  probe wrapper for future runs to set `OMP_NUM_THREADS` and
  `OPENBLAS_NUM_THREADS` from `SLURM_CPUS_PER_TASK`; current job left running.
- 2026-07-10T00:57:12Z - Ha=10 probe job `8910576` still had no completed
  timestep after several polls. Patched the probe wrapper to include `dt` and
  `nsteps` in the case tag to avoid overwriting concurrent probes. Submitted
  fallback Ha=10, `NY=160`, `PACKING=8.0`, `NSTEPS=400`, `DT=5.0e-3` job
  `8910589`; pending on priority at submission.
- 2026-07-10T00:58:20Z - Fallback Ha=10 job `8910589` reached normal timestep
  output (>200 steps observed) with `DT=5.0e-3`. Original `DT=2.0e-2` job
  `8910576` remained stuck before the first timestep and was cancelled to avoid
  wasting allocation; treat it as a timestep-stiffness/probe failure, not a
  physics comparison.
- 2026-07-10T00:59:20Z - Fallback Ha=10 packed-grid MUG probe job `8910589`
  completed successfully on `g209` (`0:0`, elapsed `00:01:16`). Run produced
  `hartmann_mug.profile` with 1445 points and `analyze_mug.py` reported:
  developed-flow nonuniformity `1.347e-7`, `Ha_pred=10`, `Ha_fit=10.01`,
  `fit/pred=1.001`, normalized-profile L2 `3.406e-4` and Linf `5.215e-4`
  vs fitted-Ha analytic, and `u_mean/u_c=0.9004` vs analytic `0.9002`.
  Pulled the case artifacts locally under `cases/hartmann/mug_ginsburg/` and
  regenerated `cases/hartmann/results/error_table.md` from stored CSVs; Ha=10
  now has both mhdFoam and MUG entries. Ha=20/50 MUG remain open.
- 2026-07-10T01:43:20Z - Continued Paper 1 high-Ha resolved-MUG extension.
  Submitted drag-off Ha=20 packed-grid Ginsburg probe job `8910691` with
  `NY=320`, `PACKING=8.0`, `DT=1.25e-3`, and `NSTEPS=1600` (final time 2.0);
  pending on priority at submission.
- 2026-07-10T01:43:56Z - Ha=20 probe job `8910691` started on `g164`.
  Input metadata confirmed `B0=2.242251959124e-02`, Hartmann-layer thickness
  `5.0e-2`, nominal 8 wall-normal cells/layer, and thread limits
  `OMP_NUM_THREADS=8`, `OPENBLAS_NUM_THREADS=8`.
- 2026-07-10T01:44:20Z - Ha=20 probe job `8910691` reached normal timestep
  output at `DT=1.25e-3`; mesh reported `hmin=1.151e-03` and early timesteps
  were stable. Result still pending until profile extraction/analyzer complete.
- 2026-07-10T01:55:28Z - Ha=20 packed-grid MUG probe job `8910691`
  completed successfully on `g164` (`0:0`, elapsed `00:09:19`). Run produced
  `hartmann_mug.profile` with 2885 points and `analyze_mug.py` reported:
  developed-flow nonuniformity `1.312e-3`, `Ha_pred=20`, `Ha_fit=20.01`,
  `fit/pred=1.0004`, normalized-profile L2 `5.548e-5` and Linf `2.697e-4`
  vs fitted-Ha analytic, and `u_mean/u_c=0.9500` vs analytic `0.9500`.
  Pulled the case artifacts locally under `cases/hartmann/mug_ginsburg/` and
  regenerated `cases/hartmann/results/error_table.md`; Ha=20 now has both
  mhdFoam and MUG entries. Ha=20 remains labeled a probe until convergence is
  checked; Ha=50 MUG remains open.
- 2026-07-10T02:10:31Z - Ha=50 bounded MUG feasibility attempts did not
  produce a profile. Job `8910721` (`NY=500`, `packing=12.0`, `DT=2e-4`,
  `NSTEPS=5000`) failed before timestepping with OFT run status 99:
  `Ground point repeated on same domain` in `mesh_global_igrnd` during square
  mesh generation. Retry job `8910734` (`NY=500`, `packing=8.0`, same `DT` and
  `NSTEPS`) was cancelled after 00:08:42 because it had metadata only, no mesh
  statistics or timestep output, and low CPU accumulation. Both failed attempt
  directories were pulled locally. Ha=50 remains blocked/open for this probe
  family; no Ha=50 MUG profile comparison was produced.
- 2026-07-10T02:27:08Z - Ha=20 `NY=400`, `packing=8.0` mesh-refinement
  repeat did not reach a profile. Job `8910763` on `g063` was cancelled after
  00:02:23 because it had metadata only, no mesh/timestep output, and low CPU
  accumulation. Retry job `8910766` with `--exclude=g063` landed on `g090` and
  was cancelled after 00:11:56; it was CPU-active but still emitted only the
  initial branch/revision lines, with no OFT initialization, mesh statistics,
  timestep output, profile, or status file. The attempt directory was pulled
  locally. Treat this as a setup/pre-initialization blocker for the `NY=400`
  refinement wrapper, not a profile comparison.
- 2026-07-10T02:34:03Z - Ha=20 timestep-refinement repeat also did not reach
  a profile. Job `8911505` (`NY=320`, `packing=8.0`, `DT=6.25e-4`,
  `NSTEPS=3200`, `--exclude=g063,g090`) landed on `g198` and was cancelled
  after 00:05:44. It wrote metadata but only the initial branch/revision
  run.log lines, with no OFT initialization, mesh statistics, timestep output,
  profile, or status file; CPU accumulation stayed low. The attempt directory
  was pulled locally. Ha=20 convergence evidence remains open beyond the
  successful single `NY=320`, `DT=1.25e-3` probe.
