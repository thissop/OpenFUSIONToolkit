#!/usr/bin/env bash
# Wall-function validation: does the Hartmann-drag closure reproduce the resolved
# duct core velocity on a COARSE (Hartmann-unresolved) mesh?
#
# Matrix per Ha (default Ha=100):
#   INSULATING (test_shercliff_mug), target Ha*u0 = 1
#     resolved      : ni=64 packing=4  drag OFF   -> reference core
#     coarse_drag   : ni=16 packing=1  duct_insulating -> should recover Ha*u0=1
#     coarse_channel: ni=16 packing=1  channel (CONTROL) -> wrong by ~Ha
#     coarse_nodrag : ni=16 packing=1  drag OFF (CONTROL) -> unresolved/wrong
#   CONDUCTING (test_hunt_mug, c=1),  target Ha^2*u0 = 1+1/c = 2
#     resolved      : ni=64 packing=4  c_wall=1 drag OFF -> reference core
#     coarse_drag   : ni=16 packing=1  duct_conducting drag_c=1 -> should recover Ha^2*u0=2
#
# Requires: OFT built at BUILD_DIR with test_shercliff_mug + test_hunt_mug
# (both carry the drag-on namelist switch; see DRAG_CLOSURE_NEXT.md).
#SBATCH --job-name=oft_wallfunc_val
#SBATCH --account=astro
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/%x-%j.out

set -euo pipefail

REPO_DIR="${REPO_DIR:?Set REPO_DIR to the staged OpenFUSIONToolkit repo}"
BUILD_DIR="${BUILD_DIR:-${REPO_DIR}/builds/build_release}"
HA="${HA:-100}"                 # sweep: submit with HA=20 / HA=100 / HA=1000
CWALL="${CWALL:-1.0}"           # conductance ratio c for the conducting case
OUT_ROOT="${REPO_DIR}/cases/shercliff/mug/wallfunc/Ha_${HA//./p}"
RESULT_DIR="${REPO_DIR}/cases/lm_mhd_papers/results"
VENV_DIR="${REPO_DIR}/.venv_lm_mhd_ginsburg_python"
SHERCLIFF_BIN="${BUILD_DIR}/tests/physics/test_shercliff_mug"
HUNT_BIN="${BUILD_DIR}/tests/physics/test_hunt_mug"
export HA CWALL OUT_ROOT

cd "${REPO_DIR}"
mkdir -p logs "${RESULT_DIR}" "${OUT_ROOT}"

if [ -f /etc/profile.d/modules.sh ]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
fi
module purge || true
module load gcc/11.2.0 openmpi/gcc/64/4.1.7a1 hdf5p/1.14.6 openblas/dynamic/0.3.7 lapack/gcc/64/3.9.0 || true
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

if [ ! -x "${VENV_DIR}/bin/python" ]; then
  /burg/opt/anaconda3-2023.09/bin/python3 -m venv "${VENV_DIR}"
fi
# shellcheck disable=SC1091
source "${VENV_DIR}/bin/activate"
python -m pip install --quiet numpy

for bin in "${SHERCLIFF_BIN}" "${HUNT_BIN}"; do
  if [ ! -x "${bin}" ]; then
    echo "wallfunc_status=FAIL"; echo "missing_executable=${bin}"; exit 2
  fi
done

# --- generate every oft.in in the matrix -----------------------------------
python - <<'PY'
from pathlib import Path
import math, os

root = Path(os.environ["OUT_ROOT"])
ha = float(os.environ["HA"]); cwall = float(os.environ["CWALL"])
a_half, eta, nu, fy, n0, den_scale = 1.0, 1.0, 1.0, 1.0e-2, 1.0, 5.98e26
mu0 = 4.0e-7 * math.pi; proton_mass = 1.67262192369e-27
rho = proton_mass * n0 * den_scale
b0 = ha * a_half * math.sqrt(mu0 * eta * rho * nu)

def mesh(ni, packing):
    return (f"&mesh_options\n meshname='test'\n cad_type=92\n nlevels=1\n nbase=1\n"
            f" grid_order=1\n fix_boundary=T\n/\n\n&cube_options\n mesh_type=1\n"
            f" ni={ni},{ni},0\n rscale=2.0,2.0,0.\n shift=-1.0,-1.0,0.\n"
            f" ref_per=F,F,F\n packing={packing},{packing},1.0\n/\n")

runtime = "&runtime_options\n ppn=1\n debug=0\n test_run=T\n use_petsc=F\n/\n\n"

def sher(ni, packing, drag):
    # drag: None | 'duct_insulating' | 'channel'
    d = " use_wall_drag=F\n" if drag is None else \
        f" use_wall_drag=T\n drag_closure='{drag}'\n drag_Ha={ha:.10e}\n"
    return runtime + mesh(ni, packing) + (
        f"&shercliff_options\n order=1\n nsteps=2000\n dt=2.5e-3\n ittarget=40\n"
        f" a_half={a_half}\n B0={b0:.12e}\n eta={eta}\n nu={nu}\n fy={fy}\n n0={n0}\n"
        f" t0=1.0\n chi=1.0e-30\n D_diff=1.0e-30\n gamma=1.67\n"
        f" den_scale={den_scale:.12e}\n lin_tol=1.0e-9\n nl_tol=1.0e-7\n pm=T\n"
        + d + "/\n")

def hunt(ni, packing, drag):
    d = " use_wall_drag=F\n" if drag is None else \
        f" use_wall_drag=T\n drag_closure='{drag}'\n drag_Ha={ha:.10e}\n drag_c={cwall:.10e}\n"
    return runtime + mesh(ni, packing) + (
        f"&hunt_options\n order=1\n nsteps=3000\n dt=2.5e-3\n ittarget=40\n"
        f" a_half={a_half}\n B0={b0:.12e}\n eta={eta}\n nu={nu}\n fy={fy}\n n0={n0}\n"
        f" t0=1.0\n chi=1.0e-30\n D_diff=1.0e-30\n gamma=1.67\n"
        f" den_scale={den_scale:.12e}\n lin_tol=1.0e-9\n nl_tol=1.0e-7\n"
        f" c_wall={cwall:.10e}\n pm=T\n" + d + "/\n")

cases = {
    "ins_resolved":       ("shercliff", sher(64, 4.0, None)),
    "ins_coarse_drag":    ("shercliff", sher(16, 1.0, "duct_insulating")),
    "ins_coarse_channel": ("shercliff", sher(16, 1.0, "channel")),
    "ins_coarse_nodrag":  ("shercliff", sher(16, 1.0, None)),
    "con_resolved":       ("hunt",      hunt(64, 4.0, None)),
    "con_coarse_drag":    ("hunt",      hunt(16, 1.0, "duct_conducting")),
}
for name, (kind, text) in cases.items():
    d = root / name; d.mkdir(parents=True, exist_ok=True)
    (d / "oft.in").write_text(text, encoding="utf-8")
    (d / "kind.txt").write_text(kind + "\n", encoding="utf-8")
print(f"B0={b0:.6e}  rho={rho:.6e}  wrote {len(cases)} cases under {root}")
PY

# --- run every case ---------------------------------------------------------
for d in "${OUT_ROOT}"/*/; do
  kind=$(cat "${d}/kind.txt")
  if [ "${kind}" = "shercliff" ]; then cp "${SHERCLIFF_BIN}" "${d}/run_bin"; else cp "${HUNT_BIN}" "${d}/run_bin"; fi
  ( cd "${d}" && ./run_bin ) > "${d}/run.log" 2>&1 || echo "RUN FAILED: ${d}" | tee -a "${d}/run.log"
done

# --- extract core velocity, compare to closure targets ----------------------
python - <<'PY' | tee "${RESULT_DIR}/wallfunc_validation_Ha_${HA//./p}.txt"
from pathlib import Path
import os, glob
root = Path(os.environ["OUT_ROOT"]); ha = float(os.environ["HA"]); cwall = float(os.environ["CWALL"])
nu, fy, a = 1.0, 1.0e-2, 1.0

def core_u0(profile):
    best, bu = 1e9, None
    for ln in open(profile):
        if ln.startswith("#"): continue
        p = ln.split()
        if len(p) < 3: continue
        x, z, v = map(float, p[:3])
        r = x*x + z*z
        if r < best: best, bu = r, v
    return bu * nu / (fy * a*a)          # nondim core u0

print(f"# Ha={ha}  c={cwall}  insulating target Ha*u0=1 ; conducting target Ha^2*u0={1+1/cwall:.3f}")
for d in sorted(root.glob("*/")):
    prof = glob.glob(str(d / "*_mug.profile"))
    if not prof: print(f"{d.name:20s}  NO PROFILE (run failed?)"); continue
    u0 = core_u0(prof[0])
    if d.name.startswith("ins"): metric, tgt = ha*u0, 1.0
    else:                        metric, tgt = ha*ha*u0, 1+1/cwall
    tag = "OK " if abs(metric-tgt)/tgt < 0.1 else "XX " if "coarse_drag" in d.name or "resolved" in d.name else "-- "
    print(f"{d.name:20s}  u0={u0:.6e}  metric={metric:.4f}  target={tgt:.4f}  {tag}")
PY

echo "wallfunc_status=DONE  results in ${RESULT_DIR}/wallfunc_validation_Ha_${HA//./p}.txt"
