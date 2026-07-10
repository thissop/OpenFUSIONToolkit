#!/usr/bin/env bash
#SBATCH --job-name=oft_mug_ha_probe
#SBATCH --account=astro
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/%x-%j.out

set -euo pipefail

REPO_DIR="${REPO_DIR:?Set REPO_DIR to the staged OpenFUSIONToolkit repo}"
BUILD_DIR="${BUILD_DIR:-${REPO_DIR}/builds/build_release}"
HA="${HA:-10}"
NY="${NY:-160}"
PACKING="${PACKING:-8.0}"
NSTEPS="${NSTEPS:-800}"
DT="${DT:-2.0e-2}"
FX="${FX:-1.0e-2}"
DEN_SCALE="${DEN_SCALE:-5.98e26}"
CASE_TAG="Ha_${HA}_ny_${NY}_pack_${PACKING//./p}_dt_${DT//./p}_nsteps_${NSTEPS}"
CASE_DIR="${REPO_DIR}/cases/hartmann/mug_ginsburg/${CASE_TAG}"
RESULT_DIR="${REPO_DIR}/cases/lm_mhd_papers/results"
VENV_DIR="${REPO_DIR}/.venv_lm_mhd_ginsburg_python"
EXECUTABLE="${BUILD_DIR}/tests/physics/test_hartmann_mug"
export CASE_DIR HA NY PACKING NSTEPS DT FX DEN_SCALE

cd "${REPO_DIR}"
mkdir -p logs "${RESULT_DIR}" "${CASE_DIR}"

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
python -m pip install --quiet pytest numpy scipy h5py matplotlib smplotlib

if [ ! -x "${EXECUTABLE}" ]; then
  echo "mug_probe_status=FAIL"
  echo "missing_executable=${EXECUTABLE}"
  exit 2
fi

cp "${EXECUTABLE}" "${CASE_DIR}/test_hartmann_mug"

python - <<'PY'
from pathlib import Path
import math
import os

case_dir = Path(os.environ["CASE_DIR"])
ha = float(os.environ["HA"])
ny = int(os.environ["NY"])
packing = float(os.environ["PACKING"])
nsteps = int(os.environ["NSTEPS"])
dt = os.environ["DT"]
fx = os.environ["FX"]
den_scale = float(os.environ["DEN_SCALE"])
a_half = 1.0
eta = 1.0
nu = 1.0
n0 = 1.0
mu0 = 4.0e-7 * math.pi
proton_mass = 1.67262192369e-27
rho = proton_mass * n0 * den_scale
b0 = ha * a_half * math.sqrt(mu0 * eta * rho * nu)

(case_dir / "oft.in").write_text(
    f"""&runtime_options
 ppn=1
 debug=0
 test_run=T
 use_petsc=F
/

&mesh_options
 meshname='test'
 cad_type=92
 nlevels=1
 nbase=1
 grid_order=1
 fix_boundary=T
/

&cube_options
 mesh_type=1
 ni=4,{ny},0
 rscale=0.5,2.0,0.
 shift=-0.25,-1.0,0.
 ref_per=T,F,F
 packing=1.0,{packing:.8g},1.0
/

&hartmann_options
 order=1
 nsteps={nsteps}
 dt={dt}
 ittarget=1000000
 a_half=1.0
 B0={b0:.12e}
 eta=1.0
 nu=1.0
 fx={fx}
 n0=1.0
 t0=1.0
 chi=1.0e-30
 D_diff=1.0e-30
 gamma=1.67
 den_scale={den_scale:.12e}
 lin_tol=1.0e-9
 nl_tol=1.0e-7
 pm=F
/
""",
    encoding="utf-8",
)

(case_dir / "probe_metadata.txt").write_text(
    "\n".join(
        [
            f"Ha={ha}",
            f"ny={ny}",
            f"packing={packing}",
            f"nsteps={nsteps}",
            f"dt={dt}",
            f"fx={fx}",
            f"den_scale={den_scale:.12e}",
            f"rho={rho:.12e}",
            f"B0={b0:.12e}",
            f"hartmann_layer_thickness={a_half / ha:.12e}",
            f"nominal_wall_cells_per_layer={ny / (2.0 * ha):.6f}",
        ]
    )
    + "\n",
    encoding="utf-8",
)
PY

{
  echo "timestamp_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "hostname=$(hostname)"
  echo "job_id=${SLURM_JOB_ID:-none}"
  echo "repo_dir=${REPO_DIR}"
  echo "build_dir=${BUILD_DIR}"
  echo "case_dir=${CASE_DIR}"
  echo "ha=${HA}"
  echo "ny=${NY}"
  echo "packing=${PACKING}"
  echo "nsteps=${NSTEPS}"
  echo "dt=${DT}"
  echo "executable=${EXECUTABLE}"
  echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
  echo "OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS}"
} | tee "${CASE_DIR}/probe_env.txt"

set +e
(
  cd "${CASE_DIR}"
  ./test_hartmann_mug
) > "${CASE_DIR}/run.log" 2>&1
run_status=$?
set -e
echo "EXIT=${run_status}" >> "${CASE_DIR}/run.log"

if [ "${run_status}" -ne 0 ]; then
  {
    echo "mug_probe_status=FAIL"
    echo "stage=run"
    echo "run_status=${run_status}"
    echo "case_dir=${CASE_DIR}"
  } | tee "${CASE_DIR}/probe_status.txt" "${RESULT_DIR}/ginsburg_mug_hartmann_probe_${CASE_TAG}_status.txt"
  exit "${run_status}"
fi

if [ ! -s "${CASE_DIR}/hartmann_mug.profile" ]; then
  {
    echo "mug_probe_status=FAIL"
    echo "stage=profile"
    echo "case_dir=${CASE_DIR}"
  } | tee "${CASE_DIR}/probe_status.txt" "${RESULT_DIR}/ginsburg_mug_hartmann_probe_${CASE_TAG}_status.txt"
  exit 3
fi

python cases/hartmann/analyze_mug.py "${CASE_DIR}/hartmann_mug.profile" 1.0 "MUG Ha=${HA}" \
  2>&1 | tee "${CASE_DIR}/analyze.log"

{
  echo "mug_probe_status=PASS"
  echo "case_dir=${CASE_DIR}"
  echo "timestamp_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} | tee "${CASE_DIR}/probe_status.txt" "${RESULT_DIR}/ginsburg_mug_hartmann_probe_${CASE_TAG}_status.txt"
