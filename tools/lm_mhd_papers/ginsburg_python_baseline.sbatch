#!/usr/bin/env bash
#SBATCH --job-name=oft_lm_py_base
#SBATCH --account=astro
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%j.out

set -euo pipefail

REPO_DIR="${REPO_DIR:?Set REPO_DIR to the staged OpenFUSIONToolkit repo}"
cd "${REPO_DIR}"

mkdir -p logs cases/lm_mhd_papers/results

RESULT_DIR="${REPO_DIR}/cases/lm_mhd_papers/results"
ENV_LOG="${RESULT_DIR}/ginsburg_python_baseline_env.txt"
TEST_LOG="${RESULT_DIR}/ginsburg_python_baseline_pytest.log"
HARTMANN_LOG="${RESULT_DIR}/ginsburg_hartmann_analytic_selftest.log"
VENV_DIR="${REPO_DIR}/.venv_lm_mhd_ginsburg_python"

if [ ! -x "${VENV_DIR}/bin/python" ]; then
  /burg/opt/anaconda3-2023.09/bin/python3 -m venv "${VENV_DIR}"
fi
# shellcheck disable=SC1091
source "${VENV_DIR}/bin/activate"
python -m pip install --quiet --upgrade pip
python -m pip install --quiet pytest numpy scipy h5py matplotlib smplotlib

{
  echo "timestamp_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "hostname=$(hostname)"
  echo "job_id=${SLURM_JOB_ID:-none}"
  echo "repo_dir=${REPO_DIR}"
  echo "branch=$(git branch --show-current 2>/dev/null || true)"
  echo "head=$(git rev-parse --short HEAD 2>/dev/null || true)"
  echo "python=$(command -v python || true)"
  python --version || true
} | tee "${ENV_LOG}"

export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="${REPO_DIR}/src/python"

python -m pytest -q \
  src/tests/physics/test_lm_mhd_coupling.py \
  src/tests/physics/test_lm_mhd_closures.py \
  src/tests/physics/test_lm_mhd_feedback.py \
  src/tests/physics/test_lm_mhd_inductionless.py \
  src/tests/physics/test_lm_mhd_full_induction_bridge.py \
  src/tests/physics/test_lm_mhd_plasma_interface.py \
  2>&1 | tee "${TEST_LOG}"

python cases/hartmann/analytic_hartmann.py 2>&1 | tee "${HARTMANN_LOG}"

python cases/lm_mhd_coupling/run_reference_matrix.py \
  2>&1 | tee "${RESULT_DIR}/ginsburg_reference_matrix.log"

python cases/lm_mhd_coupling/lm_mhd_artifact_manifest.py \
  2>&1 | tee "${RESULT_DIR}/ginsburg_artifact_manifest_refresh.log"

echo "python_baseline_status=PASS" | tee "${RESULT_DIR}/ginsburg_python_baseline_status.txt"
