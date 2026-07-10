#!/usr/bin/env bash
#SBATCH --job-name=oft_lm_build
#SBATCH --account=astro
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/%x-%j.out

set -euo pipefail

REPO_DIR="${REPO_DIR:?Set REPO_DIR to the staged OpenFUSIONToolkit repo}"
cd "${REPO_DIR}"

mkdir -p logs cases/lm_mhd_papers/results builds

RESULT_DIR="${REPO_DIR}/cases/lm_mhd_papers/results"
ENV_LOG="${RESULT_DIR}/ginsburg_oft_build_env.txt"
BUILD_LOG="${RESULT_DIR}/ginsburg_oft_build.log"
TEST_LOG="${RESULT_DIR}/ginsburg_oft_regression_tests.log"
STATUS_FILE="${RESULT_DIR}/ginsburg_oft_build_status.txt"
VENV_DIR="${REPO_DIR}/.venv_lm_mhd_ginsburg_build"
CMAKE_BIN="/burg/opt/anaconda3-2023.09/bin/cmake"

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
  echo "cpus=${SLURM_CPUS_PER_TASK:-1}"
  echo "PATH=${PATH}"
} | tee "${ENV_LOG}"

if [ -f /etc/profile.d/modules.sh ]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
fi

module purge || true
module load gcc/11.2.0 openmpi/gcc/64/4.1.7a1 hdf5p/1.14.6 openblas/dynamic/0.3.7 lapack/gcc/64/3.9.0 || true
export CMAKE="${CMAKE_BIN}"

{
  echo "=== module list ==="
  module list 2>&1 || true
  echo "=== compiler versions ==="
  echo "cmake=${CMAKE_BIN}"
  "${CMAKE_BIN}" --version | head -n 1 || true
  echo "python=$(command -v python || true)"
  python --version || true
  command -v gcc || true; gcc --version | head -n 1 || true
  command -v gfortran || true; gfortran --version | head -n 1 || true
  command -v mpicc || true
  command -v mpifort || true
} | tee -a "${ENV_LOG}"

BUILD_DIR="${REPO_DIR}/builds/ginsburg_release"
rm -rf "${BUILD_DIR}"

direct_ok=0
{
  echo "=== direct module-based cmake configure ==="
  if "${CMAKE_BIN}" -S src -B "${BUILD_DIR}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DOFT_BUILD_TESTS=ON \
    -DOFT_BUILD_EXAMPLES=OFF \
    -DOFT_BUILD_PYTHON=ON \
    -DOFT_USE_MPI=ON \
    -DOFT_USE_OpenMP=ON \
    -DOFT_PY_KERNEL=python3 \
    -DBLA_VENDOR=OpenBLAS; then
    echo "direct_config_status=0"
    if "${CMAKE_BIN}" --build "${BUILD_DIR}" --parallel "${SLURM_CPUS_PER_TASK:-8}"; then
      echo "direct_build_status=0"
      direct_ok=1
    else
      echo "direct_build_status=1"
    fi
  else
    echo "direct_config_status=1"
    echo "direct_build_status=99"
  fi
} 2>&1 | tee "${BUILD_LOG}"

if [ "${direct_ok}" -ne 1 ]; then
  echo "Direct build failed; attempting source-dependency fallback." | tee -a "${BUILD_LOG}"
  rm -rf "${BUILD_DIR}"
  (
    cd builds
    export CC=gcc
    export CXX=g++
    export FC=gfortran
    export CMAKE="${CMAKE_BIN}"
    build_libs_cmd=(
      python ../src/utilities/build_libs.py
      --build_umfpack=1
      --build_superlu=1
      --no_dl_progress
      --nthread="${SLURM_CPUS_PER_TASK:-8}"
      --build_arpack=1
      --oft_build_tests=1
      --build_mpich=1
    )
    if ! "${build_libs_cmd[@]}"; then
      if [ -f superlu-7_0_0/lib64/libsuperlu.a ] && [ ! -e superlu-7_0_0/lib ]; then
        echo "SuperLU installed into lib64; adding lib -> lib64 compatibility symlink and retrying build_libs.py."
        ln -s lib64 superlu-7_0_0/lib
      fi
      "${build_libs_cmd[@]}"
    fi
    bash config_cmake.sh
    "${CMAKE_BIN}" --build build_release --parallel "${SLURM_CPUS_PER_TASK:-8}"
  ) 2>&1 | tee -a "${BUILD_LOG}"
  BUILD_DIR="${REPO_DIR}/builds/build_release"
fi

(
  cd "${BUILD_DIR}/tests/physics"
  python -m pytest -q \
    test_alfven_2d.py \
    test_sound_2d.py \
    test_drag_decay.py
) 2>&1 | tee "${TEST_LOG}"

{
  echo "oft_build_status=PASS"
  echo "build_dir=${BUILD_DIR}"
  echo "timestamp_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} | tee "${STATUS_FILE}"
