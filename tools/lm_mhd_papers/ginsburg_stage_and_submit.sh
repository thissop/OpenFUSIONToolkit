#!/usr/bin/env bash
set -euo pipefail

REMOTE_ROOT="${REMOTE_ROOT:-/burg/astro/users/tjk2147/spf_work/oft_lm_mhd_papers}"
REMOTE_REPO="${REMOTE_ROOT}/repo"
SSH_HOST="${SSH_HOST:-ginsburg}"

cd "$(git rev-parse --show-toplevel 2>/dev/null || pwd)"

LOCAL_BRANCH="$(git branch --show-current 2>/dev/null || true)"
LOCAL_HEAD="$(git rev-parse --short HEAD 2>/dev/null || true)"
LOCAL_STATUS_COUNT="$(git status --short 2>/dev/null | wc -l | tr -d ' ')"
STAGED_TIMESTAMP_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

echo "Staging repo to ${SSH_HOST}:${REMOTE_REPO}"
ssh -o BatchMode=yes "${SSH_HOST}" "rm -rf '${REMOTE_REPO}' && mkdir -p '${REMOTE_REPO}'"

COPYFILE_DISABLE=1 tar \
  --exclude='.git' \
  --exclude='builds' \
  --exclude='.claude' \
  --exclude='deliverables/lm_mhd_for_yuchen/results' \
  --exclude='*.pyc' \
  --exclude='__pycache__' \
  -czf - . | ssh -o BatchMode=yes "${SSH_HOST}" "tar xzf - -C '${REMOTE_REPO}'"

ssh -o BatchMode=yes "${SSH_HOST}" "mkdir -p '${REMOTE_REPO}/cases/lm_mhd_papers/results' && {
  echo staged_timestamp_utc='${STAGED_TIMESTAMP_UTC}';
  echo staged_branch='${LOCAL_BRANCH}';
  echo staged_head='${LOCAL_HEAD}';
  echo staged_dirty_status_count='${LOCAL_STATUS_COUNT}';
} > '${REMOTE_REPO}/cases/lm_mhd_papers/results/staging_metadata.txt'"

echo "Submitting Python/reference baseline"
ssh -o BatchMode=yes "${SSH_HOST}" \
  "cd '${REMOTE_REPO}' && mkdir -p logs && sbatch --export=ALL,REPO_DIR='${REMOTE_REPO}' tools/lm_mhd_papers/ginsburg_python_baseline.sbatch"

echo "Submitting OFT build/regression baseline"
ssh -o BatchMode=yes "${SSH_HOST}" \
  "cd '${REMOTE_REPO}' && mkdir -p logs && sbatch --export=ALL,REPO_DIR='${REMOTE_REPO}' tools/lm_mhd_papers/ginsburg_oft_build_baseline.sbatch"

echo "Remote status:"
ssh -o BatchMode=yes "${SSH_HOST}" "squeue -u \"\$USER\" -o '%.18i %.9P %.24j %.8T %.10M %.10l %.6D %R' | head -n 20"
