#!/usr/bin/env bash
set -euo pipefail
umask 002

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

FRONTEND_PORT="${REFLEX_FRONTEND_PORT:-3100}"
BACKEND_PORT="${REFLEX_BACKEND_PORT:-8100}"
DEFAULT_REFLEX_CONDA_PREFIX="/faststorage/project/reanalyses_auh/env/talos2_env"
export PATH="${DEFAULT_REFLEX_CONDA_PREFIX}/bin:${PATH}"
REFLEX_ARGS=(
  --frontend-port "$FRONTEND_PORT"
  --backend-port "$BACKEND_PORT"
  "$@"
)

echo "Starting TALOS AARGUS TEST RUNNER"
echo "Frontend: http://localhost:${FRONTEND_PORT}"
echo "Backend:  http://localhost:${BACKEND_PORT}"

if [[ -n "${REFLEX_CONDA_ENV:-}" ]]; then
  conda run -n "$REFLEX_CONDA_ENV" reflex run "${REFLEX_ARGS[@]}"
elif [[ -n "${REFLEX_CONDA_PREFIX:-}" ]]; then
  conda run -p "$REFLEX_CONDA_PREFIX" reflex run "${REFLEX_ARGS[@]}"
elif [[ -x "${DEFAULT_REFLEX_CONDA_PREFIX}/bin/reflex" ]]; then
  conda run -p "${DEFAULT_REFLEX_CONDA_PREFIX}" reflex run "${REFLEX_ARGS[@]}"
else
  echo "[ERROR] Reflex is not installed in ${DEFAULT_REFLEX_CONDA_PREFIX}." >&2
  echo "Install it there or set REFLEX_CONDA_ENV/REFLEX_CONDA_PREFIX explicitly." >&2
  exit 1
fi
