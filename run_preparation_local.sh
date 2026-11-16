#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
ENV_DIR="${TALOS_ENV_DIR:-/faststorage/project/reanalyses_auh/env/talos2_env}"

export PATH="${ENV_DIR}/bin:${PATH}"

cd "${SCRIPT_DIR}"

# A named selector overrides the module's explicit four-CPU declaration.
CPU_CONFIG="$(mktemp "${TMPDIR:-/tmp}/talos-preparation.XXXXXX.config")"
trap 'rm -f "${CPU_CONFIG}"' EXIT
printf '%s\n' 'process { withName: ResummariseRawSubmissions { cpus = 2 } }' > "${CPU_CONFIG}"

nextflow \
  -c nextflow.config \
  -c "${CPU_CONFIG}" \
  run preparation.nf \
  --processed_annotations processed_annotations \
  --large_files large_files \
  "$@"
