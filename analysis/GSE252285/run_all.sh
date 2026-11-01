#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SCRIPT_DIR="${SCRIPT_DIR:-$repo_root/scripts/no_doublet_sensitivity}"
exec "$repo_root/scripts/no_doublet_sensitivity/run_no_doublet_robust_qc_tsne.sh" "$@"
