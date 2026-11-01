#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C
export LANG=C
PROJECT="${PROJECT:-/path/to/cochlear_analysis/GSE252285_raw_reconstruction_20260815}"
LOG="$PROJECT/logs/resume_from_cellranger.log"
exec > >(tee -a "$LOG") 2>&1

stage() {
    echo "[$(date -Is)] STAGE_START $1"
}

stage cellranger_counts
"$PROJECT/commands/run_cellranger.sh"

stage barcode_and_metrics_comparison
python3 "$PROJECT/commands/compare_barcodes_metrics.py" \
    | tee "$PROJECT/logs/compare_barcodes_metrics.log"

stage final_validation_and_summary
python3 "$PROJECT/commands/build_result_summary.py" \
    | tee "$PROJECT/logs/build_result_summary.log"

echo "[$(date -Is)] PIPELINE_SUCCESS"
