#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C
export LANG=C
PROJECT="${PROJECT:-/path/to/cochlear_analysis/GSE252285_raw_reconstruction_20260815}"
LOG="$PROJECT/logs/run_all.log"
exec > >(tee -a "$LOG") 2>&1

stage() {
    echo "[$(date -Is)] STAGE_START $1"
}

stage acquire_assets
if [[ -s "$PROJECT/checksums_assets.sha256" \
      && -x "$PROJECT/software/cellranger-7.1.0/cellranger" \
      && -d "$PROJECT/reference/refdata-gex-GRCm39-2024-A" \
      && -d "$PROJECT/reference/refdata-gex-mm10-2020-A" ]]; then
    "$PROJECT/software/cellranger-7.1.0/cellranger" --version
    echo "[$(date -Is)] STAGE_REUSE acquire_assets (completed checksum manifest present)"
else
    "$PROJECT/commands/acquire_assets.sh"
fi

stage audit_official_reference
/path/to/analysis_home/.conda/envs/RNA-seq/bin/Rscript \
    "$PROJECT/commands/summarize_official_matrix.R" \
    | tee "$PROJECT/logs/summarize_official_matrix.log"
python3 "$PROJECT/commands/audit_official_reference.py" \
    | tee "$PROJECT/logs/audit_official_reference.log"

stage download_fastqs
python3 "$PROJECT/commands/download_fastqs.py" \
    | tee "$PROJECT/logs/download_fastqs.log"

stage raw_fastq_qc_and_partition_audit
"$PROJECT/commands/run_qc.sh"

stage cellranger_counts
"$PROJECT/commands/run_cellranger.sh"

stage barcode_and_metrics_comparison
python3 "$PROJECT/commands/compare_barcodes_metrics.py" \
    | tee "$PROJECT/logs/compare_barcodes_metrics.log"

stage final_validation_and_summary
python3 "$PROJECT/commands/build_result_summary.py" \
    | tee "$PROJECT/logs/build_result_summary.log"

echo "[$(date -Is)] PIPELINE_SUCCESS"
