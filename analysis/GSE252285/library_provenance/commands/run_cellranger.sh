#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C
export LANG=C
PROJECT="${PROJECT:-/path/to/cochlear_analysis/GSE252285_raw_reconstruction_20260815}"
REAL_PROJECT="$(readlink -f "$PROJECT")"
CELLRANGER="$PROJECT/software/cellranger-7.1.0/cellranger"
REF_NAME="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["selected_reference"])' "$PROJECT/comparison/reference_feature_audit.json")"
REF="$PROJECT/reference/$REF_NAME"
LOG="$PROJECT/logs/run_cellranger.log"
exec > >(tee -a "$LOG") 2>&1

mkdir -p "$PROJECT/cellranger" "$PROJECT/logs"
python3 -c 'import json,sys; p=sys.argv[1]; d=json.load(open(p)); assert d["HHL_TECHNICAL_PARTITIONS_CONFIRMED"] == "YES"' \
    "$PROJECT/qc/run_partition_status.json"
{
    echo "timestamp=$(date -Is)"
    echo "hostname=$(hostname)"
    echo "project=$PROJECT"
    echo "real_project=$REAL_PROJECT"
    echo "localcores=32"
    echo "localmem_gb=180"
    echo "cell_calling=auto (no --expect-cells/--force-cells)"
    echo "chemistry=auto"
    echo "bam_generation=disabled via Cell Ranger 7.1.0 --no-bam (output-only storage choice)"
    echo "reference=$REF_NAME (selected by exact official feature-pair set match)"
    nproc
    free -h
    df -hT "$REAL_PROJECT"
    uptime
    "$CELLRANGER" --version
} > "$PROJECT/commands/cellranger_environment.txt"
cat "$PROJECT/commands/cellranger_environment.txt"

run_library() {
    local library="$1"
    local metrics="$PROJECT/cellranger/$library/outs/metrics_summary.csv"
    if [[ -s "$metrics" ]]; then
        echo "[$(date -Is)] validated existing completed count: $library"
    else
        echo "[$(date -Is)] count start: $library"
        local attempt
        for attempt in 1 2; do
            echo "[$(date -Is)] count attempt $attempt/2: $library"
            if (
                cd "$PROJECT/cellranger"
                "$CELLRANGER" count \
                --id="$library" \
                --transcriptome="$REF" \
                --fastqs="$PROJECT/fastq/cellranger_input/$library" \
                --sample="$library" \
                --chemistry=auto \
                --localcores=32 \
                --localmem=180 \
                --no-bam
            ) 2>&1 | tee -a "$PROJECT/logs/cellranger_${library}.log"; then
                break
            fi
            if [[ "$attempt" -eq 2 ]]; then
                echo "Cell Ranger failed twice for $library" >&2
                return 1
            fi
            echo "[$(date -Is)] resumable retry in 60 seconds: $library"
            sleep 60
        done
    fi

    for required in \
        "$PROJECT/cellranger/$library/outs/web_summary.html" \
        "$PROJECT/cellranger/$library/outs/metrics_summary.csv" \
        "$PROJECT/cellranger/$library/outs/raw_feature_bc_matrix/barcodes.tsv.gz" \
        "$PROJECT/cellranger/$library/outs/raw_feature_bc_matrix/features.tsv.gz" \
        "$PROJECT/cellranger/$library/outs/raw_feature_bc_matrix/matrix.mtx.gz" \
        "$PROJECT/cellranger/$library/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" \
        "$PROJECT/cellranger/$library/outs/filtered_feature_bc_matrix/features.tsv.gz" \
        "$PROJECT/cellranger/$library/outs/filtered_feature_bc_matrix/matrix.mtx.gz"; do
        [[ -s "$required" ]] || { echo "Missing required output: $required" >&2; return 1; }
    done

    find "$PROJECT/cellranger/$library/outs/raw_feature_bc_matrix" \
         "$PROJECT/cellranger/$library/outs/filtered_feature_bc_matrix" \
         -type f -print0 | sort -z | xargs -0 sha256sum \
         > "$PROJECT/cellranger/$library/matrix_outputs.sha256"
    sha256sum \
        "$PROJECT/cellranger/$library/outs/web_summary.html" \
        "$PROJECT/cellranger/$library/outs/metrics_summary.csv" \
        >> "$PROJECT/cellranger/$library/matrix_outputs.sha256"
    echo "[$(date -Is)] count complete and validated: $library"
}

for library in A1_Control B1_HHL C1_NIHL; do
    run_library "$library"
done

echo "[$(date -Is)] all Cell Ranger counts complete"
