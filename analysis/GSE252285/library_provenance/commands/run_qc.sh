#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C
export LANG=C
PROJECT="${PROJECT:-/path/to/cochlear_analysis/GSE252285_raw_reconstruction_20260815}"
LOG="$PROJECT/logs/run_qc.log"
FASTQC=/path/to/analysis_home/.conda/envs/RNA-seq/bin/fastqc
MULTIQC=/path/to/analysis_home/.conda/envs/RNA-seq/bin/multiqc
CELLRANGER_DIR="$PROJECT/software/cellranger-7.1.0"
exec > >(tee -a "$LOG") 2>&1

echo "[$(date -Is)] QC start"
"$FASTQC" --version
"$MULTIQC" --version
FASTQS=()
while IFS= read -r -d '' FASTQ; do
    FASTQS+=("$FASTQ")
done < <(find "$PROJECT/fastq/archive" -type f -name '*.fastq.gz' -print0 | sort -z)
if [[ "${#FASTQS[@]}" -ne 12 ]]; then
    echo "Expected 12 FASTQ files; found ${#FASTQS[@]}" >&2
    exit 1
fi

mkdir -p "$PROJECT/qc/fastqc" "$PROJECT/qc/multiqc"
"$FASTQC" --threads 4 --outdir "$PROJECT/qc/fastqc" "${FASTQS[@]}"
"$MULTIQC" --force --outdir "$PROJECT/qc/multiqc" \
    --filename GSE252285_raw_fastq_multiqc.html "$PROJECT/qc/fastqc"

WHITELIST_SOURCE="$(
    find "$CELLRANGER_DIR" -type f -name '3M-february-2018.txt.gz' -print | head -1
)"
if [[ -z "$WHITELIST_SOURCE" ]]; then
    echo "Cell Ranger v3 whitelist was not found in the official distribution" >&2
    exit 1
fi
WHITELIST="$PROJECT/metadata/3M-february-2018.txt.gz"
if [[ ! -s "$WHITELIST" ]]; then
    cp "$WHITELIST_SOURCE" "$WHITELIST"
fi
python3 "$PROJECT/commands/audit_read_structure.py" "$WHITELIST"
python3 "$PROJECT/commands/audit_run_partitions.py"

find "$PROJECT/qc/fastqc" -type f -print0 | sort -z | xargs -0 sha256sum \
    > "$PROJECT/qc/fastqc_outputs.sha256"
find "$PROJECT/qc/multiqc" -type f -print0 | sort -z | xargs -0 sha256sum \
    > "$PROJECT/qc/multiqc_outputs.sha256"
echo "[$(date -Is)] QC complete"
