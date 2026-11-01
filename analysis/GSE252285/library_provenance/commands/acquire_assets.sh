#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-/path/to/cochlear_analysis/GSE252285_raw_reconstruction_20260815}"
LOG="$PROJECT/logs/acquire_assets.log"
exec > >(tee -a "$LOG") 2>&1

download() {
    local url="$1"
    local output="$2"
    if [[ -s "$output" ]]; then
        echo "[$(date -Is)] present: $output"
        return 0
    fi
    echo "[$(date -Is)] downloading: $url"
    curl -fsSL --retry 20 --retry-delay 10 --continue-at - \
        --connect-timeout 30 -o "$output.part" "$url"
    mv "$output.part" "$output"
}

mkdir -p "$PROJECT"/{official,reference,software,logs}

BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252285/suppl"
download "$BASE/GSE252285_barcodes.tsv.gz" "$PROJECT/official/GSE252285_barcodes.tsv.gz"
download "$BASE/GSE252285_features.tsv.gz" "$PROJECT/official/GSE252285_features.tsv.gz"
download "$BASE/GSE252285_matrix.mtx.gz" "$PROJECT/official/GSE252285_matrix.mtx.gz"

gzip -t "$PROJECT/official/GSE252285_barcodes.tsv.gz"
gzip -t "$PROJECT/official/GSE252285_features.tsv.gz"
gzip -t "$PROJECT/official/GSE252285_matrix.mtx.gz"

REF_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
REF_TAR="$PROJECT/reference/refdata-gex-GRCm39-2024-A.tar.gz"
download "$REF_URL" "$REF_TAR"
gzip -t "$REF_TAR"
if [[ ! -d "$PROJECT/reference/refdata-gex-GRCm39-2024-A" ]]; then
    echo "[$(date -Is)] extracting reference"
    tar -xzf "$REF_TAR" -C "$PROJECT/reference"
fi

MM10_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz"
MM10_TAR="$PROJECT/reference/refdata-gex-mm10-2020-A.tar.gz"
download "$MM10_URL" "$MM10_TAR"
gzip -t "$MM10_TAR"
if [[ ! -d "$PROJECT/reference/refdata-gex-mm10-2020-A" ]]; then
    echo "[$(date -Is)] extracting mm10 reference"
    tar -xzf "$MM10_TAR" -C "$PROJECT/reference"
fi

CR_PAGE="https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions"
CR_URL="$(
    curl -fsSL --retry 10 --retry-delay 10 --connect-timeout 30 "$CR_PAGE" |
        python3 -c 'import re,sys; t=sys.stdin.read(); m=re.search(r"https://cf\.10xgenomics\.com/releases/cell-exp/cellranger-7\.1\.0\.tar\.gz\?[^\"<]+", t); m or sys.exit("Cell Ranger 7.1.0 signed URL not found"); print(m.group(0).replace("\\u0026", "&").replace("&amp;", "&"))'
)"
CR_TAR="$PROJECT/software/cellranger-7.1.0.tar.gz"
CR_DIR="$PROJECT/software/cellranger-7.1.0"
download "$CR_URL" "$CR_TAR"
gzip -t "$CR_TAR"
if [[ ! -x "$CR_DIR/cellranger" ]]; then
    echo "[$(date -Is)] extracting official Cell Ranger 7.1.0"
    tar -xzf "$CR_TAR" -C "$PROJECT/software"
fi

"$CR_DIR/cellranger" --version

(
    cd "$PROJECT"
    sha256sum official/* reference/refdata-gex-GRCm39-2024-A.tar.gz \
        reference/refdata-gex-mm10-2020-A.tar.gz \
        software/cellranger-7.1.0.tar.gz > checksums_assets.sha256
)

echo "[$(date -Is)] acquire_assets complete"
