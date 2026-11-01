# Current analysis objects in data collection v1.0.0

The release provides the exact selected analysis objects and the current GSE168041 sample-level pseudobulk/model records. The objects were read and checksummed for archival verification; no statistical analysis was rerun and no object labels, thresholds or sample inclusion were changed.

## Download and verify

From the repository root, using Python 3.8 or later and only its standard library:

```bash
python3 scripts/fetch_release_objects.py \
  --assets-dir data/release_downloads \
  --output-dir data/objects
```

The script downloads the six data assets listed in `manifests/RELEASE_OBJECTS.tsv`, verifies each asset's size and SHA256, and verifies the complete object after reassembly. Interrupted or corrupt downloads are not promoted to final filenames. An existing correct output is reused; an existing incorrect output is retained and reported as an error. Allow approximately 8 GB of free disk space when retaining both downloaded parts and reconstructed objects.

For assets downloaded separately:

```bash
python3 scripts/reassemble_release_objects.py \
  --manifest manifests/RELEASE_OBJECTS.tsv \
  --assets-dir data/release_downloads \
  --output-dir data/objects
```

The standalone scripts attached to the release require an explicit `--manifest RELEASE_OBJECTS.tsv` argument. Keep both scripts together if using the downloader.

## Objects and exact metadata interpretation

| Object | Cells | Features | Release layout |
|---|---:|---:|---|
| GSE252285 selected annotated Seurat object | 36,111 | 32,285 | Three parts, reconstructed into `GSE252285/GSE252285_36111cells_annotated.rds` |
| GSE168041 spiral-ganglion annotated AnnData object | 31,556 | 27,998 | `GSE168041/SGN_annotated.h5ad` |
| GSE168041 lateral-wall annotated AnnData object | 33,149 | 27,998 | `GSE168041/Lateral_wall_annotated.h5ad` |

The GSE252285 object retains the recorded metadata fields:

- Use `seurat_clusters_no_doublet` for the selected 28-cluster result at resolution 0.6. The generic `seurat_clusters` field has 37 labels from the final resolution-audit run and does not identify the selected main clustering.
- Use `broad_annotation_no_doublet_sensitivity` for the broad annotations used by the manuscript: myeloid 27,676; unresolved 4,726; fibrocyte/mesenchymal 2,997; endothelial 275; glial 264; supporting/epithelial 149; hair/sensory 24.
- The object’s `condition` column remains `UNKNOWN`. The resolved condition mapping is supplied separately in the analysis tables: `LIB3` / barcode suffix `-3` = Control (10,273 cells); `LIB1` / `-1` = HHL (16,631); `LIB2` / `-2` = NIHL (9,207). Apply that documented mapping when constructing descriptive summaries; do not infer conditions from numerical library order.

These fields were inspected without rewriting the object. The public object hashes equal the downloaded original hashes and the recorded remote hashes.

## Current GSE168041 supporting records

`GSE168041_pseudobulk_differential_support.tar.gz` contains 97 current files from `07_replicate_gate`, `08_pseudobulk` and `09_differential`, including 13 fitted DDS AnnData records. Each fitted model has eight sample observations. The archive also contains `GSE168041/RELEASE_CONTENTS.tsv` with per-file hashes. The same inventory is available as `manifests/GSE168041_SUPPORT_CONTENTS.tsv`.

Extract the verified archive into a new work directory:

```bash
mkdir -p restored_support
python3 -m tarfile -e \
  data/objects/GSE168041_pseudobulk_differential_support.tar.gz \
  restored_support
```

Archive members use generic relative names, numeric ownership 0, blank user/group names and fixed timestamps. Only deployment-specific paths in `result_file` (13 rows of `S1_DS_SUMMARY.tsv`) and `directory` (16 rows of `S1_PSEUDOBULK_AUDIT.tsv`) were changed to dataset-relative paths in the publication copies. All other table values are identical. Original and public hashes are recorded in the content manifest; the original downloads are preserved.

## Scope

The release contains the selected current results. It does not include superseded analysis releases or unexecuted proposals. Large public sequencing inputs remain available through the cited GEO/SRA accessions. These files support reproducibility of the archived work; the executed GSE160639 bulk workflow is documented separately in `analysis/GSE160639/`. Licensing and source attribution are defined in the repository license files and `THIRD_PARTY_NOTICES.md`.
