# Reproduction guide

## 1. Verify the archive

```bash
git clone https://github.com/pdx12320/cochlear-noise-transcriptomics.git
cd cochlear-noise-transcriptomics
git checkout v1.2.2
python3 scripts/verify_archive.py
```

This verifies repository-file checksums and structural invariants. It does not fit statistical models or download datasets.

## 2. Rebuild the manuscript and figures

Use an existing TeX installation with the packages listed by the LaTeX source. The files were built with Tectonic 0.17.0. `main.tex` builds the Scientific Reports manuscript, and `main_scientific_reports.tex` is an alias. The journal reference style is supplied as `scientific_reports.bst`.

```bash
cd manuscript
tectonic main.tex
tectonic supplementary.tex
tectonic cover_letter.tex
```

Figure rendering uses Python, NumPy, pandas and ReportLab, with Poppler for raster export. Exact versions are recorded in `manuscript/code/`. Run in a copy of the release if you want to preserve its delivered figure bytes:

```bash
python3 code/redraw_figures.py
python3 code/plot_cell_atlases.py
python3 code/plot_sample_robustness.py
python3 code/plot_annotation_sensitivity.py
```

## 3. GSE168041

Archived scripts are in `analysis/GSE168041/archived_run/14_methods/scripts/`. Parameters and environment versions are in the adjoining method records. The original single-cell route starts from GEO's `GSE168041_RAW.tar`; its immutable source hash and sample manifest are supplied. Each Python stage accepts a project-root argument, and the R stages accept the root as their first argument. Use a new working directory.

Data-only collection v1.0.0 supplies annotated AnnData objects and the pseudobulk/model bundle. The bundle preserves `07_replicate_gate/`, `08_pseudobulk/` and `09_differential/`. These give an entry point at sample-level differential expression without repeating QC and annotation. Create output directories `12_figures/`, `13_source_data/` and `14_methods/`, then run `phase6_pydeseq2.py --root YOUR_NEW_PROJECT_ROOT` using the documented PyDESeq2 environment.

The offline enrichment entry point reads the exact 2026.1.Mm Hallmark/Reactome membership supplied under `analysis/GSE168041/complete_results/` and checks its hash. The analysis code from ranking through all 13 models and FDR adjustment is unchanged from the original script. Its interface was adapted for local files and a new output directory. With the documented R/fgsea dependencies already installed:

```bash
python3 analysis/GSE168041/scripts/prepare_phase7_offline_inputs.py \
  manuscript/source_inputs/GSE168041 reproduced/GSE168041
Rscript --vanilla analysis/GSE168041/scripts/phase7_fgsea_offline.R \
  reproduced/GSE168041 \
  analysis/GSE168041/complete_results/MSIGDB_2026.1.Mm_Hallmark_Reactome_TERM2GENE.tsv.gz \
  analysis/GSE168041/complete_results/MSIGDB_COLLECTIONS_SNAPSHOT.tsv
```

The staging helper was executed and checked during archive preparation. The full 13-model enrichment was not rerun for this archive; its supplied complete results remain the manuscript's GSE168041 results. A new run records its own versions and outputs.

## 4. GSE252285

The existing selected-branch module remains self-contained under `analysis/GSE252285/`. Its README provides the preprocessing and descriptive-analysis commands. Its archived full-stage workflow also contains historical branches; the manuscript uses only the selected 36,111-cell branch and its stored broad labels.

The exact annotated Seurat object is available in the v1.0.0 release assets. `ucell_programs/` and `matched_null/` contain the original UCell and matched-null scripts and their parameters; `complete_results/` contains full primary and gene-length-sensitivity results. Library identities are supported by `library_provenance/`.

## 5. GSE160639

See `analysis/GSE160639/README.md` for the executed current bulk workflow, exact input hashes, fixed pathway membership, environment and entry-point command. The public release contains the adopted current results and the inputs needed to reproduce them. New output paths are required for reruns; source inputs are read-only.

## Scope of validation

The release records which stages were executed for this version and which were retrieved from earlier analyses. Successful document compilation and numerical consistency checks are reported separately from execution of upstream scientific workflows. See `VALIDATION.md`.


## Sensitivity workflows (data collection v1.1.0)

The new modules use existing recorded environments and write to new directories:

- `analysis/GSE168041/sensitivity/README.md` describes the 48 leave-one-sample-out fits, fixed original gene masks and full enrichment testing families. Its full release asset includes all frozen inputs and fit outputs.
- `analysis/GSE252285/annotation_sensitivity/README.md` describes the cluster-19 sensitivity, fixed cell set and paired 10,000-draw analyses. Its companion asset restores full null distributions and plotting inputs; original objects remain in v1.0.0.

The manuscript's stored embeddings are exported by `manuscript/code/export_gse168041_embeddings.py` from the baseline objects; no coordinates were recomputed. Figure rendering uses the standalone commands in `manuscript/Makefile` and the frozen `revision_inputs/` and `source_inputs/` tables. All 14 figure PDFs have matching SVG and 600-dpi PNG exports.

These workflows were executed in the recorded environments. Clean installation of every archived environment is a separate reproduction step; the release records identify the environments actually used and the independent numerical checks performed.
