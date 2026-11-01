# Cochlear noise transcriptomics

Analysis archive for **Cellular organization of inflammatory transcription in the mouse cochlea after acoustic trauma**.

**Dengxiang Yang** · Southern University of Science and Technology, Shenzhen, China  
Correspondence: 12412602@mail.sustech.edu.cn

This study examines cell-type inflammatory transcription using three public mouse cochlear datasets. GSE168041 supplies replicated, sample-level pseudobulk contrasts; GSE160639 supplies whole-cochlea exposure comparisons; GSE252285 supplies descriptive expression localization in 36,111 cells.

## Start here

- [Scientific Reports manuscript PDF](manuscript/main.pdf) and [LaTeX, figures and submission files](manuscript/)
- [Data sources and large analysis objects](data/README.md)
- [Reproduction guide](docs/REPRODUCING.md)
- [Analysis and validation record](docs/VALIDATION.md)
- [Release downloads](https://github.com/pdx12320/cochlear-noise-transcriptomics/releases/tag/v1.2.2)

## Organization

```text
analysis/
  GSE168041/       replicated single-cell workflow and complete enrichment results
  GSE160639/       current bulk workflow, required inputs and targeted enrichment
  GSE252285/       selected single-cell branch, library provenance and matched-null tests
manuscript/        LaTeX, PDFs, vector figures and exact panel source tables
data/              public accessions and release-object download instructions
docs/              methods, reproduction and validation documentation
manifests/         checksums and release inventories
scripts/           archive verification and asset retrieval
CITATION.cff       citation metadata
```

## Study design

| Dataset | Material | Analysis unit | Main contribution |
|---|---|---|---|
| GSE168041 | 64,705 cells from 16 dissected-compartment libraries | Four mice per condition within each compartment | Replicated cell-type pseudobulk responses |
| GSE160639 | 12 whole-cochlea libraries, 70/94/105 dB | Four samples per exposure group | Bulk differential expression and nine pathway contrasts |
| GSE252285 | 36,111 cells across Control/HHL/NIHL libraries | Descriptive pooled-library summaries | Gene and pathway localization |

The primary GSE252285 analysis uses its stored broad annotations. An explicit sensitivity analysis assigns the 319 cells of cluster 19 to fibrocyte/mesenchymal and recalculates gene and program localization. The primary labels, raw objects and sample inclusion remain preserved.

## Current archive

The current manuscript contains six main figures, eight supplementary figures, and 36 references. The Scientific Reports manuscript contains 18 pages with review line numbers; Supplementary Information contains 9 pages and the cover letter contains 1 page. The complete LaTeX package and current PDFs are available in [v1.2.2](https://github.com/pdx12320/cochlear-noise-transcriptomics/releases/tag/v1.2.2).

- All 56 focused pathway estimates retain their direction across 48 leave-one-animal-out fits. Fibrocyte TNF–NF-κB enrichment meets FDR < 0.05 in 7/8 SGN and 6/8 lateral-wall omissions; the other five focused pathway/model combinations meet the criterion in all eight.
- All eleven primary gene-localization assignments persist after cluster 19 reassignment. Apold1 is expressed in endothelial cells and the fibrocyte-associated cluster.
- Both matched-null schemes retain five positive comparisons among the six selected follow-ups; endothelial MAPK/AP-1 remains outside the FDR criterion.

The [v1.0.0 data collection](https://github.com/pdx12320/cochlear-noise-transcriptomics/releases/tag/v1.0.0) supplies the exact analysis objects and runtime sources used by this study. The [v1.1.0 sensitivity collection](https://github.com/pdx12320/cochlear-noise-transcriptomics/releases/tag/v1.1.0) supplies complete animal-omission and annotation-sensitivity records. These are data-only collections for the current paper. All release tags identify the current repository snapshot.

See [GSE168041 sensitivity](analysis/GSE168041/sensitivity/) and [GSE252285 annotation sensitivity](analysis/GSE252285/annotation_sensitivity/) for the executed scripts, parameters, result tables and asset checksums.

## Quick integrity check

```bash
python3 scripts/verify_archive.py
```

## Citation and rights

Use [CITATION.cff](CITATION.cff) and cite the source GEO studies and software listed in the manuscript.

Author-owned code is provided under the [MIT license](LICENSE). Author-created result tables and figures are provided under [CC BY 4.0](licenses/DATA_LICENSE.md). Manuscript text and external data retain their respective rights. See [third-party notices](THIRD_PARTY_NOTICES.md).
