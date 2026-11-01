# GSE168041: replicated cochlear single-cell analysis

The adopted analysis contains 31,556 spiral-ganglion-region cells and 33,149 lateral-wall cells, with eight libraries per compartment and four biological samples per condition. Thirteen cell-type/compartment pseudobulk models supply the manuscript's sample-level differential-expression and enrichment results.

## Contents

- `archived_run/00_manifest/`: biological-sample manifest.
- `archived_run/14_methods/`: QC, annotation, pseudobulk, DESeq2 and enrichment parameters, software versions and original executed scripts.
- `complete_results/`: full enrichment results, focused program summaries, QA and the frozen MSigDB 2026.1.Mm membership/collection records.
- `scripts/`: an offline enrichment entry point and staging helper using the exact frozen gene sets.

Gene-level tables, QC summaries and exact figure inputs are under `manuscript/source_inputs/GSE168041/`. The release provides both annotated AnnData objects and the pseudobulk/differential-support bundle; see [`data/README.md`](../../data/README.md) and [`docs/REPRODUCING.md`](../../docs/REPRODUCING.md) for retrieval and execution commands.

The original scripts document the full analysis route. The offline enrichment adapter replaces network retrieval with hash-verified local membership; its ranking, enrichment and multiple-testing calculations preserve the original code. Input staging and R syntax checks were executed for this archive. The delivered single-cell and enrichment values are the adopted saved analysis, verified against its source files.
