# Current bulk reconstruction and reproducibility files

The current analysis starts from the public GSE160639 count matrix, keeps the locked 12 samples and 19,673-gene universe, estimates new median-ratio size factors, and fits the specified DESeq2 models. The three new Wald tables feed the separately archived fixed nine-pathway by three-contrast enrichment workflow. Only current outputs are included here.

- `inputs/`: the three checksum-locked portable inputs; no previous statistical table is an analysis dependency.
- `scripts/`: the executed DESeq2 analysis and the publication-object metadata sanitizer.
- `results/`: four complete current gene-result tables, normalized counts, size factors, diagnostics, parameters and execution session.
- `environment/`: exact Conda lock, R versions/licenses, 48 source-archive hashes, and restore helper/instructions. The source archives are supplied as the matching Release asset.
- `docs/`: independent numerical checks and publication-object equivalence records.

The matching fitted-object Release asset contains safe publication copies. Design and dispersion-function environments were cleaned to remove captured machine-local execution arguments; numerical data, model metadata, dispersion predictions, and re-extracted results were verified identical. The original execution objects are retained privately.

The implemented workflow preserves fixed gene membership and annotations as explicit inputs and estimates fresh median-ratio size factors. GSE168041 and GSE252285 analyses are unchanged. The analysis was executed and audited; a second clean environment restore has not been executed.

Public raw counts: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160639
Runtime sources: https://cran.r-project.org/src/base/R-4/ and https://bioconductor.org/packages/release/bioc/html/DESeq2.html
