# Environment

The selected no-doublet/condition workflow used R 4.4.3 and the package versions in `software_versions.tsv`; an executed R session is in `sessionInfo.txt`. `install_packages.R` documents package sources but does not guarantee a byte-identical future environment. The original project-local runtime is intentionally excluded because it is large and platform-specific.

For the closest rerun:

1. Use R 4.4.x on Linux.
2. Install the listed CRAN and Bioconductor packages.
3. Compare installed versions to `software_versions.tsv`.
4. Set `R_LIBS_USER` to a project-local library.
5. Record a new `sessionInfo()` with the regenerated results.

Python is used only for small document/table utilities. `requirements.txt` lists the sole non-standard-library dependency.

The optional historical upstream stage1 script also contains an scDblFinder audit and therefore needs `scDblFinder`, `SingleCellExperiment` and `BiocParallel`. Those packages are not required to rerun the selected no-doublet branch from its audited raw object/QC metadata and are not listed as current required dependencies.
