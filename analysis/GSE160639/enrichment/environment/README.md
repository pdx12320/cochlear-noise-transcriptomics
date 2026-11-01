# Targeted-enrichment runtime

The executed analysis used the existing Linux x86-64 Conda environment with R 4.5.3 and fgsea 1.36.2. The 147 installed Conda package records are frozen by exact URLs, builds, MD5 and SHA-256 values. `INSTALLED_R_PACKAGES.tsv` lists all installed R packages with versions, licenses, package DESCRIPTION repositories and source-package metadata where recorded.

fgsea was installed as the Bioconda package `bioconductor-fgsea-1.36.2-r45ha27e39d_0.conda`, SHA-256 `67241d03574f4a8a8750c75a6fb6ed67b4c5bf0a29f9a7cd4b53e99e64330f7d`. Its DESCRIPTION identifies Bioconductor 3.22. R was the conda-forge package `r-base-4.5.3-h502d0c9_3.conda`. These records were obtained read-only after the current run; no package was installed or scientific analysis rerun during this export.

From this `environment/` directory, create a new local runtime:

```sh
export MAMBA_ROOT_PREFIX="$PWD/mamba_root"
micromamba create -y -p ./runtime --file CONDA_EXPLICIT.txt
micromamba run -p ./runtime Rscript --vanilla -e \
  'stopifnot(as.character(getRversion())=="4.5.3",as.character(packageVersion("fgsea"))=="1.36.2")'
```

Use this runtime with `../authorized_targeted_bulk_from_rawcounts.R` and the relative input paths in the parent analysis README. The script requires the current DESeq2 tables, locked sample manifest, frozen native mouse gene-set cache, and a new output directory; it checks their hashes before running. The current execution used one process, a fixed seed reset for each of three contrasts, and the documented 27-test BH family.

The exact lock is supplied for restoration; a second clean recreation of this environment has not been executed. SHA-256 values in the Conda records are installed-package metadata, not a claim that every package archive was downloaded and rehashed during this export. `EXECUTION_SESSIONINFO.txt` is the session recorded by the completed current enrichment.
