# Restore the bulk-analysis runtime

The executed bulk fit used Linux x86-64, R 4.6.1 and DESeq2 1.52.0. `CONDA_EXPLICIT.txt` locks the 95 realized Conda packages by URL and MD5; `CONDA_PACKAGE_SHA256.tsv` and `CONDA_PACKAGE_RECORDS.json` also record SHA-256. `INSTALLED_R_PACKAGES.tsv` records the exact R versions and licenses, and `PACKAGE_SOURCE_SHA256.tsv` locks 48 source archives. Original package archives retain their licensing metadata.

Run the following from `analysis/GSE160639/deseq2/environment/`, with micromamba 2.9.0 available. Use new runtime and library directories.

```sh
export MAMBA_ROOT_PREFIX="$PWD/mamba_root"
micromamba create -y -p ./runtime --file CONDA_EXPLICIT.txt
curl -fL -o Bulk_R_Source_Packages_20261101.tar \
  https://github.com/pdx12320/cochlear-noise-transcriptomics/releases/download/v1.0.0/Bulk_R_Source_Packages_20261101.tar
printf '%s  %s\n' \
  80ad00543a4804f5b51c6c4964fd64cfed377f747995916772611e6f41e0fc4c \
  Bulk_R_Source_Packages_20261101.tar | sha256sum -c -
tar -xf Bulk_R_Source_Packages_20261101.tar
micromamba run -p ./runtime Rscript --vanilla restore_locked_r_packages.R \
  ./package_sources ./PACKAGE_SOURCE_SHA256.tsv ./R_library
```

The [source-package Release asset](https://github.com/pdx12320/cochlear-noise-transcriptomics/releases/download/v1.0.0/Bulk_R_Source_Packages_20261101.tar) already contains the `package_sources/` directory and its SHA manifest. Extract it in the current environment directory, as shown. The optional restore helper verifies all source hashes, derives the dependency order from DESCRIPTION records, and installs exact local archives without repository access. It uses one package job and up to four compiler threads. Its R syntax was checked; a second clean runtime restore has not been executed. The executed original installation used `install_exact_bioc.R` and official CRAN/Bioconductor 3.23 sources, with exact R/DESeq2 version guards.

The analysis requires three locked inputs and a nonexistent output directory. Run its checksum and identity preflight first:

```sh
R_LIBS_USER="$PWD/R_library" R_LIBS_SITE= \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
micromamba run -p ./runtime Rscript --vanilla ../scripts/authorized_bulk_deseq2_reestimate.R \
  --preflight-only ../inputs/GSE160639_RawReadCount.txt.gz ../inputs ./new_results
```

After the preflight passes, execute the same analysis in the same directory:

```sh
R_LIBS_USER="$PWD/R_library" R_LIBS_SITE= \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
micromamba run -p ./runtime Rscript --vanilla ../scripts/authorized_bulk_deseq2_reestimate.R \
  --execute-approved ../inputs/GSE160639_RawReadCount.txt.gz ../inputs ./new_results
```

No authorization-record file is an analysis input. The recorded fit ran as one statistical process with an 8 GiB virtual-memory limit and a 1,200-second timeout, and completed in 37.927 seconds. It estimates median-ratio size factors from the fixed 19,673-gene matrix. Its three new Wald tables are the inputs to the separately documented targeted-enrichment workflow.

The explicit Conda package hashes are installed-package records. Six surviving cached archives were independently rehashed; the other 89 archives were not independently re-downloaded. All 48 frozen R source archives were independently rehashed and matched. Publication metadata uses relative paths. Runtime installation/solver-attempt logs remain outside the public package.
