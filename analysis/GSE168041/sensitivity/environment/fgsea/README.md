# Enrichment runtime

The current animal-deletion analyses used the existing Linux x86-64 R 4.5.3/fgsea 1.36.2 environment. All 147 installed Conda packages are locked by source URL, build and MD5, with SHA-256 also recorded; all installed R package versions and available source metadata are supplied.

From this directory, use a new prefix: `micromamba create -y -p ./fgsea_env --file CONDA_EXPLICIT.txt`. The R entry point is `./fgsea_env/bin/Rscript`. A second clean environment recreation was not executed during this revision.

The existing runtime was verified before all 48 fits. `EXECUTION_SESSIONINFO.txt` is taken from a completed sensitivity enrichment, with its machine-local BLAS path replaced by a generic environment token. Scientific parameters are in the executed script and analysis metadata.
