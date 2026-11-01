# Portable entry points

These copies replace deployment-specific paths with arguments and refuse an existing output directory. The scientific body beginning with `write_tsv` is unchanged; the adaptation record identifies its SHA256. The original scripts remain in the adjacent module directories.

Use the supplied object and frozen GMT files. UCell/Seurat package versions are listed in the original `outputs/SOFTWARE_VERSIONS.tsv` records. The full analyses were not rerun while adapting their file interfaces.

From the repository root, with an existing compatible R environment:

```bash
Rscript --vanilla analysis/GSE252285/portable_scripts/run_gse252285_ucell.R \
  data/objects/GSE252285/GSE252285_36111cells_annotated.rds \
  analysis/GSE252285/complete_results reproduced/GSE252285_ucell

N_ITER=10000 Rscript --vanilla analysis/GSE252285/portable_scripts/run_gse252285_matched_null.R \
  data/objects/GSE252285/GSE252285_36111cells_annotated.rds \
  analysis/GSE252285/complete_results \
  /path/to/gencode.vM23.annotation.gtf.gz reproduced/GSE252285_matched_null
```

The GENCODE reference is available from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz . Its MD5 is `b646883b82252b8831d337e9b88d1f30`. A 100-iteration benchmark can be requested with `N_ITER=100`; manuscript results use 10,000 iterations.

Use `seurat_clusters_no_doublet` for the selected 28-cluster result. The object's generic `seurat_clusters` records a separate resolution audit. Use the supplied external library-condition crosswalk; the stored generic `condition` column is unassigned.
