cran_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "data.table", "dplyr", "ggplot2",
  "patchwork", "ragg", "future", "scales", "Rtsne", "png"
)

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran)

versions <- c("R" = R.version.string, vapply(
  cran_packages,
  function(pkg) as.character(utils::packageVersion(pkg)),
  character(1)
))
write.table(
  data.frame(software = names(versions), version = unname(versions)),
  file = "installed_versions.tsv", sep = "\t", quote = FALSE, row.names = FALSE
)
