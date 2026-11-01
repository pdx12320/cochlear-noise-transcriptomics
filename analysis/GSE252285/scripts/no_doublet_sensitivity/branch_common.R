options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

project_dir <- normalizePath(
  Sys.getenv("PROJECT_DIR", Sys.getenv("GSE252285_PROJECT_ROOT", unset = ".")),
  mustWork = TRUE
)
branch_name <- "no_doublet_robust_qc_tsne_sensitivity"
branch_path <- function(section, ...) file.path(project_dir, section, branch_name, ...)

branch_sections <- c("00_manifest", "04_tables", "05_figures", "06_objects",
                     "07_logs", "08_methods", "09_qc", "10_release")
invisible(lapply(branch_sections, function(x) {
  dir.create(branch_path(x), recursive = TRUE, showWarnings = FALSE)
}))

seed <- 42L
set.seed(seed)

legacy_setting <- tolower(Sys.getenv("RUN_LEGACY_COMPARISON", "false"))
if (!legacy_setting %in% c("true", "false")) {
  stop("RUN_LEGACY_COMPARISON must be true or false")
}
run_legacy_comparison <- identical(legacy_setting, "true")

write_tsv <- function(x, path) {
  data.table::fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

theme_branch <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey75"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA)
    )
}

save_branch_dual <- function(plot, stem, width = 8, height = 6, dpi = 300) {
  pdf_path <- branch_path("05_figures", paste0(stem, ".pdf"))
  png_path <- branch_path("05_figures", paste0(stem, ".png"))
  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height,
                  units = "in", limitsize = FALSE, bg = "white")
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height,
                  units = "in", dpi = dpi, device = ragg::agg_png,
                  limitsize = FALSE, bg = "white")
  invisible(c(pdf_path, png_path))
}

stage_message <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

safe_zscore <- function(x) {
  z <- as.numeric(scale(x))
  z[!is.finite(z)] <- 0
  z
}

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  a <- sum(choose2(rowSums(tab)))
  b <- sum(choose2(colSums(tab)))
  ab <- sum(choose2(tab))
  total <- choose2(n)
  expected <- a * b / total
  denom <- 0.5 * (a + b) - expected
  if (denom == 0) return(NA_real_)
  (ab - expected) / denom
}

broad_palette <- c(
  Hair_sensory = "#D55E00",
  Supporting_epithelial = "#E69F00",
  Endothelial_vascular = "#0072B2",
  Immune_myeloid = "#CC79A7",
  Fibrocyte_mesenchymal = "#009E73",
  Neuronal = "#56B4E9",
  Glial = "#F0E442",
  Other = "#7A7A7A",
  UNRESOLVED = "#000000"
)

top20_fixed <- c(
  "Dusp1", "Fos", "Nfkbia", "Arrdc2", "Apold1", "Sik1", "Per1", "Ccn1",
  "Cebpb", "Fosb", "Rasd1", "Gadd45g", "Depp1", "Serpine1", "Errfi1",
  "Nr4a1", "Sgk1", "Rgs1", "Egr1", "Jun"
)

canonical_panels <- list(
  Hair_sensory = c("Pou4f3", "Myo7a", "Otof", "Slc26a5", "Ocm", "Pvalb"),
  Supporting_epithelial = c("Sox2", "Gjb2", "Gjb6", "Krt8", "Krt18", "Epcam", "Krt19"),
  Endothelial_vascular = c("Pecam1", "Cdh5", "Kdr", "Emcn", "Vwf", "Esam", "Klf2"),
  Immune_myeloid = c("Ptprc", "Lyz2", "Tyrobp", "Aif1", "C1qa", "C1qb", "C1qc", "Cd68", "Ctss", "Fcerg1"),
  Fibrocyte_mesenchymal = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Pdgfra", "Col6a1"),
  Neuronal = c("Tubb3", "Snap25", "Syt1", "Rbfox3", "Nefh", "Prph"),
  Glial = c("Plp1", "Mpz", "Sox10", "S100b", "Fabp7")
)
