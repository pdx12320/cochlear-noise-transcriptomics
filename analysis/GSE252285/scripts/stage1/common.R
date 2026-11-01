options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

project_dir <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
p <- function(...) file.path(project_dir, ...)

required_dirs <- c(
  "00_manifest", "01_input", "02_processed", "03_scripts", "04_tables",
  "05_figures", "06_objects", "07_logs", "08_methods", "09_qc", "10_release"
)
invisible(lapply(file.path(project_dir, required_dirs), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

seed <- 42L
set.seed(seed)

write_tsv <- function(x, path) {
  data.table::fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

save_dual <- function(plot, stem, width = 8, height = 6, dpi = 220) {
  pdf_path <- p("05_figures", paste0(stem, ".pdf"))
  png_path <- p("05_figures", paste0(stem, ".png"))
  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height,
                  units = "in", limitsize = FALSE)
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height,
                  units = "in", dpi = dpi, device = ragg::agg_png,
                  limitsize = FALSE)
  invisible(c(pdf_path, png_path))
}

barcode_suffix <- function(x) {
  has_suffix <- grepl("-[0-9]+$", x)
  out <- rep("NO_SUFFIX", length(x))
  out[has_suffix] <- sub("^.*(-[0-9]+)$", "\\1", x[has_suffix])
  out
}

robust_upper <- function(x, q = 0.995, nmads = 6) {
  x <- x[is.finite(x)]
  qlim <- unname(stats::quantile(x, q, na.rm = TRUE, type = 8))
  med <- stats::median(x, na.rm = TRUE)
  madv <- stats::mad(x, center = med, constant = 1.4826, na.rm = TRUE)
  madlim <- if (is.finite(madv) && madv > 0) med + nmads * madv else Inf
  limit <- min(qlim, madlim)
  if (!is.finite(limit)) limit <- qlim
  list(limit = limit, median = med, mad = madv, q995 = qlim, mad_limit = madlim)
}

theme_sc <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = element_blank(),
                   strip.background = element_rect(fill = "grey95"))
}

stage_message <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}
