suppressPackageStartupMessages(library(data.table))

project_dir <- Sys.getenv("GSE252285_PROJECT_ROOT", unset = ".")
branch <- "no_doublet_robust_qc_tsne_sensitivity"
analysis <- "condition_descriptive_top20_resolved_geo_mapping"
out_path <- function(section, ...) file.path(project_dir, section, branch, analysis, ...)

before_path <- out_path("00_manifest", "INPUT_OBJECT_SHA256_BEFORE.tsv")
after_path <- out_path("00_manifest", "INPUT_OBJECT_SHA256_AFTER.tsv")
validation_path <- out_path("09_qc", "CONDITION_DESCRIPTIVE_VALIDATION.tsv")
mapping_path <- out_path("00_manifest", "GEO_CONDITION_MAPPING_AUDIT.tsv")

before <- readLines(before_path, warn = FALSE)
after <- readLines(after_path, warn = FALSE)
validation <- fread(validation_path)
mapping <- fread(mapping_path)
checks <- data.table(
  check = c("input_object_sha256_unchanged", "all_pipeline_validations_pass",
            "mapping_resolved_yes", "mapping_complete_three_rows",
            "no_rds_written_in_analysis_outputs"),
  pass = c(
    identical(before, after),
    nrow(validation) > 0L && all(validation$pass),
    nrow(mapping) == 3L && all(mapping$CONDITION_MAPPING_RESOLVED == "YES"),
    nrow(mapping) == 3L && all(mapping$mapping_complete),
    !any(grepl("\\.rds$", unlist(lapply(
      c("00_manifest", "04_tables", "05_figures", "07_logs", "08_methods", "09_qc", "10_release"),
      function(section) list.files(out_path(section), recursive = TRUE, full.names = TRUE)
    )), ignore.case = TRUE))
  )
)
fwrite(checks, out_path("09_qc", "RELEASE_CONSISTENCY_FINAL.tsv"), sep = "\t")
if (any(!checks$pass)) stop("Final release consistency check failed")

roots <- vapply(c("00_manifest", "04_tables", "05_figures", "07_logs",
                  "08_methods", "09_qc", "10_release"),
                function(section) out_path(section), character(1))
files <- sort(unique(c(
  unlist(lapply(roots, function(root) list.files(root, recursive = TRUE, full.names = TRUE))),
  file.path(project_dir, "03_scripts", branch, "17_condition_descriptive_top20.R"),
  file.path(project_dir, "03_scripts", branch, "18_finalize_condition_release.R")
)))
files <- files[file.exists(files) & !file.info(files)$isdir]
manifest_path <- out_path("00_manifest", "OUTPUT_MANIFEST_FINAL.tsv")
files <- setdiff(files, manifest_path)
sha256 <- vapply(files, function(path) {
  sub(" .*", "", system2("sha256sum", shQuote(path), stdout = TRUE))
}, character(1))
manifest <- data.table(
  relative_path = sub(paste0("^", project_dir, "/"), "", files),
  bytes = as.numeric(file.info(files)$size),
  sha256 = sha256
)
fwrite(manifest, manifest_path, sep = "\t")
cat("FINAL_RELEASE_CONSISTENCY=PASS\n")
