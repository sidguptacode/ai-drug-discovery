#!/usr/bin/env Rscript
# =============================================================================
# step_5.R — Export cell-type metadata for the Python pipeline
# Reads:  OUT_DIR/step4_seurat_annotated.rds
#         REQUIRED: OUT_DIR/step4_adjudication_labels.csv
#         REQUIRED: OUT_DIR/step4_adjudication_report.json
# Writes: OUT_DIR/cell_type_metadata.csv   <- input for step_6.py
#         OUT_DIR/seurat_integrated.rds
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
  library(jsonlite)
})

cfg_path <- Sys.getenv("PIPELINE_STEP_CONFIG", "config.yml")
if (grepl("\\.json$", cfg_path, ignore.case = TRUE)) {
  cfg <- jsonlite::read_json(cfg_path, simplifyVector = FALSE)
} else {
  cfg <- yaml::read_yaml(cfg_path)
}
OUT_DIR <- cfg$out_dir

cat(sprintf("====== step_5.R | Export | Seurat v%s ======\n", packageVersion("Seurat")))

step4_path <- file.path(OUT_DIR, "step4_seurat_annotated.rds")
if (!file.exists(step4_path)) stop("Missing: ", step4_path, "\n  Run step_4.R first.")
adj_path <- file.path(OUT_DIR, "step4_adjudication_labels.csv")
rep_path <- file.path(OUT_DIR, "step4_adjudication_report.json")
if (!file.exists(adj_path) || !file.exists(rep_path)) {
  missing <- c()
  if (!file.exists(adj_path)) missing <- c(missing, basename(adj_path))
  if (!file.exists(rep_path)) missing <- c(missing, basename(rep_path))
  stop(
    "Step 4 final outputs not completed: missing ",
    paste(missing, collapse = " and "),
    " in OUT_DIR. After step 4, the reflecting agent must write adjudication labels and ",
    "report (see pipeline skill) before step 5.\n  OUT_DIR: ",
    OUT_DIR,
    call. = FALSE
  )
}

cat("  Loading step4_seurat_annotated.rds...\n")
seurat_int <- readRDS(step4_path)
cat(sprintf("  %d spots | %d communities | %d samples\n",
            ncol(seurat_int),
            length(unique(seurat_int$community)),
            length(unique(seurat_int$sample))))

cat("  Applying final labels from step4_adjudication_labels.csv (all clusters)...\n")
adj <- read.csv(adj_path, stringsAsFactors = FALSE, check.names = FALSE)
req <- c("cluster", "adjudication_label")
miss <- setdiff(req, names(adj))
if (length(miss) > 0) {
  stop("step4_adjudication_labels.csv is missing required columns: ",
       paste(miss, collapse = ", "))
}
comm_char <- as.character(seurat_int$community)
n_applied <- 0L
for (i in seq_len(nrow(adj))) {
  cl <- trimws(as.character(adj$cluster[i]))
  if (!nzchar(cl)) next
  lab <- trimws(as.character(adj$adjudication_label[i]))
  if (!nzchar(lab)) {
    stop("Empty adjudication_label for cluster ", cl,
         " — every cluster must have a final label before step 5.")
  }
  hit <- comm_char == cl
  if (!any(hit)) {
    warning("No spots found for adjudication cluster ", cl, " — skipping.")
    next
  }
  seurat_int$cell_type_label[hit] <- lab
  n_applied <- n_applied + sum(hit)
}
communities <- sort(unique(as.character(seurat_int$community)))
clusters_in_csv <- unique(trimws(as.character(adj$cluster)))
clusters_in_csv <- clusters_in_csv[nzchar(clusters_in_csv)]
missing_cov <- setdiff(communities, clusters_in_csv)
if (length(missing_cov) > 0) {
  stop("step4_adjudication_labels.csv is missing rows for communities: ",
       paste(missing_cov, collapse = ", "), call. = FALSE)
}
cat(sprintf("  Set cell_type_label for %d spot(s) from adjudication CSV.\n", n_applied))

export_df <- data.frame(
  barcode         = colnames(seurat_int),
  sample          = seurat_int$sample,
  community       = seurat_int$community,
  cell_type_label = seurat_int$cell_type_label
)
write.csv(export_df, file.path(OUT_DIR, "cell_type_metadata.csv"), row.names = FALSE)
saveRDS(seurat_int, file.path(OUT_DIR, "seurat_integrated.rds"))

cat(sprintf("  Spots: %d | Communities: %d | Samples: %d\n",
            nrow(export_df),
            length(unique(seurat_int$community)),
            length(unique(seurat_int$sample))))
cat("  Saved: cell_type_metadata.csv\n")
cat("  Saved: seurat_integrated.rds\n")
cat("\n====== step_5.R COMPLETE - proceed to step_6.py ======\n")
