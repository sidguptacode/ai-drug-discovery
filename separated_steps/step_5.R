#!/usr/bin/env Rscript
# =============================================================================
# step_5.R — Export cell-type metadata for the Python pipeline
# Reads:  OUT_DIR/step4_seurat_annotated.rds
#         Optional: OUT_DIR/step4_adjudication_labels.csv (Step 4b overrides)
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
SAMPLES <- if (is.list(cfg$samples)) unlist(cfg$samples) else cfg$samples

cat(sprintf("====== step_5.R | Export | Seurat v%s ======\n", packageVersion("Seurat")))

step4_path <- file.path(OUT_DIR, "step4_seurat_annotated.rds")
if (!file.exists(step4_path)) stop("Missing: ", step4_path, "\n  Run step_4.R first.")
cat("  Loading step4_seurat_annotated.rds...\n")
seurat_int <- readRDS(step4_path)
cat(sprintf("  %d spots | %d communities | %d samples\n",
            ncol(seurat_int),
            length(unique(seurat_int$community)),
            length(unique(seurat_int$sample))))

# Optional: apply Step 4b adjudication overrides (step4_adjudication_labels.csv)
adj_path <- file.path(OUT_DIR, "step4_adjudication_labels.csv")
if (file.exists(adj_path)) {
  cat("  Found step4_adjudication_labels.csv — applying overrides where override_applied is true...\n")
  adj <- read.csv(adj_path, stringsAsFactors = FALSE, check.names = FALSE)
  req <- c("cluster", "adjudication_label", "override_applied")
  miss <- setdiff(req, names(adj))
  if (length(miss) > 0) {
    stop("step4_adjudication_labels.csv is missing required columns: ",
         paste(miss, collapse = ", "))
  }
  parse_override <- function(x) {
    if (is.logical(x)) return(x)
    v <- toupper(trimws(as.character(x)))
    v %in% c("TRUE", "T", "1", "YES", "Y")
  }
  comm_char <- as.character(seurat_int$community)
  n_override <- 0L
  for (i in seq_len(nrow(adj))) {
    if (!parse_override(adj$override_applied[i])) next
    cl <- trimws(as.character(adj$cluster[i]))
    if (!nzchar(cl)) next
    lab <- trimws(as.character(adj$adjudication_label[i]))
    if (!nzchar(lab)) {
      warning("override_applied is true but adjudication_label is empty for cluster ",
              cl, " — skipping.")
      next
    }
    hit <- comm_char == cl
    if (!any(hit)) {
      warning("No spots found for adjudication cluster ", cl, " — skipping.")
      next
    }
    seurat_int$cell_type_label[hit] <- lab
    n_override <- n_override + sum(hit)
  }
  cat(sprintf("  Adjudication: updated %d spot label(s) across override_applied clusters.\n",
              n_override))
} else {
  cat("  No step4_adjudication_labels.csv — exporting step 4 labels as-is.\n")
}

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
