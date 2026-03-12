#!/usr/bin/env Rscript
# =============================================================================
# step_5.R — Export cell-type metadata for the Python pipeline
# Reads:  OUT_DIR/step4_seurat_annotated.rds
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
