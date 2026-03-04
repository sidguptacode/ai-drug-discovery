#!/usr/bin/env Rscript
# =============================================================================
# step_1.R — QC
# Reads:  DATA_DIR/<sample>/<sample>_filtered_feature_bc_matrix.h5
#         DATA_DIR/<sample>/<sample>_tissue_positions_list.csv
# Writes: OUT_DIR/step1_seurat_list.rds
#         OUT_DIR/step1_qc.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(yaml)
})

options(repos = c(CRAN = "https://cloud.r-project.org/"))

cfg      <- yaml::read_yaml("/scratch/baderlab/sgupta/ai-drug-discovery/config.yml")
DATA_DIR <- cfg$data_dir
OUT_DIR  <- cfg$out_dir
SAMPLES  <- cfg$samples
QC       <- cfg$qc

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

cat(sprintf("====== step_1.R | QC | Seurat v%s ======\n", packageVersion("Seurat")))
cat(sprintf("  nCount filter: [%d, %d]", QC$ncount_min, QC$ncount_max))
if (!is.null(QC$mt_cutoff)) cat(sprintf(" | MT cutoff: %.0f%%", QC$mt_cutoff))
cat("\n")

load_and_qc <- function(samp) {
  h5_path    <- file.path(DATA_DIR, samp, paste0(samp, "_filtered_feature_bc_matrix.h5"))
  coord_path <- file.path(DATA_DIR, samp, paste0(samp, "_tissue_positions_list.csv"))

  counts <- Read10X_h5(h5_path, use.names = TRUE, unique.features = TRUE)
  so     <- CreateSeuratObject(counts = counts, project = samp,
                               min.cells = 0, min.features = 0)
  so$sample <- samp
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

  coords <- tryCatch({
    df <- read.csv(coord_path, header = FALSE,
                   col.names = c("barcode","in_tissue","array_row",
                                 "array_col","pxl_row","pxl_col"))
    df[df$in_tissue == 1, ]
  }, error = function(e) {
    df <- read.csv(coord_path)
    if ("in_tissue" %in% colnames(df)) df[df$in_tissue == 1, ] else df
  })
  rownames(coords) <- coords$barcode
  shared <- intersect(colnames(so), rownames(coords))
  so <- so[, shared]
  so$pxl_row <- as.numeric(coords[shared, "pxl_row"])
  so$pxl_col <- as.numeric(coords[shared, "pxl_col"])

  n_before    <- ncol(so)
  filter_expr <- sprintf("nCount_RNA >= %d & nCount_RNA <= %d",
                         QC$ncount_min, QC$ncount_max)
  if (!is.null(QC$mt_cutoff))
    filter_expr <- paste0(filter_expr, sprintf(" & percent.mt <= %g", QC$mt_cutoff))
  keep <- with(so@meta.data, eval(parse(text = filter_expr)))
  so   <- so[, which(keep)]

  cat(sprintf("  %s: %d -> %d spots\n", samp, n_before, ncol(so)))
  so
}

seurat_list        <- lapply(SAMPLES, load_and_qc)
names(seurat_list) <- SAMPLES

saveRDS(seurat_list, file.path(OUT_DIR, "step1_seurat_list.rds"))
cat("  Saved: step1_seurat_list.rds\n")

pdf(file.path(OUT_DIR, "step1_qc.pdf"), width = 14, height = 5)
for (samp in SAMPLES) {
  so <- seurat_list[[samp]]
  df <- so@meta.data
  p1 <- VlnPlot(so, features = c("nCount_RNA","nFeature_RNA","percent.mt"),
                ncol = 3, pt.size = 0) &
    theme(axis.title.x = element_blank())
  p2 <- ggplot(df, aes(pxl_col, -pxl_row, colour = log10(nCount_RNA + 1))) +
    geom_point(size = 0.4) + scale_colour_viridis_c("log10(nCount)") +
    coord_fixed() + theme_bw() + labs(title = paste(samp, "- nCount"))
  p3 <- ggplot(df, aes(pxl_col, -pxl_row, colour = percent.mt)) +
    geom_point(size = 0.4) + scale_colour_viridis_c("% MT") +
    coord_fixed() + theme_bw() + labs(title = paste(samp, "- % MT"))
  print((p1 / (p2 + p3)) + plot_annotation(title = samp))
}
dev.off()
cat("  Saved: step1_qc.pdf\n")

cat(sprintf("\n====== step_1.R COMPLETE | %d samples ======\n", length(SAMPLES)))
