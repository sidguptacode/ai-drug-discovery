#!/usr/bin/env Rscript
# =============================================================================
# step_1.R — QC
# Reads:  (h5)  DATA_DIR/<sample>/<sample>_filtered_feature_bc_matrix.h5
#               DATA_DIR/<sample>/<sample>_tissue_positions_list.csv
#         (mtx, flat)   DATA_DIR/<sample>_{matrix.mtx,barcodes.tsv,features.tsv}
#                       DATA_DIR/<sample>_tissue_positions_list.csv
#         (mtx, subdir) DATA_DIR/<sample>_{matrix.mtx,barcodes.tsv,features.tsv}
#                       DATA_DIR/<sample>_spatial/spatial/tissue_positions_list.csv
# Writes: OUT_DIR/step1_seurat_list.rds
#         OUT_DIR/step1_qc.pdf
#
# Usage:  Rscript step_1.R [config.yml]
#         Falls back to config_dipg.yml if no argument given.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(yaml)
})

options(repos = c(CRAN = "https://cloud.r-project.org/"))

args     <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else
              "/scratch/baderlab/sgupta/workflows_march/mar9_ai_drug_discovery/config_dipg.yml"

cfg      <- yaml::read_yaml(cfg_path)
DATA_DIR <- cfg$data_dir
OUT_DIR  <- cfg$out_dir
SAMPLES  <- cfg$samples
QC       <- cfg$qc

DATA_FORMAT    <- if (!is.null(cfg$data_format))    cfg$data_format    else "h5"
SPATIAL_LAYOUT <- if (!is.null(cfg$spatial_layout)) cfg$spatial_layout else "subdir"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

cat(sprintf("====== step_1.R | QC | Seurat v%s ======\n", packageVersion("Seurat")))
cat(sprintf("  format=%s | spatial=%s\n", DATA_FORMAT, SPATIAL_LAYOUT))
cat(sprintf("  nCount filter: [%d, %d]", QC$ncount_min, QC$ncount_max))
if (!is.null(QC$mt_cutoff)) cat(sprintf(" | MT cutoff: %.0f%%", QC$mt_cutoff))
cat("\n")

# Returns named list: coord_path, sf_path
get_spatial_paths <- function(samp) {
  if (DATA_FORMAT == "h5") {
    list(
      coord_path = file.path(DATA_DIR, samp, paste0(samp, "_tissue_positions_list.csv")),
      sf_path    = file.path(DATA_DIR, samp, paste0(samp, "_scalefactors_json.json"))
    )
  } else if (SPATIAL_LAYOUT == "subdir") {
    list(
      coord_path = file.path(DATA_DIR, paste0(samp, "_spatial"), "spatial",
                             "tissue_positions_list.csv"),
      sf_path    = file.path(DATA_DIR, paste0(samp, "_spatial"), "spatial",
                             "scalefactors_json.json")
    )
  } else {
    # flat: all files alongside MTX with same prefix
    list(
      coord_path = file.path(DATA_DIR, paste0(samp, "_tissue_positions_list.csv")),
      sf_path    = file.path(DATA_DIR, paste0(samp, "_scalefactors_json.json"))
    )
  }
}

read_tissue_positions <- function(coord_path) {
  tryCatch({
    df <- read.csv(coord_path, header = FALSE,
                   col.names = c("barcode","in_tissue","array_row",
                                 "array_col","pxl_row","pxl_col"))
    if (is.character(df$in_tissue)) stop("header row detected")
    df[df$in_tissue == 1, ]
  }, error = function(e) {
    df <- read.csv(coord_path)
    colnames(df) <- tolower(colnames(df))
    if ("in_tissue" %in% colnames(df)) df[df$in_tissue == 1, ] else df
  })
}

load_and_qc <- function(samp) {
  sp <- get_spatial_paths(samp)

  if (DATA_FORMAT == "h5") {
    h5_path <- file.path(DATA_DIR, samp, paste0(samp, "_filtered_feature_bc_matrix.h5"))
    counts  <- Read10X_h5(h5_path, use.names = TRUE, unique.features = TRUE)
  } else {
    mtx_path  <- file.path(DATA_DIR, paste0(samp, "_matrix.mtx"))
    bc_path   <- file.path(DATA_DIR, paste0(samp, "_barcodes.tsv"))
    feat_path <- file.path(DATA_DIR, paste0(samp, "_features.tsv"))
    counts    <- ReadMtx(mtx            = mtx_path,
                         cells          = bc_path,
                         features       = feat_path,
                         feature.column = 2,
                         unique.features = TRUE)
  }

  so <- CreateSeuratObject(counts = counts, project = samp,
                           min.cells = 0, min.features = 0)
  so$sample <- samp
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

  coords <- read_tissue_positions(sp$coord_path)
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

  cat(sprintf("  %s: %d -> %d spots (after in_tissue + QC filters)\n",
              samp, n_before, ncol(so)))
  so
}

step1_ckpt <- file.path(OUT_DIR, "step1_seurat_list.rds")
if (file.exists(step1_ckpt)) {
  cat(sprintf("  [RESUME] Loading existing %s — skipping QC re-run.\n", step1_ckpt))
  seurat_list <- readRDS(step1_ckpt)
  names(seurat_list) <- SAMPLES
} else {
  seurat_list        <- lapply(SAMPLES, load_and_qc)
  names(seurat_list) <- SAMPLES
  saveRDS(seurat_list, step1_ckpt)
  cat("  Saved: step1_seurat_list.rds\n")
}

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
