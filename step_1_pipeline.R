#!/usr/bin/env Rscript
# =============================================================================
# pipeline.R — ST LR pipeline (dataset-agnostic)
# All parameters read from config.yaml
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(hdf5r)
})

try_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Package '%s' not installed — skipping.", pkg))
    return(FALSE)
  }
  library(pkg, character.only = TRUE); TRUE
}

has_yaml     <- try_library("yaml")
has_clustree <- try_library("clustree")
has_enrichr  <- try_library("enrichR")

if (!has_yaml) stop("yaml package required. Install with: install.packages('yaml')")

# ── Load config ───────────────────────────────────────────────────────────────
cfg <- yaml::read_yaml("/scratch/baderlab/sgupta/aidrugdiscovery_st_analysis/config.yml")

DATASET   <- cfg$dataset_name
DATA_DIR  <- cfg$data_dir
OUT_DIR   <- cfg$out_dir
SAMPLES   <- cfg$samples
SPECIES   <- cfg$species

QC        <- cfg$qc
INTEG     <- cfg$integration
CLUST     <- cfg$clustering
ANNOT     <- cfg$annotation

dir.create(OUT_DIR, showWarnings = FALSE)
set.seed(42)

SEURAT_MAJOR <- as.integer(strsplit(as.character(packageVersion("Seurat")), "\\.")[[1]][1])
cat(sprintf("====== pipeline.R | %s | Seurat v%s ======\n",
            DATASET, packageVersion("Seurat")))
cat(sprintf("  Samples: %d | Species: %s\n", length(SAMPLES), SPECIES))

# =============================================================================
# STEP 1: QC
# =============================================================================
cat("\n================ STEP 1: QC ================\n")
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

  n_before <- ncol(so)
  filter_expr <- sprintf("nCount_RNA >= %d & nCount_RNA <= %d",
                         QC$ncount_min, QC$ncount_max)
  if (!is.null(QC$mt_cutoff))
    filter_expr <- paste0(filter_expr, sprintf(" & percent.mt <= %g", QC$mt_cutoff))
  keep <- with(so@meta.data, eval(parse(text = filter_expr)))
  so <- so[, which(keep)]

  cat(sprintf("  %s: %d → %d spots\n", samp, n_before, ncol(so)))
  so
}

seurat_list <- lapply(SAMPLES, load_and_qc)
names(seurat_list) <- SAMPLES

pdf(file.path(OUT_DIR, "step1_qc.pdf"), width = 14, height = 5)
for (samp in SAMPLES) {
  so  <- seurat_list[[samp]]
  df  <- so@meta.data
  p1  <- VlnPlot(so, features = c("nCount_RNA","nFeature_RNA","percent.mt"),
                 ncol = 3, pt.size = 0) &
    theme(axis.title.x = element_blank())
  p2  <- ggplot(df, aes(pxl_col, -pxl_row, colour = log10(nCount_RNA + 1))) +
    geom_point(size = 0.4) + scale_colour_viridis_c("log10(nCount)") +
    coord_fixed() + theme_bw() + labs(title = paste(samp, "— nCount"))
  p3  <- ggplot(df, aes(pxl_col, -pxl_row, colour = percent.mt)) +
    geom_point(size = 0.4) + scale_colour_viridis_c("% MT") +
    coord_fixed() + theme_bw() + labs(title = paste(samp, "— % MT"))
  print((p1 / (p2 + p3)) + plot_annotation(title = samp))
}
dev.off()
cat("  Saved: outputs/step1_qc.pdf\n")

# =============================================================================
# STEP 2: Integration
# =============================================================================
cat("\n================ STEP 2: Integration ================\n")
cat(sprintf("  nfeatures=%d | n_pcs=%d | n_dims=%d\n",
            INTEG$n_features, INTEG$n_pcs, INTEG$n_dims))
options(future.globals.maxSize = 2000 * 1024^2)  # 2 GB for SCTransform parallelism

seurat_list <- lapply(seurat_list, function(so) {
  SCTransform(so, vst.flavor = "v2", verbose = FALSE,
              variable.features.n = INTEG$n_features,
              return.only.var.genes = FALSE)
})
features    <- SelectIntegrationFeatures(object.list = seurat_list,
                                          nfeatures = INTEG$n_features)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
anchors     <- FindIntegrationAnchors(object.list = seurat_list,
                                      normalization.method = "SCTransform",
                                      anchor.features = features, verbose = FALSE)
seurat_int  <- IntegrateData(anchorset = anchors,
                              normalization.method = "SCTransform", verbose = FALSE)
DefaultAssay(seurat_int) <- "integrated"
seurat_int  <- RunPCA(seurat_int, npcs = INTEG$n_pcs, verbose = FALSE)
seurat_int  <- RunUMAP(seurat_int, dims = 1:INTEG$n_dims, verbose = FALSE, seed.use = 42)

pdf(file.path(OUT_DIR, "step2_integration.pdf"), width = 10, height = 7)
print(ElbowPlot(seurat_int, ndims = INTEG$n_pcs) + ggtitle("PCA Elbow"))
print(DimPlot(seurat_int, reduction = "umap", group.by = "sample", pt.size = 0.3) +
        ggtitle("UMAP — by sample (batch check)") + theme_bw())
dev.off()
cat("  Saved: outputs/step2_integration.pdf\n")

# =============================================================================
# STEP 3: Clustering
# =============================================================================
cat("\n================ STEP 3: Clustering ================\n")
cat(sprintf("  Resolutions tested: %s | Chosen: %g\n",
            paste(CLUST$resolutions, collapse=", "), CLUST$chosen_resolution))

seurat_int <- FindNeighbors(seurat_int, dims = 1:INTEG$n_dims, verbose = FALSE)
for (res in CLUST$resolutions) {
  seurat_int <- FindClusters(seurat_int, resolution = res, verbose = FALSE)
}

if (has_clustree) {
  pdf(file.path(OUT_DIR, "step3_clustree.pdf"), width = 12, height = 10)
  print(clustree(seurat_int, prefix = "integrated_snn_res.") +
          ggtitle("Clustree — inspect to confirm chosen_resolution in config.yaml"))
  dev.off()
  cat("  Saved: outputs/step3_clustree.pdf\n")
  cat("  >>> If chosen_resolution looks wrong, update config.yaml and rerun.\n")
}

res_col <- paste0("integrated_snn_res.", CLUST$chosen_resolution)
Idents(seurat_int) <- res_col
seurat_int$community <- paste0("C", as.integer(seurat_int[[res_col]][[1]]))

n_comm <- length(unique(seurat_int$community))
cat(sprintf("  Resolution %g → %d communities\n", CLUST$chosen_resolution, n_comm))
print(sort(table(seurat_int$community), decreasing = TRUE))

prop_df <- as.data.frame(
  table(Community = seurat_int$community, Sample = seurat_int$sample)
) %>% group_by(Sample) %>% mutate(Proportion = Freq / sum(Freq)) %>% ungroup()

pdf(file.path(OUT_DIR, "step3_clusters.pdf"), width = 14, height = 7)
print(DimPlot(seurat_int, group.by = "community", label = TRUE,
              pt.size = 0.3, repel = TRUE) +
        ggtitle(sprintf("UMAP — communities (res = %g)", CLUST$chosen_resolution)) +
        theme_bw())
print(ggplot(prop_df, aes(x = Sample, y = Proportion, fill = Community)) +
        geom_bar(stat = "identity") + theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Community proportions per sample"))
dev.off()
cat("  Saved: outputs/step3_clusters.pdf\n")

# =============================================================================
# STEP 4: Marker genes + EnrichR annotation
# =============================================================================
cat("\n================ STEP 4: Marker genes & annotation ================\n")
cat(sprintf("  FindAllMarkers: min_pct=%g, logfc>=%g | EnrichR top %d genes\n",
            ANNOT$min_pct, ANNOT$logfc_threshold, ANNOT$n_marker_genes))

DefaultAssay(seurat_int) <- "SCT"
Idents(seurat_int) <- "community"
if (SEURAT_MAJOR >= 5) {
  tryCatch(seurat_int <- PrepSCTFindMarkers(seurat_int, verbose = FALSE),
           error = function(e) message("  [WARN] PrepSCTFindMarkers: ", e$message))
}

all_markers <- FindAllMarkers(seurat_int,
                               only.pos        = TRUE,
                               min.pct         = ANNOT$min_pct,
                               logfc.threshold = ANNOT$logfc_threshold,
                               test.use        = "wilcox",
                               verbose         = FALSE)
sig_markers <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC))
write.csv(sig_markers, file.path(OUT_DIR, "step4_markers.csv"), row.names = FALSE)
cat(sprintf("  Found %d significant markers across %d communities.\n",
            nrow(sig_markers), length(unique(sig_markers$cluster))))

all_communities  <- sort(unique(seurat_int$community))
cell_type_labels <- setNames(all_communities, all_communities)

if (has_enrichr) {
  cat(sprintf("  Querying EnrichR: %s\n", paste(ANNOT$enrichr_dbs, collapse=", ")))
  tryCatch({
    setEnrichrSite("Enrichr")
    for (comm in sort(unique(sig_markers$cluster))) {
      top_genes <- sig_markers %>%
        filter(cluster == comm) %>%
        arrange(desc(avg_log2FC)) %>%
        slice_head(n = ANNOT$n_marker_genes) %>%
        pull(gene)
      if (length(top_genes) < 5) next
      tryCatch({
        res     <- enrichr(top_genes, ANNOT$enrichr_dbs)
        top_hit <- res[[ANNOT$primary_db]] %>%
          arrange(Adjusted.P.value) %>% slice_head(n = 1) %>% pull(Term)
        comm_key <- paste0("C", comm)
        cell_type_labels[comm_key] <- if (length(top_hit) > 0 && !is.na(top_hit))
          top_hit else comm_key
        for (db in names(res)) {
          write.csv(res[[db]],
                    file.path(OUT_DIR, sprintf("step4_enrichr_%s_%s.csv",
                                               comm_key,
                                               gsub("[^A-Za-z0-9]","_",db))),
                    row.names = FALSE)
        }
      }, error = function(e) cat(sprintf("    [WARN] EnrichR C%s: %s\n", comm, e$message)))
    }
    cat("\n  Labels (review step4_enrichr_*.csv before proceeding):\n")
    for (comm in names(cell_type_labels))
      cat(sprintf("    %-6s → %s\n", comm, cell_type_labels[comm]))
  }, error = function(e) {
    cat(sprintf("  [WARN] EnrichR failed: %s — using community IDs.\n", e$message))
  })
} else {
  cat("  enrichR not installed — labels default to community IDs.\n")
}

seurat_int$cell_type_label <- cell_type_labels[seurat_int$community]

top5       <- sig_markers %>% group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>% ungroup()
plot_genes <- unique(top5$gene)

pdf(file.path(OUT_DIR, "step4_annotation.pdf"),
    width = max(16, length(plot_genes) * 0.4), height = 8)
print(DotPlot(seurat_int, features = plot_genes, group.by = "community", dot.min = 0.05) +
        RotatedAxis() + scale_colour_gradient2(low="blue",mid="white",high="red") +
        labs(title="Top 5 marker genes per community", x=NULL, y="Community") +
        theme(axis.text.x = element_text(size=6, face="italic")))
print(DimPlot(seurat_int, reduction="umap", group.by="cell_type_label",
              label=TRUE, label.size=3, pt.size=0.3, repel=TRUE) +
        ggtitle("UMAP — cell type labels") + theme_bw() +
        theme(legend.text = element_text(size=7)))
dev.off()
cat("  Saved: outputs/step4_annotation.pdf\n")

# =============================================================================
# STEP 5: Export
# =============================================================================
cat("\n================ STEP 5: Export ================\n")

export_df <- data.frame(
  barcode         = colnames(seurat_int),
  sample          = seurat_int$sample,
  community       = seurat_int$community,
  cell_type_label = seurat_int$cell_type_label
)
write.csv(export_df, file.path(OUT_DIR, "cell_type_metadata.csv"), row.names = FALSE)
saveRDS(seurat_int,  file.path(OUT_DIR, "seurat_integrated.rds"))

cat(sprintf("  Spots: %d | Communities: %d | Samples: %d\n",
            nrow(export_df), n_comm, length(SAMPLES)))
cat("  Saved: outputs/cell_type_metadata.csv, outputs/seurat_integrated.rds\n")
cat("\n====== pipeline.R COMPLETE — proceed to pipeline.py ======\n")