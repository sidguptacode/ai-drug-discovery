#!/usr/bin/env Rscript
# =============================================================================
# pipeline_cleaned.R — Steps 1–5 of the ST LR pipeline
#
# Step 1 : QC                    (load H5 + spatial coords, filter spots)
# Step 2 : Integration           (SCTransform per sample, RPCA anchor integration)
# Step 3 : Clustering            (FindClusters, clustree, community assignment)
# Step 4 : Annotation            (FindAllMarkers, EnrichR, CNS-preference labelling)
# Step 5 : Export                (cell_type_metadata.csv → input for pipeline.py)
#
# All parameters read from config.yml
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(yaml)
})

try_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("  Package '%s' not installed — skipping.", pkg))
    return(FALSE)
  }
  library(pkg, character.only = TRUE); TRUE
}

options(repos = c(CRAN = "https://cloud.r-project.org/"))
options(future.globals.maxSize = Inf)

if (!requireNamespace("harmony", quietly = TRUE)) {
  message("Installing harmony from CRAN...")
  install.packages("harmony")
}

has_harmony  <- try_library("harmony")
has_clustree <- try_library("clustree")
has_enrichr  <- try_library("enrichR")

if (!has_harmony) stop("harmony package required. Install with: install.packages('harmony')")

# Timestamped progress helper (useful in SLURM logs)
ts <- function(msg) {
  cat(sprintf("  [%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  invisible(Sys.time())
}

# ── Load config ────────────────────────────────────────────────────────────────
cfg     <- yaml::read_yaml("/scratch/baderlab/sgupta/ai-drug-discovery/config.yml")
DATASET <- cfg$dataset_name
DATA_DIR <- cfg$data_dir
OUT_DIR  <- cfg$out_dir
SAMPLES  <- cfg$samples
QC       <- cfg$qc
INTEG    <- cfg$integration
CLUST    <- cfg$clustering
ANNOT    <- cfg$annotation

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

cat(sprintf("====== pipeline_cleaned.R | %s | Seurat v%s ======\n",
            DATASET, packageVersion("Seurat")))
cat(sprintf("  Samples: %d | Out: %s\n", length(SAMPLES), OUT_DIR))

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
  so   <- so[, which(keep)]

  cat(sprintf("  %s: %d → %d spots\n", samp, n_before, ncol(so)))
  so
}

seurat_list <- lapply(SAMPLES, load_and_qc)
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
    coord_fixed() + theme_bw() + labs(title = paste(samp, "— nCount"))
  p3 <- ggplot(df, aes(pxl_col, -pxl_row, colour = percent.mt)) +
    geom_point(size = 0.4) + scale_colour_viridis_c("% MT") +
    coord_fixed() + theme_bw() + labs(title = paste(samp, "— % MT"))
  print((p1 / (p2 + p3)) + plot_annotation(title = samp))
}
dev.off()
cat("  Saved: step1_qc.pdf\n")

# =============================================================================
# STEP 2: Integration (SCTransform + RPCA)
# =============================================================================
cat("\n================ STEP 2: Integration (SCTransform + RPCA) ================\n")
cat(sprintf("  nfeatures=%d | n_pcs=%d | n_dims=%d\n",
            INTEG$n_features, INTEG$n_pcs, INTEG$n_dims))

# ── Part 1: SCTransform + PrepSCTIntegration + PCA per sample ─────────────────
n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_workers) || n_workers < 1L)
  n_workers <- max(1L, parallel::detectCores(logical = FALSE) %/% 2L)
cat(sprintf("  Parallel workers: %d\n", n_workers))
if (.Platform$OS.type == "unix") {
  future::plan(future::multicore, workers = n_workers)
} else {
  future::plan(future::multisession, workers = n_workers)
}

sct_flavor <- if (requireNamespace("glmGamPoi", quietly = TRUE)) {
  cat("  SCTransform flavor: v2\n"); "v2"
} else {
  cat("  SCTransform flavor: v1\n"); "v1"
}

cat(sprintf("  Running SCTransform on %d samples...\n", length(seurat_list)))
seurat_list <- lapply(seq_along(seurat_list), function(i) {
  samp <- names(seurat_list)[i]
  so   <- seurat_list[[i]]
  cat(sprintf("    [%d/%d] SCTransform: %s (%d spots)...\n",
              i, length(seurat_list), samp, ncol(so)))
  SCTransform(so,
              vst.flavor            = sct_flavor,
              variable.features.n   = INTEG$n_features,
              ncells                = min(ncol(so), 5000L),
              vars.to.regress       = NULL,
              return.only.var.genes = FALSE,
              verbose               = FALSE)
})
names(seurat_list) <- SAMPLES

cat(sprintf("  Selecting %d integration features...\n", INTEG$n_features))
integration_features <- SelectIntegrationFeatures(
  object.list = seurat_list,
  nfeatures   = INTEG$n_features
)

cat("  Running PrepSCTIntegration...\n")
seurat_list <- PrepSCTIntegration(
  object.list     = seurat_list,
  anchor.features = integration_features
)

cat(sprintf("  Running PCA on each sample (npcs=%d)...\n", INTEG$n_dims))
seurat_list <- lapply(seurat_list, function(so)
  RunPCA(so, npcs = INTEG$n_dims, verbose = FALSE))
names(seurat_list) <- SAMPLES

# ── Part 2: FindIntegrationAnchors + IntegrateData + UMAP + Neighbors ─────────
if (requireNamespace("pbapply", quietly = TRUE))
  pbapply::pboptions(type = "txt", char = "=")

t0 <- ts("Starting FindIntegrationAnchors (RPCA, SCT)...")
anchors <- FindIntegrationAnchors(
  object.list          = seurat_list,
  normalization.method = "SCT",
  anchor.features      = integration_features,
  reduction            = "rpca",
  dims                 = 1:INTEG$n_dims,
  k.filter             = NA,
  verbose              = TRUE
)
ts(sprintf("FindIntegrationAnchors done. Elapsed: %.1f min",
           as.numeric(difftime(Sys.time(), t0, units = "mins"))))

t0 <- ts("Starting IntegrateData...")
seurat_int <- IntegrateData(
  anchorset            = anchors,
  normalization.method = "SCT",
  dims                 = 1:INTEG$n_dims,
  verbose              = TRUE
)
ts(sprintf("IntegrateData done. Elapsed: %.1f min",
           as.numeric(difftime(Sys.time(), t0, units = "mins"))))

DefaultAssay(seurat_int) <- "integrated"
cat(sprintf("  Running PCA (npcs=%d)...\n", INTEG$n_pcs))
seurat_int <- RunPCA(seurat_int, npcs = INTEG$n_pcs, verbose = FALSE)

cat(sprintf("  Running UMAP (dims=1:%d)...\n", INTEG$n_dims))
seurat_int <- RunUMAP(seurat_int, dims = 1:INTEG$n_dims, verbose = FALSE, seed.use = 42)

cat(sprintf("  Building neighbour graph (dims=1:%d)...\n", INTEG$n_dims))
seurat_int <- FindNeighbors(seurat_int, dims = 1:INTEG$n_dims, verbose = FALSE)

pdf(file.path(OUT_DIR, "step2_integration.pdf"), width = 10, height = 7)
tryCatch(
  print(ElbowPlot(seurat_int, ndims = INTEG$n_pcs) +
          ggtitle("PCA Elbow (SCTransform + RPCA)")),
  error = function(e) message("  [WARN] ElbowPlot: ", e$message)
)
tryCatch(
  print(DimPlot(seurat_int, reduction = "umap", group.by = "sample", pt.size = 0.3) +
          ggtitle("UMAP — by sample (batch check)") + theme_bw()),
  error = function(e) message("  [WARN] DimPlot: ", e$message)
)
dev.off()
cat("  Saved: step2_integration.pdf\n")

saveRDS(seurat_int, file.path(OUT_DIR, "step2_seurat_integrated.rds"))
cat("  Saved: step2_seurat_integrated.rds\n")

# =============================================================================
# STEP 3: Clustering
# =============================================================================
cat("\n================ STEP 3: Clustering ================\n")
cat(sprintf("  Resolutions tested: %s | Chosen: %g\n",
            paste(CLUST$resolutions, collapse=", "), CLUST$chosen_resolution))

# Detect the SNN graph — "integrated_snn" for SCTransform, "RNA_snn" for Harmony
snn_graph <- grep("_snn$", names(seurat_int@graphs), value = TRUE)[1]
if (is.na(snn_graph)) stop("No SNN graph found. Check that FindNeighbors ran.")
cat(sprintf("  Using graph: %s\n", snn_graph))

for (res in CLUST$resolutions)
  seurat_int <- FindClusters(seurat_int, graph.name = snn_graph,
                             resolution = res, verbose = FALSE)

clustree_prefix <- paste0(snn_graph, "_res.")

if (has_clustree && length(CLUST$resolutions) >= 2) {
  pdf(file.path(OUT_DIR, "step3_clustree.pdf"), width = 12, height = 10)
  print(clustree(seurat_int, prefix = clustree_prefix) +
          ggtitle("Clustree — inspect to confirm chosen_resolution in config.yml"))
  dev.off()
  cat("  Saved: step3_clustree.pdf\n")
  cat("  >>> If chosen_resolution looks wrong, update config.yml and rerun.\n")
}

res_col <- paste0(clustree_prefix, CLUST$chosen_resolution)
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
cat("  Saved: step3_clusters.pdf\n")

saveRDS(seurat_int, file.path(OUT_DIR, "step3_seurat_clustered.rds"))
cat("  Saved: step3_seurat_clustered.rds\n")

# =============================================================================
# STEP 4: Marker genes + EnrichR annotation
# =============================================================================
cat("\n================ STEP 4: Marker genes & annotation ================\n")
cat(sprintf("  FindAllMarkers: min_pct=%g | logfc>=%g | EnrichR top %d genes\n",
            ANNOT$min_pct, ANNOT$logfc_threshold, ANNOT$n_marker_genes))

# Switch to RNA for marker detection; JoinLayers collapses per-sample layers
# so Wilcoxon can run across all spots. NormalizeData creates the data layer
# if SCTransform was used (which never calls NormalizeData on the RNA assay).
DefaultAssay(seurat_int) <- "RNA"
cat("  Joining RNA layers...\n")
seurat_int <- JoinLayers(seurat_int)
cat("  Normalizing RNA assay...\n")
seurat_int <- NormalizeData(seurat_int, normalization.method = "LogNormalize",
                            scale.factor = 10000, verbose = FALSE)
Idents(seurat_int) <- "community"

cat("  Running FindAllMarkers...\n")
all_markers <- FindAllMarkers(seurat_int,
                              only.pos        = TRUE,
                              min.pct         = ANNOT$min_pct,
                              logfc.threshold = ANNOT$logfc_threshold,
                              test.use        = "wilcox",
                              verbose         = FALSE)

if (nrow(all_markers) == 0 || !"p_val_adj" %in% colnames(all_markers))
  stop("FindAllMarkers returned no results.")

sig_markers <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC))
write.csv(sig_markers, file.path(OUT_DIR, "step4_markers.csv"), row.names = FALSE)
cat(sprintf("  Found %d significant markers across %d communities.\n",
            nrow(sig_markers), length(unique(sig_markers$cluster))))

# ── CNS-preference label picker ───────────────────────────────────────────────
# Prefers the best-ranked EnrichR term matching a CNS tissue region OR a
# CNS-specific cell type. Falls back to the overall best hit with a
# "[non-CNS]" prefix so ambiguous labels are visible in output.
pick_cns_label <- function(db_df) {
  if (nrow(db_df) == 0) return(NA_character_)

  cns_regions <- paste(c(
    "Brain", "Cortex", "Cerebellum", "Cerebellar", "Hippocampus",
    "Spinal Cord", "Spinal", "\\bCNS\\b", "Brainstem", "Striatum",
    "Thalamus", "Hypothalamus", "Midbrain", "\\bPons\\b", "Medulla",
    "Olfactory", "Amygdala", "Prefrontal", "White Matter", "Grey Matter",
    "Gray Matter", "Frontal Lobe", "Temporal Lobe", "Ventricle",
    "Substantia Nigra", "Basal Ganglia", "Cerebrum"
  ), collapse = "|")

  cns_celltypes <- paste(c(
    "Astrocyte", "Oligodendrocyte", "\\bOPC\\b", "Oligodendrocyte Precursor",
    "Microglia", "\\bNeuron\\b", "\\bNeuronal\\b", "Neural Progenitor",
    "Radial Glia", "\\bGlioblast\\b", "Ependymal", "Choroid Plexus",
    "Bergmann", "Purkinje", "Granule Cell", "Motor Neuron", "Interneuron",
    "Schwann", "\\bGlia\\b", "\\bGlial\\b", "Neuroepithelial", "Tanycyte",
    "Pericyte"
  ), collapse = "|")

  cns_pattern <- paste(cns_regions, cns_celltypes, sep = "|")
  cns_hits    <- db_df[grepl(cns_pattern, db_df$Term, ignore.case = TRUE), ]

  if (nrow(cns_hits) > 0)
    return(cns_hits %>% arrange(Adjusted.P.value) %>% slice_head(n = 1) %>% pull(Term))

  best <- db_df %>% arrange(Adjusted.P.value) %>% slice_head(n = 1) %>% pull(Term)
  paste0("[non-CNS] ", best)
}

# ── EnrichR annotation ────────────────────────────────────────────────────────
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
        res <- enrichr(top_genes, ANNOT$enrichr_dbs)
        # Remove mouse entries
        res <- lapply(res, function(df)
          df[!grepl("Mus musculus|\\bMouse\\b|\\bmouse\\b", df$Term), ])
        comm_key <- as.character(comm)
        top_hit  <- pick_cns_label(res[[ANNOT$primary_db]])
        cell_type_labels[comm_key] <- if (!is.na(top_hit)) top_hit else comm_key
        for (db in names(res))
          write.csv(res[[db]],
                    file.path(OUT_DIR, sprintf("step4_enrichr_%s_%s.csv",
                                               comm_key,
                                               gsub("[^A-Za-z0-9]","_",db))),
                    row.names = FALSE)
      }, error = function(e) cat(sprintf("    [WARN] EnrichR %s: %s\n", comm, e$message)))
    }
    cat("\n  Labels (review step4_enrichr_*.csv before proceeding):\n")
    for (comm in names(cell_type_labels))
      cat(sprintf("    %-6s → %s\n", comm, cell_type_labels[comm]))
  }, error = function(e)
    cat(sprintf("  [WARN] EnrichR failed: %s — using community IDs.\n", e$message)))
} else {
  cat("  enrichR not installed — labels default to community IDs.\n")
}

seurat_int$cell_type_label <- unname(cell_type_labels[seurat_int$community])

top5       <- sig_markers %>% group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>% ungroup()
plot_genes <- unique(top5$gene)

pdf(file.path(OUT_DIR, "step4_annotation.pdf"),
    width = max(16, length(plot_genes) * 0.4), height = 8)
tryCatch(
  print(DotPlot(seurat_int, features = plot_genes, group.by = "community",
                dot.min = 0.05) +
          RotatedAxis() +
          scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
          labs(title = "Top 5 marker genes per community", x = NULL, y = "Community") +
          theme(axis.text.x = element_text(size = 6, face = "italic"))),
  error = function(e) message("  [WARN] DotPlot: ", e$message)
)
tryCatch(
  print(DimPlot(seurat_int, reduction = "umap", group.by = "cell_type_label",
                label = TRUE, label.size = 3, pt.size = 0.3, repel = TRUE) +
          ggtitle("UMAP — cell type labels") + theme_bw() +
          theme(legend.text = element_text(size = 7))),
  error = function(e) message("  [WARN] DimPlot: ", e$message)
)
dev.off()
cat("  Saved: step4_annotation.pdf\n")

saveRDS(seurat_int, file.path(OUT_DIR, "step4_seurat_annotated.rds"))
cat("  Saved: step4_seurat_annotated.rds\n")

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
cat("  Saved: cell_type_metadata.csv, seurat_integrated.rds\n")
cat("\n====== pipeline_cleaned.R COMPLETE — proceed to pipeline.py ======\n")
