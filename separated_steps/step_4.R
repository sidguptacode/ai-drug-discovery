#!/usr/bin/env Rscript
# =============================================================================
# step_4.R — Marker genes + EnrichR cell-type annotation
# Reads:  OUT_DIR/step3_seurat_clustered.rds
# Writes: OUT_DIR/step4_seurat_annotated.rds
#         OUT_DIR/step4_markers.csv
#         OUT_DIR/step4_enrichr_<community>_<db>.csv
#         OUT_DIR/step4_annotation.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(yaml)
})

options(repos = c(CRAN = "https://cloud.r-project.org/"))
options(future.globals.maxSize = Inf)

try_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("  Package '%s' not installed - skipping.", pkg))
    return(FALSE)
  }
  library(pkg, character.only = TRUE); TRUE
}

has_enrichr <- try_library("enrichR")

cfg     <- yaml::read_yaml("/scratch/baderlab/sgupta/ai-drug-discovery/config.yml")
OUT_DIR <- cfg$out_dir
ANNOT   <- cfg$annotation

set.seed(42)
cat(sprintf("====== step_4.R | Annotation | Seurat v%s ======\n", packageVersion("Seurat")))
cat(sprintf("  FindAllMarkers: min_pct=%g | logfc>=%g | EnrichR top %d genes\n",
            ANNOT$min_pct, ANNOT$logfc_threshold, ANNOT$n_marker_genes))

step3_path <- file.path(OUT_DIR, "step3_seurat_clustered.rds")
if (!file.exists(step3_path)) stop("Missing: ", step3_path, "\n  Run step_3.R first.")
cat("  Loading step3_seurat_clustered.rds...\n")
seurat_int <- readRDS(step3_path)
cat(sprintf("  %d spots | communities: %s\n",
            ncol(seurat_int), paste(sort(unique(seurat_int$community)), collapse=" ")))

# Switch to RNA assay; JoinLayers collapses per-sample layers so Wilcoxon can
# run across all spots. NormalizeData creates the data layer (absent when
# SCTransform was used without a subsequent NormalizeData call).
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
      cat(sprintf("    %-6s -> %s\n", comm, cell_type_labels[comm]))
  }, error = function(e)
    cat(sprintf("  [WARN] EnrichR failed: %s - using community IDs.\n", e$message)))
} else {
  cat("  enrichR not installed - labels default to community IDs.\n")
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
          ggtitle("UMAP - cell type labels") + theme_bw() +
          theme(legend.text = element_text(size = 7))),
  error = function(e) message("  [WARN] DimPlot: ", e$message)
)
dev.off()
cat("  Saved: step4_annotation.pdf\n")

saveRDS(seurat_int, file.path(OUT_DIR, "step4_seurat_annotated.rds"))
cat("  Saved: step4_seurat_annotated.rds\n")

cat(sprintf("\n====== step_4.R COMPLETE | %d communities annotated ======\n",
            length(unique(seurat_int$community))))
