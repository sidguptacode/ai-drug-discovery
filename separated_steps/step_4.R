#!/usr/bin/env Rscript
# =============================================================================
# step_4.R — Marker genes + EnrichR cell-type annotation
# Reads:  OUT_DIR/step3_seurat_clustered.rds
# Writes: OUT_DIR/step4_seurat_annotated.rds
#         OUT_DIR/step4_markers.csv
#         OUT_DIR/step4_enrichr_<community>_<db>.csv
#         OUT_DIR/step4_annotation.pdf
#         OUT_DIR/step4_annotation_scores.csv
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(yaml)
  library(jsonlite)
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

cfg_path <- Sys.getenv("PIPELINE_STEP_CONFIG", "config.yml")
if (grepl("\\.json$", cfg_path, ignore.case = TRUE)) {
  cfg <- jsonlite::read_json(cfg_path, simplifyVector = FALSE)
} else {
  cfg <- yaml::read_yaml(cfg_path)
}
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

markers_path <- file.path(OUT_DIR, "step4_markers.csv")
if (file.exists(markers_path)) {
  cat("  Resuming: step4_markers.csv found, skipping FindAllMarkers.\n")
  sig_markers <- read.csv(markers_path)
  cat(sprintf("  Loaded %d markers across %d communities.\n",
              nrow(sig_markers), length(unique(sig_markers$cluster))))
} else {
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
  write.csv(sig_markers, markers_path, row.names = FALSE)
  cat(sprintf("  Found %d significant markers across %d communities.\n",
              nrow(sig_markers), length(unique(sig_markers$cluster))))
}

# ── Label filter presets (CNS vs Blank) ───────────────────────────────────────
# CNS = disqualify mouse/species, prefer CNS regions+cell types, fallback prefix "[non-CNS]"
# Blank = no disqualify, no prefer, no prefix (raw best EnrichR term; for testing)
DEFAULT_DISQUALIFY_CNS <- c("Mus musculus", "\\bMouse\\b", "\\bmouse\\b")
DEFAULT_PREFER_CNS <- c(
  "Brain", "Cortex", "Cerebellum", "Cerebellar", "Hippocampus",
  "Spinal Cord", "Spinal", "\\bCNS\\b", "Brainstem", "Striatum",
  "Thalamus", "Hypothalamus", "Midbrain", "\\bPons\\b", "Medulla",
  "Olfactory", "Amygdala", "Prefrontal", "White Matter", "Grey Matter",
  "Gray Matter", "Frontal Lobe", "Temporal Lobe", "Ventricle",
  "Substantia Nigra", "Basal Ganglia", "Cerebrum",
  "Astrocyte", "Oligodendrocyte", "\\bOPC\\b", "Oligodendrocyte Precursor",
  "Microglia", "\\bNeuron\\b", "\\bNeuronal\\b", "Neural Progenitor",
  "Radial Glia", "\\bGlioblast\\b", "Ependymal", "Choroid Plexus",
  "Bergmann", "Purkinje", "Granule Cell", "Motor Neuron", "Interneuron",
  "Schwann", "\\bGlia\\b", "\\bGlial\\b", "Neuroepithelial", "Tanycyte",
  "Pericyte"
)
DEFAULT_FALLBACK_PREFIX_CNS <- "[non-CNS] "

# Coerce config value to character vector (handles JSON list or YAML vector)
config_char_vec <- function(x) {
  if (is.null(x)) return(NULL)
  as.character(unlist(x))
}

preset <- tolower(trimws(if (is.null(ANNOT$label_filter_preset) || !nzchar(trimws(ANNOT$label_filter_preset))) "CNS" else ANNOT$label_filter_preset))
disqualify_list <- config_char_vec(ANNOT$label_disqualify_patterns)
if (is.null(disqualify_list)) disqualify_list <- if (preset == "blank") character(0) else DEFAULT_DISQUALIFY_CNS
prefer_list <- config_char_vec(ANNOT$label_prefer_patterns)
if (is.null(prefer_list)) prefer_list <- if (preset == "blank") character(0) else DEFAULT_PREFER_CNS
fallback_prefix <- if (!is.null(ANNOT$label_fallback_prefix)) trimws(as.character(ANNOT$label_fallback_prefix)[1]) else if (preset == "blank") "" else DEFAULT_FALLBACK_PREFIX_CNS

disqualify_regex <- if (length(disqualify_list) > 0) paste(disqualify_list, collapse = "|") else NULL
cat(sprintf("  Label filter: preset=%s | disqualify=%d patterns | prefer=%d patterns\n",
            preset, length(disqualify_list), length(prefer_list)))

# ── Preferred-label picker ────────────────────────────────────────────────────
# prefer_patterns = character vector of regex; fallback_prefix = string (can be "")
# If no prefer patterns: return best row by p-value with no prefix.
# Else: prefer terms matching any pattern; if none match, use best overall and prefix.
# Returns list(term = character, adj_pvalue = numeric)
pick_preferred_label <- function(db_df, prefer_patterns, fallback_prefix) {
  if (nrow(db_df) == 0) return(list(term = NA_character_, adj_pvalue = NA_real_))
  row <- db_df %>% arrange(Adjusted.P.value) %>% slice_head(n = 1)
  best_term <- row %>% pull(Term)
  best_pval <- row %>% pull(Adjusted.P.value)
  if (length(prefer_patterns) == 0) {
    return(list(term = best_term, adj_pvalue = best_pval))
  }
  prefer_regex <- paste(prefer_patterns, collapse = "|")
  preferred <- db_df[grepl(prefer_regex, db_df$Term, ignore.case = TRUE), ]
  if (nrow(preferred) > 0) {
    row <- preferred %>% arrange(Adjusted.P.value) %>% slice_head(n = 1)
    return(list(term = row %>% pull(Term), adj_pvalue = row %>% pull(Adjusted.P.value)))
  }
  out_term <- if (nzchar(fallback_prefix)) paste0(fallback_prefix, best_term) else best_term
  return(list(term = out_term, adj_pvalue = best_pval))
}

# ── EnrichR annotation ────────────────────────────────────────────────────────
all_communities  <- sort(unique(seurat_int$community))
cell_type_labels <- setNames(all_communities, all_communities)
annotation_adj_pvalue <- setNames(rep(NA_real_, length(all_communities)), all_communities)

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
        if (!is.null(disqualify_regex))
          res <- lapply(res, function(df) df[!grepl(disqualify_regex, df$Term), ])
        comm_key <- as.character(comm)
        pick     <- pick_preferred_label(res[[ANNOT$primary_db]], prefer_list, fallback_prefix)
        cell_type_labels[comm_key] <- if (!is.na(pick$term)) pick$term else comm_key
        if (!is.na(pick$adj_pvalue)) annotation_adj_pvalue[comm_key] <- pick$adj_pvalue
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

# ── Per-cluster and aggregate annotation quality (EnrichR adjusted P-value) ───
cat("\n  Per-cluster annotation P-value (EnrichR adjusted P-value for chosen label):\n")
for (comm in names(annotation_adj_pvalue)) {
  pv <- annotation_adj_pvalue[comm]
  pv_str <- if (is.na(pv)) "NA" else sprintf("%.4g", pv)
  cat(sprintf("    cluster %-6s  adjusted_pvalue = %s\n", comm, pv_str))
}
pv_valid <- annotation_adj_pvalue[!is.na(annotation_adj_pvalue)]
if (length(pv_valid) > 0) {
  agg_mean   <- mean(pv_valid)
  agg_median <- median(pv_valid)
  agg_frac_sig <- mean(pv_valid < 0.05)
  cat(sprintf("\n  Aggregate annotation quality:\n"))
  cat(sprintf("    mean(adjusted_pvalue)     = %.4g  (lower = better)\n", agg_mean))
  cat(sprintf("    median(adjusted_pvalue)   = %.4g  (lower = better)\n", agg_median))
  cat(sprintf("    fraction with p < 0.05    = %.2f  (higher = better)\n", agg_frac_sig))
} else {
  agg_mean <- agg_median <- NA_real_
  agg_frac_sig <- NA_real_
  cat("\n  No EnrichR p-values available - aggregate scores omitted.\n")
}

annotation_scores_df <- data.frame(
  cluster = c(names(annotation_adj_pvalue), "aggregate_mean", "aggregate_median", "fraction_p_under_0.05"),
  label   = c(unname(cell_type_labels[names(annotation_adj_pvalue)]), "", "", ""),
  adjusted_pvalue = c(unname(annotation_adj_pvalue), if (length(pv_valid) > 0) c(agg_mean, agg_median, agg_frac_sig) else c(NA, NA, NA))
)
write.csv(annotation_scores_df, file.path(OUT_DIR, "step4_annotation_scores.csv"), row.names = FALSE)
cat(sprintf("\n  Saved: step4_annotation_scores.csv\n"))

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
