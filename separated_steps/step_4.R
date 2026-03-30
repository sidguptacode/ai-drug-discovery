#!/usr/bin/env Rscript
# =============================================================================
# step_4.R — Marker genes + EnrichR cell-type annotation
# Reads:  OUT_DIR/step3_seurat_clustered.rds
# Writes: OUT_DIR/step4_seurat_annotated.rds
#         OUT_DIR/step4_markers.csv          (full FindAllMarkers; never prior-filtered)
#         OUT_DIR/step4_filtered_markers.csv (prior-filtered rows only; audit; EnrichR may use
#           unfiltered markers per cluster when the prior matches no genes in that cluster)
#         OUT_DIR/step4_enrichr_marker_source.csv — per cluster: prior_filtered vs unfiltered_fallback
#         OUT_DIR/step4_disease_marker_prior.json — built before this step by
#           scripts/build_step4_disease_prior.py (Open Targets associations); optional hand-curate.
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
markers_filtered_path <- file.path(OUT_DIR, "step4_filtered_markers.csv")

if (file.exists(markers_path)) {
  cat("  Resuming: step4_markers.csv found, skipping FindAllMarkers.\n")
  markers_unfiltered <- read.csv(markers_path, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded %d markers across %d communities.\n",
              nrow(markers_unfiltered), length(unique(markers_unfiltered$cluster))))
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

  markers_unfiltered <- all_markers %>%
    filter(p_val_adj < 0.05) %>%
    arrange(cluster, desc(avg_log2FC))
  write.csv(markers_unfiltered, markers_path, row.names = FALSE)
  cat(sprintf("  Found %d significant markers across %d communities.\n",
              nrow(markers_unfiltered), length(unique(markers_unfiltered$cluster))))
}

# ── Prior-filtered markers (EnrichR, collation, DotPlot, Seurat labels) ───────
load_step4_prior_genes <- function(path) {
  if (!file.exists(path)) return(character(0))
  doc <- tryCatch(
    jsonlite::read_json(path, simplifyVector = TRUE),
    error = function(e) {
      warning("Could not read prior JSON: ", e$message); NULL
    }
  )
  if (is.null(doc)) return(character(0))
  g <- doc$genes
  if (is.null(g)) return(character(0))
  g <- unique(toupper(trimws(as.character(g))))
  g[nzchar(g)]
}

prior_path <- file.path(OUT_DIR, "step4_disease_marker_prior.json")
if (!is.null(ANNOT$disease_marker_prior_json)) {
  p <- trimws(as.character(ANNOT$disease_marker_prior_json)[1])
  if (nzchar(p))
    prior_path <- if (startsWith(p, "/")) p else file.path(OUT_DIR, p)
}

prior_genes <- load_step4_prior_genes(prior_path)
if (length(prior_genes) > 0) {
  gene_u <- toupper(as.character(markers_unfiltered$gene))
  sig_markers <- markers_unfiltered[gene_u %in% prior_genes, , drop = FALSE]
  cat(sprintf("  Disease prior (%s): %d gene symbols -> %d / %d marker rows (prior-filtered).\n",
              basename(prior_path), length(prior_genes),
              nrow(sig_markers), nrow(markers_unfiltered)))
  if (nrow(sig_markers) == 0)
    cat("  Disease prior: no marker rows matched the prior; all clusters will use unfiltered markers for EnrichR (per-cluster fallback).\n")
} else {
  sig_markers <- markers_unfiltered
  cat("  Disease prior: no usable gene list (missing, empty, or invalid JSON); ",
      "using all marker rows for EnrichR.\n", sep = "")
}
write.csv(sig_markers, markers_filtered_path, row.names = FALSE)
cat(sprintf("  Saved: %s\n", basename(markers_filtered_path)))

# Rows to use for EnrichR per cluster: prior-filtered if any rows exist for that cluster, else unfiltered
markers_for_cluster <- function(comm, sig_tbl, unf_tbl) {
  comm <- as.character(comm)
  s <- sig_tbl %>% filter(cluster == comm)
  if (nrow(s) > 0) return(s)
  unf_tbl %>% filter(cluster == comm)
}

cluster_has_prior_rows <- function(comm, sig_tbl) {
  comm <- as.character(comm)
  nrow(sig_tbl %>% filter(cluster == comm)) > 0
}

enrichr_min_genes <- if (!is.null(ANNOT$enrichr_min_genes)) as.integer(ANNOT$enrichr_min_genes)[1] else 5L

# Coerce config value to character vector (handles JSON list or YAML vector)
config_char_vec <- function(x) {
  if (is.null(x)) return(NULL)
  as.character(unlist(x))
}

preset <- tolower(trimws(if (is.null(ANNOT$label_filter_preset) || !nzchar(trimws(ANNOT$label_filter_preset))) "blank" else ANNOT$label_filter_preset))
disqualify_list <- config_char_vec(ANNOT$label_disqualify_patterns)
if (is.null(disqualify_list)) disqualify_list <- character(0)
prefer_list <- config_char_vec(ANNOT$label_prefer_patterns)
if (is.null(prefer_list)) prefer_list <- character(0)
fallback_prefix <- if (!is.null(ANNOT$label_fallback_prefix)) trimws(as.character(ANNOT$label_fallback_prefix)[1]) else ""

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
  enrichr_source_rows <- list()
  tryCatch({
    setEnrichrSite("Enrichr")
    comm_list <- sort(unique(as.character(seurat_int$community)))
    for (comm in comm_list) {
      m_cl <- markers_for_cluster(comm, sig_markers, markers_unfiltered)
      if (length(prior_genes) == 0L) {
        src <- "no_prior_all_markers"
      } else if (cluster_has_prior_rows(comm, sig_markers)) {
        src <- "prior_filtered"
      } else {
        src <- "unfiltered_fallback"
      }
      top_genes <- m_cl %>%
        arrange(desc(avg_log2FC)) %>%
        slice_head(n = ANNOT$n_marker_genes) %>%
        pull(gene)
      enrichr_source_rows[[length(enrichr_source_rows) + 1L]] <- data.frame(
        cluster = as.character(comm),
        enrichr_marker_source = src,
        n_genes_for_enrichr = length(top_genes),
        stringsAsFactors = FALSE
      )
      if (length(top_genes) < enrichr_min_genes) next
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
    if (length(enrichr_source_rows) > 0) {
      write.csv(bind_rows(enrichr_source_rows),
                file.path(OUT_DIR, "step4_enrichr_marker_source.csv"),
                row.names = FALSE)
      cat(sprintf("  Saved: step4_enrichr_marker_source.csv\n"))
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

# ── C1 collation JSON (markers + top EnrichR labels) ─────────────────────────
parse_overlap_num <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(as.character(x))) return(NA_real_)
  parts <- strsplit(as.character(x), "/", fixed = TRUE)[[1]]
  as.numeric(parts[1])
}

build_cluster_collation <- function(cluster_id = "C1",
                                    marker_top_n = 15,
                                    label_top_n = 2) {
  cluster_markers <- markers_for_cluster(cluster_id, sig_markers, markers_unfiltered) %>%
    arrange(desc(avg_log2FC), p_val_adj)
  if (nrow(cluster_markers) == 0) {
    return(list(
      cluster = cluster_id,
      error = sprintf("No markers found for cluster %s", cluster_id)
    ))
  }

  key_markers <- cluster_markers %>%
    slice_head(n = marker_top_n) %>%
    transmute(
      gene = as.character(gene),
      avg_log2FC = as.numeric(avg_log2FC),
      p_val_adj = as.numeric(p_val_adj)
    )

  top_marker_genes <- unique(as.character((cluster_markers %>% slice_head(n = 30))$gene))
  db_list <- as.character(unlist(ANNOT$enrichr_dbs))
  top_labels <- list()
  evidence_genes_all <- character(0)

  for (db in db_list) {
    db_stub <- gsub("[^A-Za-z0-9]", "_", db)
    db_path <- file.path(OUT_DIR, sprintf("step4_enrichr_%s_%s.csv", cluster_id, db_stub))
    if (!file.exists(db_path)) {
      top_labels[[db]] <- list()
      next
    }
    db_df <- read.csv(db_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(db_df) == 0 || !all(c("Term", "Adjusted.P.value", "Combined.Score", "Genes", "Overlap") %in% colnames(db_df))) {
      top_labels[[db]] <- list()
      next
    }
    db_df <- db_df %>% arrange(Adjusted.P.value, desc(Combined.Score))
    top_df <- db_df %>% slice_head(n = label_top_n)
    labels_this_db <- lapply(seq_len(nrow(top_df)), function(i) {
      row <- top_df[i, , drop = FALSE]
      genes <- strsplit(as.character(row$Genes), ";", fixed = TRUE)[[1]]
      genes <- trimws(genes)
      genes <- genes[nzchar(genes)]
      evidence_genes_all <<- c(evidence_genes_all, genes)
      list(
        label = as.character(row$Term),
        adjusted_p_value = as.numeric(row$Adjusted.P.value),
        combined_score = as.numeric(row$Combined.Score),
        overlap = as.character(row$Overlap),
        overlap_count = parse_overlap_num(row$Overlap),
        marker_evidence_genes = unname(genes)
      )
    })
    top_labels[[db]] <- labels_this_db
  }

  recurrent_df <- data.frame(
    gene = names(sort(table(evidence_genes_all), decreasing = TRUE)),
    count_across_top_labels = as.integer(sort(table(evidence_genes_all), decreasing = TRUE)),
    stringsAsFactors = FALSE
  )
  recurrent_df <- recurrent_df %>% filter(count_across_top_labels >= 2)
  recurrent_and_top_markers <- recurrent_df %>% filter(gene %in% top_marker_genes)

  list(
    cluster = cluster_id,
    key_marker_genes = key_markers,
    top_candidate_labels_by_database = top_labels,
    markers_mainly_driving_label_assignment = list(
      recurrent_across_top_labels = recurrent_df,
      recurrent_and_also_top_C1_markers = recurrent_and_top_markers
    )
  )
}

c1_collation <- build_cluster_collation(cluster_id = "C1", marker_top_n = 15, label_top_n = 2)
c1_collation_path <- file.path(OUT_DIR, "step4_c1_collation.json")
jsonlite::write_json(c1_collation, c1_collation_path, pretty = TRUE, auto_unbox = TRUE, na = "null")
cat(sprintf("  Saved: %s\n", basename(c1_collation_path)))

seurat_int$cell_type_label <- unname(cell_type_labels[seurat_int$community])

top5 <- bind_rows(lapply(sort(unique(as.character(seurat_int$community))), function(cx) {
  markers_for_cluster(cx, sig_markers, markers_unfiltered) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 5) %>%
    ungroup()
}))
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
