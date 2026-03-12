#!/usr/bin/env Rscript
# =============================================================================
# step_3.R — Clustering
# Reads:  OUT_DIR/step2_seurat_integrated.rds
# Writes: OUT_DIR/step3_seurat_clustered.rds
#         OUT_DIR/step3_clustree.pdf   (if clustree installed and >1 resolution)
#         OUT_DIR/step3_clusters.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(yaml)
  library(jsonlite)
})

options(repos = c(CRAN = "https://cloud.r-project.org/"))

try_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("  Package '%s' not installed - skipping.", pkg))
    return(FALSE)
  }
  library(pkg, character.only = TRUE); TRUE
}

has_clustree <- try_library("clustree")

cfg_path <- Sys.getenv("PIPELINE_STEP_CONFIG", "config.yml")
if (grepl("\\.json$", cfg_path, ignore.case = TRUE)) {
  cfg <- jsonlite::read_json(cfg_path, simplifyVector = FALSE)
} else {
  cfg <- yaml::read_yaml(cfg_path)
}
OUT_DIR <- cfg$out_dir
CLUST   <- cfg$clustering

set.seed(42)
cat(sprintf("====== step_3.R | Clustering | Seurat v%s ======\n", packageVersion("Seurat")))
cat(sprintf("  Resolutions tested: %s | Chosen: %g\n",
            paste(CLUST$resolutions, collapse=", "), CLUST$chosen_resolution))

step2_path <- file.path(OUT_DIR, "step2_seurat_integrated.rds")
if (!file.exists(step2_path)) stop("Missing: ", step2_path, "\n  Run step_2.R first.")
cat("  Loading step2_seurat_integrated.rds...\n")
seurat_int <- readRDS(step2_path)
cat(sprintf("  %d spots | %d samples\n",
            ncol(seurat_int), length(unique(seurat_int$sample))))

# Detect the SNN graph — "integrated_snn" (SCTransform) or "RNA_snn" (Harmony)
snn_graph <- grep("_snn$", names(seurat_int@graphs), value = TRUE)[1]
if (is.na(snn_graph)) stop("No SNN graph found. Check that FindNeighbors ran in step_2.R.")
cat(sprintf("  Using graph: %s\n", snn_graph))

for (res in CLUST$resolutions)
  seurat_int <- FindClusters(seurat_int, graph.name = snn_graph,
                             resolution = res, verbose = FALSE)

clustree_prefix <- paste0(snn_graph, "_res.")

if (has_clustree && length(CLUST$resolutions) >= 2) {
  pdf(file.path(OUT_DIR, "step3_clustree.pdf"), width = 12, height = 10)
  print(clustree(seurat_int, prefix = clustree_prefix) +
          ggtitle("Clustree - inspect to confirm chosen_resolution in config.yml"))
  dev.off()
  cat("  Saved: step3_clustree.pdf\n")
  cat("  >>> If chosen_resolution looks wrong, update config.yml and rerun.\n")
}

res_col <- paste0(clustree_prefix, CLUST$chosen_resolution)
Idents(seurat_int) <- res_col
seurat_int$community <- paste0("C", as.integer(seurat_int[[res_col]][[1]]))

n_comm <- length(unique(seurat_int$community))
cat(sprintf("  Resolution %g -> %d communities\n", CLUST$chosen_resolution, n_comm))
print(sort(table(seurat_int$community), decreasing = TRUE))

prop_df <- as.data.frame(
  table(Community = seurat_int$community, Sample = seurat_int$sample)
) %>% group_by(Sample) %>% mutate(Proportion = Freq / sum(Freq)) %>% ungroup()

pdf(file.path(OUT_DIR, "step3_clusters.pdf"), width = 14, height = 7)
print(DimPlot(seurat_int, group.by = "community", label = TRUE,
              pt.size = 0.3, repel = TRUE) +
        ggtitle(sprintf("UMAP - communities (res = %g)", CLUST$chosen_resolution)) +
        theme_bw())
print(ggplot(prop_df, aes(x = Sample, y = Proportion, fill = Community)) +
        geom_bar(stat = "identity") + theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Community proportions per sample"))
dev.off()
cat("  Saved: step3_clusters.pdf\n")

saveRDS(seurat_int, file.path(OUT_DIR, "step3_seurat_clustered.rds"))
cat("  Saved: step3_seurat_clustered.rds\n")

cat(sprintf("\n====== step_3.R COMPLETE | %d communities ======\n", n_comm))
