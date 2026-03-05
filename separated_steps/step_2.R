#!/usr/bin/env Rscript
# =============================================================================
# step_2.R — Integration (SCTransform per sample + RPCA anchor integration)
# Reads:  OUT_DIR/step1_seurat_list.rds
# Writes: OUT_DIR/step2_seurat_integrated.rds
#         OUT_DIR/step2_integration.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(yaml)
  library(jsonlite)
})

options(repos = c(CRAN = "https://cloud.r-project.org/"))
options(future.globals.maxSize = Inf)

if (!requireNamespace("harmony", quietly = TRUE))
  install.packages("harmony")

ts <- function(msg) {
  cat(sprintf("  [%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  invisible(Sys.time())
}

cfg_path <- Sys.getenv("PIPELINE_STEP_CONFIG", "config.yml")
if (grepl("\\.json$", cfg_path, ignore.case = TRUE)) {
  cfg <- jsonlite::read_json(cfg_path, simplifyVector = FALSE)
} else {
  cfg <- yaml::read_yaml(cfg_path)
}
OUT_DIR <- cfg$out_dir
SAMPLES <- if (is.list(cfg$samples)) unlist(cfg$samples) else cfg$samples
INTEG   <- cfg$integration

set.seed(42)
cat(sprintf("====== step_2.R | Integration | Seurat v%s ======\n", packageVersion("Seurat")))
cat(sprintf("  nfeatures=%d | n_pcs=%d | n_dims=%d\n",
            INTEG$n_features, INTEG$n_pcs, INTEG$n_dims))

step1_path <- file.path(OUT_DIR, "step1_seurat_list.rds")
if (!file.exists(step1_path)) stop("Missing: ", step1_path, "\n  Run step_1.R first.")
cat("  Loading step1_seurat_list.rds...\n")
seurat_list <- readRDS(step1_path)
cat(sprintf("  Loaded %d samples.\n", length(seurat_list)))

# ── Parallel workers for SCTransform ─────────────────────────────────────────
n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_workers) || n_workers < 1L)
  n_workers <- max(1L, parallel::detectCores(logical = FALSE) %/% 2L)
cat(sprintf("  Parallel workers: %d\n", n_workers))
if (.Platform$OS.type == "unix") {
  future::plan(future::multicore, workers = n_workers)
} else {
  future::plan(future::multisession, workers = n_workers)
}

# ── SCTransform per sample ────────────────────────────────────────────────────
sct_flavor <- if (requireNamespace("glmGamPoi", quietly = TRUE)) "v2" else "v1"
cat(sprintf("  SCTransform flavor: %s\n", sct_flavor))

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
  object.list = seurat_list, nfeatures = INTEG$n_features)

cat("  Running PrepSCTIntegration...\n")
seurat_list <- PrepSCTIntegration(
  object.list = seurat_list, anchor.features = integration_features)

cat(sprintf("  Running PCA on each sample (npcs=%d)...\n", INTEG$n_dims))
seurat_list <- lapply(seurat_list, function(so)
  RunPCA(so, npcs = INTEG$n_dims, verbose = FALSE))
names(seurat_list) <- SAMPLES

# ── FindIntegrationAnchors + IntegrateData ────────────────────────────────────
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

# ── PCA + UMAP + Neighbors on integrated assay ────────────────────────────────
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
          ggtitle("UMAP - by sample (batch check)") + theme_bw()),
  error = function(e) message("  [WARN] DimPlot: ", e$message)
)
dev.off()
cat("  Saved: step2_integration.pdf\n")

saveRDS(seurat_int, file.path(OUT_DIR, "step2_seurat_integrated.rds"))
cat("  Saved: step2_seurat_integrated.rds\n")

cat(sprintf("\n====== step_2.R COMPLETE | %d spots | %d samples ======\n",
            ncol(seurat_int), length(unique(seurat_int$sample))))
