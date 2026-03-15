# Step specifications

All steps read config from `config.yml` at repo root, or from the path in `PIPELINE_STEP_CONFIG`. Paths below are relative to `OUT_DIR` (and `DATA_DIR` for raw data) defined in config.

## Waiting for RDS writes

R steps (1–5) write large `.rds` files. The process may return or log "complete" before the file is fully written to disk (buffering/sync). **Before starting the step that consumes an RDS** (e.g. step 5 reads `step4_seurat_annotated.rds`), ensure the file has finished growing: monitor its size (e.g. `stat -c%s OUT_DIR/stepN_*.rds`) every 2–3 seconds until the value is unchanged for at least 2–3 consecutive checks. Do not run the next step until the RDS size is stable; otherwise the reader may see a truncated or corrupt file.

## .data layout

`data_dir` in config (e.g. `.data`) must contain one folder per sample. Per sample (e.g. `GSM8608424_VLP73_A1`):

- `<sample>_filtered_feature_bc_matrix.h5` — required for step 1 and step 6
- `<sample>_tissue_positions_list.csv` — required; step 6 expects columns that can provide full-res pixel coords (e.g. `pxl_row_in_fullres` / `pxl_col_in_fullres`) for spatial/stlearn
- `<sample>_scalefactors_json.json` — required for step 6

---

## Step 1 — QC

- **Purpose**: Load H5 count matrices and tissue positions per sample, filter spots by nCount (and optionally % mitochondrial).
- **Libraries**: Seurat, ggplot2, patchwork, dplyr, yaml (R).
- **Inputs**: For each sample in config: `DATA_DIR/<sample>/<sample>_filtered_feature_bc_matrix.h5`, `DATA_DIR/<sample>/<sample>_tissue_positions_list.csv`.
- **Required outputs**: `step1_seurat_list.rds`, `step1_qc.pdf`.

---

## Step 2 — Integration

- **Purpose**: SCTransform per sample, RPCA anchor integration, IntegrateData, then PCA, UMAP, and neighbour graph on integrated assay.
- **Libraries**: Seurat, ggplot2, yaml (R); optional harmony.
- **Inputs**: `step1_seurat_list.rds`.
- **Required outputs**: `step2_seurat_integrated.rds`, `step2_integration.pdf`.

---

## Step 3 — Clustering

- **Purpose**: FindClusters at config resolutions, assign community IDs; optional clustree for resolution choice.
- **Libraries**: Seurat, ggplot2, dplyr, yaml (R); optional clustree.
- **Inputs**: `step2_seurat_integrated.rds`.
- **Required outputs**: `step3_seurat_clustered.rds`, `step3_clusters.pdf`; `step3_clustree.pdf` if clustree is installed and more than one resolution is used.

---

## Step 4 — Annotation

- **Purpose**: FindAllMarkers per community, EnrichR on top marker genes, label selection/filtering per cluster.
- **Libraries**: Seurat, ggplot2, dplyr, yaml (R); optional enrichR.
- **Inputs**: `step3_seurat_clustered.rds`.
- **Required outputs**: `step4_seurat_annotated.rds`, `step4_markers.csv`, `step4_enrichr_<community>_<db>.csv` (per community and EnrichR DB), `step4_annotation.pdf`, `step4_annotation_scores.csv`.

---

## Step 5 — Export

- **Purpose**: Export cell-type metadata and integrated Seurat object for the Python pipeline.
- **Libraries**: Seurat, yaml (R).
- **Inputs**: `step4_seurat_annotated.rds`.
- **Required outputs**: `cell_type_metadata.csv`, `seurat_integrated.rds`.

---

## Step 6 — Load samples (Python)

- **Purpose**: Load per-sample AnnData from raw H5 and spatial files, attach community and cell_type_label from step 5, set up spatial coordinates for stlearn.
- **Libraries**: scanpy, pandas, yaml (Python).
- **Inputs**: `cell_type_metadata.csv`; per sample: `DATA_DIR/<sample>/<sample>_filtered_feature_bc_matrix.h5`, `*_tissue_positions_list.csv`, `*_scalefactors_json.json`.
- **Required outputs**: `step6_<sample>.h5ad` (one per sample).

---

## Step 7 — Preprocess (Python)

- **Purpose**: Filter cells/genes, normalize total, log1p transform.
- **Libraries**: scanpy, matplotlib, yaml (Python).
- **Inputs**: `step6_<sample>.h5ad` (each sample).
- **Required outputs**: `step7_<sample>.h5ad` (each sample), `step7_count_distributions.pdf`.

---

## Step 8 — LR scoring (Python)

- **Purpose**: LR scoring with stlearn (e.g. CellPhoneDB); per-spot LR scores and p-values.
- **Libraries**: stlearn, scanpy, yaml (Python).
- **Inputs**: `step7_<sample>.h5ad` (each sample).
- **Required outputs**: `step8_<sample>.h5ad`, `<sample>_lr_scores.csv`, `<sample>_lr_pvals.csv`, `step8_<sample>_lr_spatial.pdf` (each sample).

---

## Step 9 — CCI (Python)

- **Purpose**: Cell–cell interaction analysis with stlearn (run_cci using cell_type_label); retain significant interactions.
- **Libraries**: stlearn, scanpy, yaml (Python).
- **Inputs**: `step8_<sample>.h5ad` (each sample).
- **Required outputs**: `<sample>_cci_results.csv` (one per sample).

---

## Step 10 — Aggregate and rank (Python)

- **Purpose**: Aggregate CCI results across samples and rank LR pairs (e.g. by recurrence and mean score).
- **Libraries**: pandas, yaml, matplotlib (Python).
- **Inputs**: `<sample>_cci_results.csv` (all samples).
- **Required outputs**: `GROUND_TRUTH_lr_pairs_ranked.csv`, `step10_ranked_lr_pairs.pdf`.
