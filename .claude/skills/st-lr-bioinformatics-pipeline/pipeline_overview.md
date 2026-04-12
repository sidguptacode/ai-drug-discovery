# Pipeline overview

## Goal

Transform **raw spatial transcriptomics data** (Visium-style: per-sample H5 + tissue positions + scalefactors) into **extracted, targetable ligand–receptor pairs** with sender and receiver cell-type labels.

## Data flow

1. **.data** — Per-sample folders with `*_filtered_feature_bc_matrix.h5`, `*_tissue_positions_list.csv`, `*_scalefactors_json.json`.
2. **R pipeline (steps 1–5)** — QC → integration → clustering → **step 4** (intermediate markers + EnrichR + `step4_collated_candidates.json`) → **reflecting agent** writes **`step4_adjudication_labels.csv`** and **`step4_adjudication_report.json`** → **step 5** exports `cell_type_metadata.csv` and integrated Seurat object. Step 5 **requires** the two adjudication files and applies **`adjudication_label`** to every cluster.
3. **Python pipeline (steps 6–10)** — Load samples with cell-type metadata → preprocess → LR scoring (stlearn) → CCI → aggregate and rank LR pairs.

Final outputs of interest: `GROUND_TRUTH_lr_pairs_ranked.csv` and per-sample CCI/LR outputs in `out_dir`.

## Key steps

### Step 1–2: QC and integration

- **Step 1** does spot-level QC (nCount, optional %MT), loads spatial coordinates, and writes per-sample Seurat objects and `step1_qc.pdf`.
- **Step 2** runs SCTransform per sample, RPCA-based integration, then PCA/UMAP/neighbours on the integrated assay; writes `step2_integration.pdf`.

**Check**: Inspect `step1_qc.pdf` and `step2_integration.pdf` for even spatial and count distributions and for batch mixing (UMAP by sample). Adjust `config.yml` QC and integration parameters if needed.

### Step 3: Clustering

- **Step 3** runs FindClusters at the resolutions listed in config and assigns a single “chosen” resolution to define communities (e.g. C1, C2, …).

**Check**: Use `step3_clustree.pdf` and `step3_clusters.pdf` to confirm the number of communities and sample proportions. Set `clustering.chosen_resolution` in config accordingly.

### Step 4: Annotation

- **Step 4** runs FindAllMarkers, queries EnrichR, writes per-DB CSVs and **`step4_collated_candidates.json`** (mechanical top-k collation). It does **not** assign final cell-type names; the Seurat object keeps `cell_type_label` equal to the community id until step 5.

**Check**: Review enrichment and collation for plausibility. The **reflecting agent** must then write adjudication CSV + JSON before step 5. Tune `annotation.*` (e.g. `enrichr_dbs`, `label_disqualify_patterns`, disease prior) and re-run step 4 if evidence is poor.

### Step 4b: Collation + adjudication

- **Collation** is produced automatically as **`step4_collated_candidates.json`**.
- **Adjudication** is **mandatory** before step 5: **`step4_adjudication_labels.csv`** and **`step4_adjudication_report.json`** ([adjudication_artifacts.md](adjudication_artifacts.md)), **3-tier** policy.
- **Step 5** fails without both files and sets **`cell_type_label`** from **`adjudication_label`** for all clusters.

### Step 8: LR scoring

- **Step 8** uses stlearn (e.g. CellPhoneDB) to compute per-spot LR scores and p-values. These feed into step 9 (CCI) and step 10 (ranked LR pairs).

**Check**: After the full pipeline, use reflection guidelines to sanity-check top LR pairs (e.g. literature/disease relevance).

## Config

A single **config.yml** at repo root (or path in `PIPELINE_STEP_CONFIG`) drives the pipeline:

- **Dataset**: `dataset_name`, `disease`, `species`
- **Paths**: `data_dir`, `out_dir`
- **Samples**: list of sample IDs (folder names under `data_dir`)
- **Per-step blocks**: `qc`, `integration`, `clustering`, `annotation`, `preprocessing`, `lr_scoring`, `cci`

Change these to run the same pipeline on another dataset or to tune QC, clustering, annotation, or LR/CCI behaviour.
