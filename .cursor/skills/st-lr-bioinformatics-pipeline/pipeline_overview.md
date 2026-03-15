# Pipeline overview

## Goal

Transform **raw spatial transcriptomics data** (Visium-style: per-sample H5 + tissue positions + scalefactors) into **extracted, targetable ligand–receptor pairs** with sender and receiver cell-type labels.

## Data flow

1. **.data** — Per-sample folders with `*_filtered_feature_bc_matrix.h5`, `*_tissue_positions_list.csv`, `*_scalefactors_json.json`.
2. **R pipeline (steps 1–5)** — QC → integration → clustering → annotation → export of `cell_type_metadata.csv` and integrated Seurat object.
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

- **Step 4** runs FindAllMarkers, queries EnrichR with top marker genes per community, and assigns cell-type labels. Outputs include marker tables, EnrichR CSVs, and `step4_annotation_scores.csv`.

**Check**: Review labels and EnrichR scores for biological plausibility. **Cross-check labels with dataset identity** in config: species and (where applicable) tissue/region implied by disease or sample; labels should not contradict these. If they do, use `label_prefer_patterns` or `label_disqualify_patterns` and re-run step 4. Tune `annotation.*` (e.g. `primary_db`, `label_filter_preset`) if labels are poor.

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
