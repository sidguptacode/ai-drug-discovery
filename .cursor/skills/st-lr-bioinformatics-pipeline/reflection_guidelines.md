# Reflection guidelines

Use these checkpoints to reflect on pipeline outputs before trusting results or changing code. Apply them after the relevant steps (e.g. after 1→2, after 4, after 8→10).

---

## Checkpoint 1: QC effectiveness

**When**: After steps 1 and 2.

**What to do**:

- Open **step1_qc.pdf**. For each sample, check:
  - nCount_RNA, nFeature_RNA, and percent.mt distributions (violin plots).
  - Spatial plots (nCount and % MT): is tissue coverage even, or did filtering remove large regions?
- Optionally open **step2_integration.pdf**: UMAP by sample — do batches mix reasonably, or are there strong sample-driven clusters?

**Reflect**:

- Is the data still sensibly distributed, or do spots/samples look over-filtered or biased?
- If QC is too strict (e.g. too many spots lost), consider relaxing `qc.ncount_min` / `qc.ncount_max` or adjusting `qc.mt_cutoff` in config. If distributions look odd, consider tightening or revisiting filters.

---

## Checkpoint 2: Annotation labels and scores

**When**: After step 4.

**What to do**:

- Open **step4_annotation_scores.csv**: per-cluster EnrichR label and adjusted P-value; aggregate mean/median and fraction with p < 0.05.
- Open **step4_enrichr_*.csv** (and **step4_annotation.pdf**) to see which terms were considered and the chosen labels.

**Reflect**:

- Do the labels make sense in the context of the **disease**, **data sample**, and **user’s target** (e.g. cell types expected in the tissue)?
- Are EnrichR scores good enough (e.g. many clusters with adjusted P < 0.05)? If labels are generic or scores are poor, consider:
  - Changing `annotation.primary_db` or `annotation.enrichr_dbs`
  - Using `label_filter_preset`, `label_prefer_patterns`, or `label_disqualify_patterns` in config to steer or filter terms

---

## Checkpoint 3: LR pair sanity (literature)

**When**: After step 10 (or when inspecting final LR results).

**What to do**:

- Open **GROUND_TRUTH_lr_pairs_ranked.csv** and take the top 10–20 LR pairs (and their sender/receiver cell types).
- For a subset of these, do a **quick literature / web search**: are they known to be relevant to the disease, tissue, or pathway of interest?

**Reflect**:

- Are the top pairs biologically plausible for the dataset and disease?
- If many top pairs look irrelevant or noisy, consider tightening annotation (checkpoint 2), LR/CCI parameters in config (e.g. `lr_scoring.min_spots`, `cci.significance_cutoff`), or revisiting clustering/annotation before re-running from the appropriate step.
