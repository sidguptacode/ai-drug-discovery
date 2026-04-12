# Reflection guidelines

Use these checkpoints to reflect on pipeline outputs before trusting results or changing code. Apply them after the relevant steps (e.g. after 1→2, after 4, after 8→10).

**Agent obligation**: At each checkpoint, perform **every** item below. Do not proceed to the next step until all checks are satisfied or you have applied a remediation (config change + re-run) and documented it in `runs/<run_id>/run_info.json`. Before reflecting, read `config.yml` for **dataset identity** (`species`, `disease`, `dataset_name`) and use them when evaluating outputs.

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

**When**: After step 4; **adjudication CSV + JSON must exist before step 5** (step 5 hard-fails if they are missing).

**What to do**:

- Read **config.yml** (or merged config) for **dataset identity**: `species`, `disease`, and (if present) any tissue or sample context. Many diseases imply a body region or tissue (e.g. brain, lung, liver); use that when judging labels.
- Open **`step4_collated_candidates.json`** (mechanical aggregation) and **`step4_enrichr_*.csv`** / **step4_annotation.pdf** to see marker support and candidate terms per database.
- **Write (reflecting agent)**: **`step4_adjudication_labels.csv`** (include **`label_evidence`** per cluster, citing genes/terms from step 4) and **`step4_adjudication_report.json`** (per-cluster **`rationale`**, **`key_marker_genes`**, **`candidates`** as in the skill); there is no `step4_annotation_scores.csv` from R anymore.

**Label–config consistency (mandatory)**:

- Cross-check **adjudicated** labels against config: **species** and **tissue/region** plausibility for the disease or sample.
- If any label **contradicts** dataset identity, remediate (e.g. `annotation.label_disqualify_patterns`, different `enrichr_dbs`, marker/prior tweaks), re-run step 4, re-adjudicate, and document in `runs/<run_id>/run_info.json`.

### Adjudication quality (before step 5)

- **Inputs**: Merged config identity plus **`step4_collated_candidates.json`** and raw EnrichR CSVs as needed ([adjudication_artifacts.md](adjudication_artifacts.md)).
- **3-tier policy**: For each cluster, `tier` ∈ `tier1_strong_consensus` | `tier2_partial_consensus` | `tier3_conflicted`, rationale matches tier; tier 3 ⇒ **`needs_review: true`** in JSON unless explicitly waived and documented.
- **Coverage**: Every cluster must have a **non-empty** `adjudication_label` in the CSV; **`chosen_label`** in JSON must match **`adjudication_label`**; **`override_applied`** must match in CSV and JSON (`override_applied` is audit-only — step 5 applies all adjudicated labels).
- **Pass / fail**:
  - **Pass**: Validator OK (`validate_step4b_artifacts.py`), labels consistent with dataset identity, `needs_review` documented where needed.
  - **Fail**: Missing file, empty labels, JSON/CSV mismatch, or biological contradiction — **do not run step 5**.
- **Remediation**: Revise adjudication artifacts, or adjust `annotation.*` / prior context, re-run step 4, then re-adjudicate.

**Reflect**:

- Do the **final** labels make sense for the **disease**, **sample**, and expected cell types?
- If enrichment is weak or generic, tune `annotation.enrichr_dbs`, prior genes, or `label_disqualify_patterns`, then re-run step 4 and re-adjudicate.

---

## Checkpoint 3: LR pair sanity (literature)

**When**: After step 10 (or when inspecting final LR results).

**What to do**:

- Open **GROUND_TRUTH_lr_pairs_ranked.csv** and take the top 10–20 LR pairs (and their sender/receiver cell types).
- For a subset of these, do a **quick literature / web search**: are they known to be relevant to the disease, tissue, or pathway of interest?

**Reflect**:

- Are the top pairs biologically plausible for the dataset and disease?
- If many top pairs look irrelevant or noisy, consider tightening annotation (checkpoint 2), LR/CCI parameters in config (e.g. `lr_scoring.min_spots`, `cci.significance_cutoff`), or revisiting clustering/annotation before re-running from the appropriate step.
