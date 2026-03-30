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

**When**: After step 4 (and **before step 5** if adjudication artifacts are used).

**What to do**:

- Read **config.yml** (or merged config) for **dataset identity**: `species`, `disease`, and (if present) any tissue or sample context. Many diseases imply a body region or tissue (e.g. brain, lung, liver); use that when judging labels.
- Open **step4_annotation_scores.csv**: per-cluster EnrichR label and adjusted P-value; aggregate mean/median and fraction with p < 0.05.
- Open **step4_enrichr_*.csv** (and **step4_annotation.pdf**) to see which terms were considered and the chosen labels.

**Label–config consistency (mandatory)**:

- Cross-check assigned labels against what is known from config: **species** (labels should match the dataset species when the database mixes species), and **tissue/region** (labels should be plausible for the disease or sample context—e.g. brain disease → brain-relevant cell types; lung sample → lung/airway plausible).
- If any label **contradicts** dataset identity (wrong species, or tissue/body region that does not fit the disease or sample), do not treat the run as acceptable. Remediate (e.g. `label_prefer_patterns` / `label_disqualify_patterns` to align with species and context, or change `enrichr_dbs`), re-run from step 4, and document in `runs/<run_id>/run_info.json`.

### Annotation adjudication (when Step 4b is used)

If **`outputs/<run_id>/step4_adjudication_labels.csv`** and **`step4_adjudication_report.json`** exist (or you intend to create them), complete this **before running step 5**:

- **Inputs**: Confirm adjudication used merged config identity plus top-k terms from each `step4_enrichr_<community>_<db>.csv` (see [adjudication_artifacts.md](adjudication_artifacts.md)).
- **3-tier policy**: For each cluster, verify `tier` is one of `tier1_strong_consensus` | `tier2_partial_consensus` | `tier3_conflicted`, and that the written rationale matches the tier:
  - **Tier 1**: Strong cross-DB convergence — specific label justified.
  - **Tier 2**: Partial agreement — prefer conservative lineage/parent label if specific terms disagree across DBs.
  - **Tier 3**: Conflicted or weak — broad safe label; **`needs_review` must be true** in JSON unless explicitly waived and documented.
- **Overrides**: Rows with `override_applied=true` must have a non-empty `adjudication_label` and must remain consistent with species and disease context.
- **Pass / fail**:
  - **Pass**: All clusters accounted for (or skips documented), adjudicated labels consistent with dataset identity, and any `needs_review` clusters listed for follow-up. Per-cluster **`override_applied`** in `step4_adjudication_report.json` matches the CSV (so overrides are unambiguous in both artifacts).
  - **Fail**: Contradiction with config identity, missing clusters, or JSON/CSV mismatch for `chosen_label` vs `adjudication_label`, or **`override_applied` mismatch** between JSON and CSV. **Do not run step 5** until remediated (revise adjudication files, or drop overrides and use step 4 labels only by removing/renaming the adjudication CSV or setting all `override_applied` to false).
- **Remediation**: Update adjudication artifacts and re-validate; or adjust `annotation.*` in `runs/<run_id>/config_overrides.yml` and re-run step 4, then re-adjudicate. Record in `runs/<run_id>/run_info.json`.

**Reflect**:

- Do the labels make sense in the context of the **disease**, **data sample**, and **user’s target** (e.g. cell types expected in the tissue)?
- Are EnrichR scores good enough (e.g. many clusters with adjusted P < 0.05)? If labels are generic or scores are poor, consider:
  - Changing `annotation.primary_db` or `annotation.enrichr_dbs`
  - Using `label_filter_preset`, `label_prefer_patterns`, or `label_disqualify_patterns` in config to steer or filter terms
- After adjudication, confirm **step 5** will apply overrides: see [adjudication_artifacts.md](adjudication_artifacts.md) validation checklist.

---

## Checkpoint 3: LR pair sanity (literature)

**When**: After step 10 (or when inspecting final LR results).

**What to do**:

- Open **GROUND_TRUTH_lr_pairs_ranked.csv** and take the top 10–20 LR pairs (and their sender/receiver cell types).
- For a subset of these, do a **quick literature / web search**: are they known to be relevant to the disease, tissue, or pathway of interest?

**Reflect**:

- Are the top pairs biologically plausible for the dataset and disease?
- If many top pairs look irrelevant or noisy, consider tightening annotation (checkpoint 2), LR/CCI parameters in config (e.g. `lr_scoring.min_spots`, `cci.significance_cutoff`), or revisiting clustering/annotation before re-running from the appropriate step.
