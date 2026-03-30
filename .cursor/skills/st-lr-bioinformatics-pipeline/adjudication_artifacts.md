# Annotation adjudication artifacts (Step 4b)

Artifacts live under **`outputs/<run_id>/`** (same as `OUT_DIR`). Canonical step 4 outputs are unchanged; these files **override** labels at step 5 when `override_applied` is true.

---

## `step4_adjudication_labels.csv`

UTF-8 CSV with header row. **Required columns**:

| Column | Type | Description |
|--------|------|-------------|
| `cluster` | string | Community id as in Seurat `community` (e.g. `C1`, `C2`). Must match `as.character(seurat_int$community)` for spots to update. |
| `original_label` | string | Label from step 4 (`cell_type_label` before override); for audit. |
| `adjudication_label` | string | Final chosen label for this cluster (may equal `original_label` if no wording change). |
| `override_applied` | boolean | `true` / `false` (or `TRUE`/`FALSE`). Only rows parsed as true cause step 5 to replace `cell_type_label` for that cluster. |
| `tier` | string | One of: `tier1_strong_consensus`, `tier2_partial_consensus`, `tier3_conflicted`. |
| `confidence` | number | 0–1 inclusive (recommended); may be empty only if documented in JSON. |

**Optional columns** (recommended for audit): `needs_review` (boolean), `notes` (string).

**Semantics**:

- If `override_applied=false`, step 5 leaves `cell_type_label` as in the RDS (step 4) for that cluster; `adjudication_label` may still document the adjudicator’s preferred wording for the report.
- If `override_applied=true`, step 5 sets all spots in `cluster` to `adjudication_label`.

---

## `step4_adjudication_report.json`

Single JSON object. **Minimum top-level keys**:

| Key | Type | Description |
|-----|------|-------------|
| `schema_version` | string | e.g. `"1.0"`. |
| `run_id` | string | Pipeline run id. |
| `dataset_name`, `species`, `disease` | string | Copied from merged config at adjudication time. |
| `created_at` | string | ISO-8601 timestamp. |
| `adjudication_method` | string | e.g. `"llm"`, `"manual"`, `"hybrid"`. |
| `model` | string or null | If LLM: model id; else `null`. |
| `prompt_id` | string or null | Optional hash or id of prompt template for reproducibility. |
| `clusters` | array | One object per cluster (see below). |
| `overridden_clusters` | array of string (optional) | Cluster ids where `override_applied` is `true` — convenience for quick filtering without scanning `clusters`. |

**Per-cluster object in `clusters`** (minimum):

| Key | Type | Description |
|-----|------|-------------|
| `cluster` | string | Same as CSV `cluster`. |
| `override_applied` | boolean | **Required.** Must match the CSV row for this cluster: `true` if step 5 will replace `cell_type_label` with `chosen_label` / `adjudication_label`; `false` if the Seurat object keeps the step-4 label. Same meaning as CSV `override_applied`. |
| `tier` | string | Same enum as CSV `tier`. |
| `confidence` | number | 0–1. |
| `needs_review` | boolean | `true` if label should be manually verified (e.g. tier 3). |
| `chosen_label` | string | Final label (should match CSV `adjudication_label`). |
| `rationale` | string | Short explanation. |
| `candidates` | array | Top-k terms per DB: `{ "db": "...", "term": "...", "adjusted_pvalue": number, "genes": "..." }`. |
| `tier_rationale` | string | Why this tier was assigned. |

Additional keys are allowed (e.g. `web_searches_performed`).

**Consistency**: For each cluster, `override_applied` in JSON must equal `override_applied` in `step4_adjudication_labels.csv`. If you omit `overridden_clusters` at the top level, derive it as the list of `cluster` values where `override_applied` is `true`.

---

## 3-tier policy (reference)

- **Tier 1 (`tier1_strong_consensus`)**: Strong cross-database convergence on lineage or specific type; choose best-supported specific label.
- **Tier 2 (`tier2_partial_consensus`)**: Partial agreement; prefer conservative lineage/parent label to avoid overfitting rare tissue-specific terms.
- **Tier 3 (`tier3_conflicted`)**: Conflicting or weak enrichment; broad safe label + `needs_review: true`.

---

## Validation checklist (rollout)

Use after creating artifacts and before **step 5** (example run id: `5ebbc253-9584-41db-9530-1337776b41c1`):

1. **Files**: `step4_adjudication_labels.csv` and `step4_adjudication_report.json` exist under `outputs/<run_id>/`.
2. **CSV**: Every community from `step4_annotation_scores.csv` (excluding aggregate rows) has a matching `cluster` row, or JSON documents skips.
3. **Overrides**: For each row with `override_applied=true`, `adjudication_label` is non-empty and `cluster` exists in the Seurat object’s `community`.
4. **JSON–CSV**: For each cluster, `chosen_label` matches `adjudication_label` in the CSV, and **`override_applied` matches** in both files (or discrepancies noted). If present, `overridden_clusters` lists exactly those clusters with `override_applied=true`.
5. **Step 5**: Run step 5 with `PIPELINE_RUN_ID=<run_id>`; open `cell_type_metadata.csv` and confirm `cell_type_label` reflects overrides for overridden clusters.
6. **Downstream**: Re-run or verify step 6+ read the updated metadata; spot-check one h5ad `obs` for `cell_type_label` if needed.
7. **Documentation**: Update `runs/<run_id>/run_info.json` with adjudication completion and any `needs_review` clusters.

---

## Example minimal rows (CSV)

```csv
cluster,original_label,adjudication_label,override_applied,tier,confidence
C1,Epithelial Cell Undefined Human,Epithelial (luminal-like) Human,true,tier2_partial_consensus,0.72
C2,Endothelial Cell Artery Human,Endothelial Cell Artery Human,false,tier1_strong_consensus,0.91
```

## Example minimal `clusters[]` entry (JSON)

```json
{
  "cluster": "C1",
  "override_applied": true,
  "tier": "tier2_partial_consensus",
  "confidence": 0.72,
  "needs_review": false,
  "chosen_label": "Epithelial (luminal-like) Human",
  "rationale": "PanglaoDB and CellMarker both support epithelial programs; specific subtype terms disagree; conservative label.",
  "candidates": [
    { "db": "CellMarker_2024", "term": "Epithelial Cell Undefined Human", "adjusted_pvalue": 0.000162, "genes": "MUC1;ITGB4" }
  ],
  "tier_rationale": "Partial cross-DB agreement on epithelial lineage only."
}
```
