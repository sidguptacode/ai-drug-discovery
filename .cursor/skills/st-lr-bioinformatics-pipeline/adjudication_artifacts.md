# Annotation adjudication artifacts (after step 4)

Step **4** writes **intermediate** tables only (markers, per-DB EnrichR CSVs, `step4_collated_candidates.json`). It sets `cell_type_label` in the Seurat object to the **community id** (placeholder). **Final** cell-type names and scores come from the **reflecting agent**, which must write the two files below **before step 5**.

**Mechanical aggregation** is already in `step4_collated_candidates.json` (top marker genes + top-k EnrichR rows per database per cluster). The agent’s job is **not** to redo that collation, but to interpret it with dataset identity and the 3-tier policy, then emit the tabular scores/labels and the JSON report.

---

## `step4_adjudication_labels.csv`

UTF-8 CSV with header row. **Required columns**:

| Column | Type | Description |
|--------|------|-------------|
| `cluster` | string | Community id as in Seurat `community` (e.g. `C1`, `C2`). Must cover **every** community in the object (step 5 errors if any are missing). |
| `original_label` | string | Placeholder from step 4 — typically the same as `cluster` (audit only). |
| `adjudication_label` | string | **Final** cell-type label for export and steps 6–10. **Must be non-empty for every row** (step 5 applies this to all spots in the cluster). |
| `override_applied` | boolean | **Audit**: `true` if `adjudication_label` differs from `original_label`; `false` if unchanged. Step 5 does **not** use this to gate application — all `adjudication_label` values are applied. |
| `tier` | string | One of: `tier1_strong_consensus`, `tier2_partial_consensus`, `tier3_conflicted`. |
| `confidence` | number | 0–1 inclusive (recommended); may be empty only if documented in JSON. |
| `label_evidence` | string | **Required.** Plain-language justification that **cites concrete evidence** from step 4: at least one **HGNC gene symbol** from markers (e.g. from `step4_collated_candidates.json` → `key_marker_genes` or `step4_markers.csv`) **and/or** a named EnrichR **term** and database. Generic text with no symbols/terms is not acceptable. |

**Optional columns** (recommended): `adjusted_pvalue` (primary evidence P-value), `needs_review` (boolean), `notes` (string).

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
| `overridden_clusters` | array of string (optional) | Cluster ids where `override_applied` is `true` in the CSV — convenience list. |

**Per-cluster object in `clusters`** (minimum):

| Key | Type | Description |
|-----|------|-------------|
| `cluster` | string | Same as CSV `cluster`. |
| `override_applied` | boolean | **Required.** Must match the CSV row for this cluster. |
| `tier` | string | Same enum as CSV `tier`. |
| `confidence` | number | 0–1. |
| `needs_review` | boolean | `true` if label should be manually verified (e.g. tier 3). |
| `chosen_label` | string | Final label (**must match** CSV `adjudication_label`). |
| `rationale` | string | **Required.** Grounded explanation (see below). Must **not** be empty. |
| `key_marker_genes` | array of string | **Recommended** (validator may require ≥1). HGNC symbols from step 4 that support the label (subset of top markers / driver genes). |
| `candidates` | array | Top-k terms per DB: `{ "db": "...", "term": "...", "adjusted_pvalue": number, "genes": "..." }`. |
| `tier_rationale` | string | **Required.** Why this tier was assigned (ties to cross-DB agreement). |

**Evidence rule for `rationale`**: The text must explicitly connect the chosen label to **marker-level and/or enrichment-level** evidence from step 4 outputs (e.g. “high `KRT8`/`EPCAM` in `key_marker_genes` + CellMarker `Epithelial…` term”). Restating the label alone is insufficient.

**Consistency**: For each cluster, `override_applied` in JSON must equal `override_applied` in `step4_adjudication_labels.csv`, and `chosen_label` must equal `adjudication_label`. CSV `label_evidence` should align with JSON `rationale` + `key_marker_genes` / `candidates`.

---

## 3-tier policy (reference)

- **Tier 1 (`tier1_strong_consensus`)**: Strong cross-database convergence on lineage or specific type; choose best-supported specific label.
- **Tier 2 (`tier2_partial_consensus`)**: Partial agreement; prefer conservative lineage/parent label to avoid overfitting rare tissue-specific terms.
- **Tier 3 (`tier3_conflicted`)**: Conflicting or weak enrichment; broad safe label + `needs_review: true`.

---

## Validation checklist

Use after creating artifacts and before **step 5**:

1. **Files**: `step4_adjudication_labels.csv` and `step4_adjudication_report.json` exist under `OUT_DIR`.
2. **CSV**: Every cluster in `step4_collated_candidates.json` (or `step4_enrichr_marker_source.csv` if collation is absent) has a matching `cluster` row with non-empty `adjudication_label` and substantive `label_evidence` (citing markers/terms), or skips are explicitly documented in the JSON.
3. **JSON–CSV**: For each cluster, `chosen_label` matches `adjudication_label`, `override_applied` matches, and `rationale` / `tier_rationale` are grounded (validator enforces minimum length plus `key_marker_genes` and/or `candidates`).
4. **Validator**: Run `scripts/validate_step4b_artifacts.py --config <merged_or_base_config>` (or `run_step4b_validate.sh`).
5. **Step 5**: Run step 5 with `PIPELINE_RUN_ID=<run_id>`; confirm `cell_type_metadata.csv` carries the adjudicated labels.
6. **Documentation**: Update `runs/<run_id>/run_info.json` with adjudication completion and any `needs_review` clusters.

---

## Example minimal rows (CSV)

```csv
cluster,original_label,adjudication_label,override_applied,tier,confidence,label_evidence
C1,C1,Epithelial (luminal-like) Human,true,tier2_partial_consensus,0.72,"Top markers include KRT8 EPCAM; CellMarker and PanglaoDB terms converge on epithelial lineage though subtype differs."
C2,C2,Endothelial Cell Artery Human,true,tier1_strong_consensus,0.91,"PECAM1 VWF high; CellMarker_2024 term Artery Endothelial matches Human_Gene_Atlas."
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
  "key_marker_genes": ["KRT8", "EPCAM", "MUC1"],
  "rationale": "Top markers (KRT8, EPCAM) match epithelial programs; CellMarker 'Epithelial Cell Undefined Human' and PanglaoDB luminal epithelial terms partially agree; subtype conflicted so conservative luminal-like label.",
  "candidates": [
    { "db": "CellMarker_2024", "term": "Epithelial Cell Undefined Human", "adjusted_pvalue": 0.000162, "genes": "MUC1;ITGB4" }
  ],
  "tier_rationale": "Partial cross-DB agreement on epithelial lineage only."
}
```
