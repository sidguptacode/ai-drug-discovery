# ST-LR Bioinformatics Pipeline Skill

Run and iteratively improve the spatial transcriptomics ligand-receptor pipeline (`separated_steps`). Use when running the ST-LR pipeline, debugging steps, or interpreting QC, annotation, or LR pair outputs.

## Initialization

Before running or changing the pipeline, read:

1. [pipeline_overview.md](pipeline_overview.md)
2. [step_specs.md](step_specs.md)
3. [reflection_guidelines.md](reflection_guidelines.md)
4. [adjudication_artifacts.md](adjudication_artifacts.md) — when performing or reviewing Step 4b (annotation adjudication).

## Disease Prior Obligation (Step 4)

There is **no** automated Open Targets / GraphQL step. The **controller agent** must create **`OUT_DIR/step4_disease_marker_prior.json`** (under `outputs/<run_id>/` when using `PIPELINE_RUN_ID`) **before** `run_step_4.sh` whenever top-level **`disease`** is set in merged config. Step 4 **stops** if that file is missing or if **`genes`** is absent or empty.

When curating **`genes`** (HGNC symbols for the species in config), reason explicitly in three buckets — **most weight on (3)**:

1. **Disease-linked** — markers and pathways directly associated with the configured disease.
2. **Locus / anatomical context** — genes reflecting expected tissue or regional identity aligned with the dataset.
3. **Canonical cell types** — established markers for cell types expected in this spatial transcriptomics analysis (broad coverage across expected lineages); this bucket should carry the **largest** share of the list.

Use your knowledge and **web search** as needed to verify symbols and marker lists. Shape: [step4_disease_marker_prior.example.json](step4_disease_marker_prior.example.json). Override path: **`annotation.disease_marker_prior_json`** (relative to `out_dir` or absolute).

## Static vs Run-Specific Files

- **Static (read-only for run-specific changes):**
  - `config.yml` — Base dataset identity, paths, and default parameters. Keep `run_id: null`.
  - `separated_steps/` — Canonical pipeline code. Do not edit for run-specific behaviour.

- **Run-specific (editable):**
  - `runs/<run_id>/config_overrides.yml` — Holds `run_id`, `samples`, and parameter overrides. Document any script overrides here using `script_overrides:`.
  - `runs/<run_id>/steps/` — Run-specific script overwrites (e.g. `step_4.R`). `run_step_N.sh` uses these when present.
  - `runs/<run_id>/run_info.json` — Run metadata, notes, reflection outcomes.

```yaml
# Example script_overrides in config_overrides.yml
script_overrides:
  step_4: "Custom annotation for DIPG; uses runs/<run_id>/steps/step_4.R"
```

## Execution Environment

**Never run pipeline scripts locally.** Always use the helper script:

```bash
bash scripts/run_on_comps0.sh <script>
```

This SSHes to `cchin@comps0` and runs via `srun --partition gpunodes -c 2 --mem=64G -t 60 --pty`. `PIPELINE_RUN_ID` and `PIPELINE_STEP_CONFIG` are forwarded automatically.

Examples:
```bash
bash scripts/run_on_comps0.sh scripts/set_run_id.sh
bash scripts/run_on_comps0.sh run_step_1.sh
PIPELINE_RUN_ID=<run_id> bash scripts/run_on_comps0.sh run_step_4.sh
```

If the repo path on comps0 differs: `REPO_ON_COMPS0=/other/path bash scripts/run_on_comps0.sh run_step_1.sh`

## Step Order

Run **one step at a time**. `run_pipeline.sh` exists only to print this workflow — do not use it to run all steps.

| Step | Name | Notes |
|------|------|-------|
| 1 | QC | R |
| 2 | Integration | R |
| 3 | Clustering | R |
| 4 | Annotation | R — intermediate artifacts only; no final labels |
| 4b | Adjudication | **Agent** writes `step4_adjudication_labels.csv` + `step4_adjudication_report.json` |
| 5 | Export | R — fails if adjudication files missing |
| 6 | Load samples | Python |
| 7 | Preprocess | Python |
| 8 | LR scoring | Python |
| 9 | CCI | Python |
| 10 | Aggregate & rank | Python |

## Step 4 → Adjudication → Step 5

**Step 4** writes intermediate outputs only: `step4_collated_candidates.json`, `step4_enrichr_*.csv`, `step4_enrichr_marker_source.csv`, `step4_annotation.pdf`. It does **not** assign final cell-type names.

**Reflecting agent (after step 4, required before step 5):**

1. Read merged config for dataset identity (`species`, `disease`, `dataset_name`).
2. Open `step4_collated_candidates.json` and raw `step4_enrichr_*.csv`.
3. Apply 3-tier reasoning (see [step_specs.md](step_specs.md) Step 4b and [reflection_guidelines.md](reflection_guidelines.md) Checkpoint 2).
4. Write both artifacts under `OUT_DIR`:
   - **`step4_adjudication_labels.csv`** — final labels per cluster with `label_evidence` citing marker genes and/or EnrichR terms.
   - **`step4_adjudication_report.json`** — per-cluster `rationale`, `tier_rationale`, `needs_review`, audit fields.

**Step 5** hard-fails if either file is missing. Validate with:
```bash
bash run_step4b_validate.sh
# or
python scripts/validate_step4b_artifacts.py --config <merged_or_base_config>
```

Record model id / prompt version in `step4_adjudication_report.json` when an LLM is used.

## Execution Loop

1. **Start of run (once):** `bash scripts/run_on_comps0.sh scripts/set_run_id.sh` — creates `runs/<run_id>/` with `config_overrides.yml`, `run_info.json`, `logs/`.
2. Read pipeline docs and dataset identity from merged config (`species`, `disease`, `dataset_name`).
3. **Before step 4:** if `disease` is set in merged config, write `outputs/<run_id>/step4_disease_marker_prior.json`.
4. Run one step: `PIPELINE_RUN_ID=<run_id> bash scripts/run_on_comps0.sh run_step_N.sh`
5. Inspect outputs in `outputs/<run_id>/` and logs in `runs/<run_id>/logs/step_N.log`.

   > **RDS buffering:** R steps may signal completion before the `.rds` file is fully written. Poll file size every 2–3 seconds until stable for 2–3 checks before starting the consuming step.

6. **Reflection (mandatory)** at checkpoints — see [reflection_guidelines.md](reflection_guidelines.md):
   - After steps 1–2: QC effectiveness
   - After step 4: annotation labels (must produce adjudication files before step 5)
   - After step 10: LR pair sanity
7. If checks pass, proceed to next step. If not: apply run-specific changes in `runs/<run_id>/config_overrides.yml`, re-run relevant steps, do not advance until reflection passes.
8. Update `runs/<run_id>/run_info.json` with `steps_completed`, `notes`, and `reflection_outcome` per checkpoint.

## Run Directory Layout

```
runs/<run_id>/
  config_overrides.yml   # run_id, samples, parameter overrides, script_overrides
  config_merged.yml      # auto-generated; do not edit
  run_info.json          # agent-maintained metadata
  logs/
    step_1.log … step_10.log
  steps/                 # optional run-specific script overwrites

outputs/<run_id>/
  step4_disease_marker_prior.json   # agent writes before step 4
  step4_collated_candidates.json    # written by step 4
  step4_adjudication_labels.csv     # agent writes (step 4b)
  step4_adjudication_report.json    # agent writes (step 4b)
  cell_type_metadata.csv            # written by step 5
  GROUND_TRUTH_lr_pairs_ranked.csv  # written by step 10
  ... (rds, h5ad, csv, pdf per step)
```

## File Restrictions

- **Do not edit** `config.yml` or `separated_steps/*` for run-specific behaviour.
- **Do not modify** `.data/*` (raw data) or `evaluation/` (if present).
- Run-specific tuning → `runs/<run_id>/config_overrides.yml`
- Run-specific code → `runs/<run_id>/steps/` (document in `config_overrides.yml`)

## Supporting Documents

- [pipeline_overview.md](pipeline_overview.md) — end-to-end data flow
- [step_specs.md](step_specs.md) — per-step inputs, outputs, libraries
- [reflection_guidelines.md](reflection_guidelines.md) — mandatory quality checkpoints
- [adjudication_artifacts.md](adjudication_artifacts.md) — CSV/JSON schemas and validation checklist
- [step4_disease_marker_prior.example.json](step4_disease_marker_prior.example.json) — disease prior template
