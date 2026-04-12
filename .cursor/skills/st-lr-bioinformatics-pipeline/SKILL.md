---
name: st-lr-bioinformatics-pipeline
description: Run and iteratively improve the spatial transcriptomics ligand-receptor pipeline (separated_steps). Use when running the ST-LR pipeline, debugging steps, or interpreting QC, annotation, or LR pair outputs.
---

# ST-LR Bioinformatics Pipeline Skill

## Initialization

Before running or changing the pipeline, read the following documents (from this skill directory):

1. [pipeline_overview.md](pipeline_overview.md)
2. [step_specs.md](step_specs.md)
3. [reflection_guidelines.md](reflection_guidelines.md)
4. [adjudication_artifacts.md](adjudication_artifacts.md) — when performing or reviewing Step 4b (annotation adjudication).

Items 1–3 describe the pipeline structure, step definitions, and how to reason about outputs. Item 4 defines adjudication CSV/JSON schemas and validation.

**Disease prior (Step 4) — controller agent obligation:** There is **no** automated Open Targets / GraphQL step. The **controller agent** must create **`OUT_DIR/step4_disease_marker_prior.json`** (under `outputs/<run_id>/` when using `PIPELINE_RUN_ID`) **before** `run_step_4.sh` whenever top-level **`disease`** is set in merged config. Step 4 **stops** if that file is missing or if **`genes`** is absent or empty.

When curating **`genes`** (HGNC symbols for the species in config), reason explicitly in three buckets — **most weight on (3)**:

1. **Disease-linked** — markers and pathways directly associated with the configured disease.
2. **Locus / anatomical context** — genes reflecting expected tissue or regional identity (e.g. organ- or niche-relevant markers) aligned with the dataset.
3. **Canonical cell types** — established markers for cell types you expect in this kind of spatial transcriptomics analysis (broad coverage across expected lineages); this bucket should carry the **largest** share of the list.

Use your knowledge and **web search** as needed to verify symbols and marker lists. Optional JSON fields (`disease_prior_tissues`, `disease_prior_cell_types`, `notes`, `source`, etc.) are for audit only; **marker filtering uses `genes` only.** Shape: [step4_disease_marker_prior.example.json](step4_disease_marker_prior.example.json). Override path: **`annotation.disease_marker_prior_json`** (relative to `out_dir` or absolute).

## Static files vs run-specific edits

- **Static (canonical; do not edit for a specific run)**  
  - **config.yml** — Base dataset identity, paths, and default parameters. Keep `run_id: null` and generic defaults here.  
  - **separated_steps/** — **Canonical pipeline code.** This is the reference implementation. Do not rewrite or edit these files to implement run-specific behaviour (e.g. a different annotation step or a step done separately). Treat them as read-only for run-specific changes.

- **Run-specific (edits that depend on context/data/outcomes for that run)**  
  - **runs/<run_id>/config_overrides.yml** — Created when the run is started. Holds `run_id`, `samples`, and any parameter overrides. Steps use **merged** config (base + this file) when `PIPELINE_RUN_ID` is set. **Must also document any run-specific script overrides** (see below).  
  - **runs/<run_id>/steps/** — **Run-specific script overwrites.** If a step needs different logic for this run, place the modified script here (e.g. `runs/<run_id>/steps/step_4.R`). `run_step_N.sh` uses `runs/<run_id>/steps/step_N.*` when present instead of `separated_steps/step_N.*`. Do not edit `separated_steps/` for that run—write the override here and document it in `config_overrides.yml`.  
  - **runs/<run_id>/run_info.json** — Run metadata, notes, reflection outcomes. Safe to update during/after the run.

**Documenting script overrides in config_overrides.yml**  
When a run uses scripts in `runs/<run_id>/steps/`, add a `script_overrides` section to that run’s `config_overrides.yml` so the run is self-describing. Example:

```yaml
# In runs/<run_id>/config_overrides.yml
script_overrides:
  step_4: "Custom annotation for DIPG; uses runs/<run_id>/steps/step_4.R"
  step_6: "Separate load step for this cohort; see steps/step_6.py"
```

This does not drive behaviour (the runner uses file presence under `runs/<run_id>/steps/`); it explains which steps are overridden and why.

Do **not** change `config.yml` or `separated_steps/*` to tune for one run. Put run-specific parameters in `runs/<run_id>/config_overrides.yml`, run-specific code in `runs/<run_id>/steps/`, and document overrides in `config_overrides.yml`. Run steps with `PIPELINE_RUN_ID=<run_id>`.

## Core Commands

**Config**: Base config is `config.yml` at repo root. When `PIPELINE_RUN_ID` is set and `runs/<run_id>/config_overrides.yml` exists, steps use the **merged** config (base + overrides). Override the single config path with `PIPELINE_STEP_CONFIG` if needed.

**Environment**: R steps (1–5) run with the default R. Python steps (6–10) use the conda env `ai-drug-discovery` (their scripts activate it).

### Where to run (mandatory)

**Do not run pipeline scripts locally.** Whenever the agent runs any of the internal scripts (`scripts/set_run_id.sh`, `run_step_1.sh` … `run_step_10.sh`), it must:

1. **SSH to comps0**: `ssh cchin@comps0`
2. **Run the script via srun** with: `--partition gpunodes -c 2 --mem=64G -t 60 --pty`

Use the helper script from the **local** repo (so the agent runs this locally):

```bash
bash scripts/run_on_comps0.sh <script>
```

Examples:

- `bash scripts/run_on_comps0.sh scripts/set_run_id.sh`
- `bash scripts/run_on_comps0.sh run_step_1.sh`
- For a **specific run** (use merged config and that run’s logs): `PIPELINE_RUN_ID=<run_id> bash scripts/run_on_comps0.sh run_step_4.sh`

This SSH’s to `cchin@comps0`, then runs `srun --partition gpunodes -c 2 --mem=64G -t 60 --pty bash <script>` in the project directory on comps0. `PIPELINE_RUN_ID` and `PIPELINE_STEP_CONFIG` are forwarded to the remote when set. If the repo path on comps0 differs, set `REPO_ON_COMPS0` before calling (e.g. `REPO_ON_COMPS0=/other/path bash scripts/run_on_comps0.sh run_step_1.sh`).

**Run one step at a time** (via the helper above). Do not run all steps in a row. Steps in order:

- Step 1 — QC
- Step 2 — Integration
- Step 3 — Clustering
- Step 4 — Annotation (intermediate artifacts only; no final cell-type labels in R)
- Step 4b — Mechanical collation is in step 4 (`step4_collated_candidates.json`); **reflecting agent** writes final labels + report (mandatory before step 5)
- Step 5 — Export (fails if adjudication files are missing)
- Step 6 — Load samples (Python)
- Step 7 — Preprocess
- Step 8 — LR scoring
- Step 9 — CCI
- Step 10 — Aggregate and rank

After each step, check outputs and decide whether to proceed or re-run (see Execution loop).

**Note**: `run_pipeline.sh` exists only to print this workflow; it exits with an error and must not be used to execute all steps in sequence.

### Step 4 → reflecting agent → step 5 (mandatory adjudication path)

**Step 4** (`separated_steps/step_4.R`) writes **intermediate** outputs only: marker tables, `step4_enrichr_<community>_<db>.csv`, `step4_enrichr_marker_source.csv`, **`step4_collated_candidates.json`** (mechanical top-k collation per cluster and DB), and `step4_annotation.pdf`. It does **not** assign final cell-type names from EnrichR; `cell_type_label` in `step4_seurat_annotated.rds` is set to the **community id** until step 5.

**Reflecting agent (controller), after step 4 — required before step 5**:

1. Read merged config for **dataset identity** (`species`, `disease`, `dataset_name`) and open **`step4_collated_candidates.json`** plus raw `step4_enrichr_*.csv` as needed.
2. Apply **3-tier** reasoning ([step_specs.md](step_specs.md) Step 4b, [reflection_guidelines.md](reflection_guidelines.md) Checkpoint 2).
3. Write **both** artifacts under **`OUT_DIR`** (see [adjudication_artifacts.md](adjudication_artifacts.md)):
   - **`step4_adjudication_labels.csv`** — final labels per cluster, required **`label_evidence`** column (plain text citing marker genes and/or EnrichR terms from step 4).
   - **`step4_adjudication_report.json`** — per-cluster **`rationale`** and **`tier_rationale`** (grounded in markers/`key_marker_genes` and/or **`candidates`**), `needs_review`, audit fields.

**Step 5** **requires** those two files. If either is missing, it stops with: **“Step 4 final outputs not completed”**. It applies **`adjudication_label`** to **all** clusters (not conditional on `override_applied`; that flag is audit-only).

**Audit**: Record model id / prompt version in `step4_adjudication_report.json` when an LLM is used; collation inputs are reproducible from `step4_collated_candidates.json` and EnrichR CSVs.

**Validation helper**: Run `scripts/validate_step4b_artifacts.py --config <merged_or_base_config>` (or `run_step4b_validate.sh`). Optional: `--forbid-regex`, `--require-approval` + `step4_adjudication_approved.flag`. The validator does not pick labels.

### Run ID and logging

At the **start of a run**, create a run ID (run on comps0 via the helper):

```bash
bash scripts/run_on_comps0.sh scripts/set_run_id.sh
```

This generates a UUID, creates `runs/<run_id>/` with `logs/`, `config_overrides.yml` (run_id + samples from base), and `run_info.json`. It does **not** modify `config.yml`. To have step logs and merged config for this run, run steps with **PIPELINE_RUN_ID** set:

```bash
PIPELINE_RUN_ID=<run_id> bash scripts/run_on_comps0.sh run_step_1.sh
```

Each `run_step_N.sh` then merges `config.yml` with `runs/<run_id>/config_overrides.yml` (and sets `out_dir` to **outputs/<run_id>** for that run), uses that merged config, writes stdout/stderr to `runs/<run_id>/logs/step_N.log`, and writes pipeline artifacts (rds, h5ad, csv, pdf) to **outputs/<run_id>/**.

- **outputs/<run_id>/** — Pipeline outputs for this run (step RDS/h5ad/CSV/PDF). When `PIPELINE_RUN_ID` is set, the merged config sets `out_dir` to this path. Other folders under `outputs/` (e.g. `ovarian_v1_visium_format`) are not run-specific.
- **runs/<run_id>/config_overrides.yml** — Run-specific overrides (run_id, samples, annotation filters, etc.) and optional `script_overrides` describing any scripts in `steps/`. Edit this for that run; do not edit `config.yml` for run-specific tuning.
- **runs/<run_id>/steps/** — Optional. Run-specific script overwrites (e.g. `step_4.R`, `step_6.py`). When present, `run_step_N.sh` uses these instead of `separated_steps/`. Document in `config_overrides.yml` with `script_overrides`.
- **runs/<run_id>/config_merged.yml** — Generated by step scripts when `PIPELINE_RUN_ID` is set (base + overrides). Do not edit by hand.
- **runs/<run_id>/logs/step_1.log … step_10.log** — Per-step execution logs.
- **runs/<run_id>/run_info.json** — Run metadata; the agent can update this after each step or at the end of the run.

If you run steps **without** `PIPELINE_RUN_ID`, the single config is `config.yml` (and no step log file is written unless `config.yml` has a non-null `run_id` from legacy use).

## Execution Loop

Do not run all steps in a row. For each step:

1. **Start of run (once)**: Run `bash scripts/run_on_comps0.sh scripts/set_run_id.sh` to create `runs/<run_id>/` (config_overrides.yml, run_info.json, logs/). Do not edit `config.yml` for this run.
2. Read the pipeline documentation (overview, step_specs, reflection_guidelines) as needed. **Before any reflection**, read **dataset identity** from the effective config (merged config when using `PIPELINE_RUN_ID`, else `config.yml`): `species`, `disease`, `dataset_name`.
3. Run **one** step via comps0 with that run’s config: `PIPELINE_RUN_ID=<run_id> bash scripts/run_on_comps0.sh run_step_N.sh`. That step’s output is written to `runs/<run_id>/logs/step_N.log`. **Before step 4**, if merged config has non-empty `disease`, write `outputs/<run_id>/step4_disease_marker_prior.json` (see disease-prior obligation above); step 4 will not run without it.
4. Inspect that step’s outputs in `outputs/<run_id>/` (when using `PIPELINE_RUN_ID`) or `out_dir` from config, and step log in `runs/<run_id>/logs/step_N.log`.
5. **Reflection (mandatory)** at the relevant checkpoint (after 1–2, after 4, after 8–10): Open [reflection_guidelines.md](reflection_guidelines.md) and perform **every** check for that checkpoint. At Checkpoint 2, after step 4, the agent must **produce adjudication CSV + JSON** (final labels) consistent with dataset identity; step 5 will not run without them. If labels contradict dataset identity, remediate (config + re-run step 4, or revised adjudication files) before step 5.
6. If all checkpoint checks pass, proceed to the **next** step. If any check fails, apply **run-specific** changes in `runs/<run_id>/config_overrides.yml` (e.g. annotation filters), or shared code changes in `separated_steps/*`; do not edit `config.yml` for run-specific tuning. Re-run the relevant step(s) with `PIPELINE_RUN_ID=<run_id>` and do not advance until the re-run passes reflection. Update `runs/<run_id>/run_info.json` with `steps_completed`, `notes`, and optionally a short `reflection_outcome` per checkpoint (e.g. `passed`, `action`, `reason`).
7. Repeat: run one step → check outputs → perform full reflection for that checkpoint → next step or re-run. Never batch-run the full pipeline without evaluation between steps.

## File Restrictions

- **Static / canonical (do not edit for run-specific behaviour)**: `config.yml` (base config), `separated_steps/*` (pipeline code). For run-specific parameter tuning, edit `runs/<run_id>/config_overrides.yml`. For run-specific step logic, add scripts under `runs/<run_id>/steps/` and document them in `config_overrides.yml` (`script_overrides`); do not modify `separated_steps/*` for that run.
- **Run-specific (editable for that run)**: `runs/<run_id>/config_overrides.yml`, `runs/<run_id>/steps/*`, `runs/<run_id>/run_info.json`.
- **Read-only**: `.data/*` (raw data). Do not modify. If `evaluation/` exists, treat it as read-only.

## Additional Resources

- Step-by-step specs and required outputs: [step_specs.md](step_specs.md)
- Pipeline flow and key steps: [pipeline_overview.md](pipeline_overview.md)
- Reflection checkpoints: [reflection_guidelines.md](reflection_guidelines.md)
- Adjudication artifact schemas and validation: [adjudication_artifacts.md](adjudication_artifacts.md)
