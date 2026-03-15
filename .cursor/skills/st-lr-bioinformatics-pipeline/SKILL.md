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

These describe the pipeline structure, step definitions, and how to reason about outputs.

## Core Commands

**Config**: `config.yml` at repo root. Override with `PIPELINE_STEP_CONFIG` if using a different path.

**Environment**: R steps (1–5) run with the default R. Python steps (6–10) use the conda env `ai-drug-discovery` (their scripts activate it).

### Where to run (mandatory)

**Do not run pipeline scripts locally.** Whenever the agent runs any of the internal scripts (`scripts/set_run_id.sh`, `run_step_1.sh` … `run_step_10.sh`), it must:

1. **SSH to comps0**: `ssh sarya@comps0`
2. **Run the script via srun** with: `--partition gpunodes -c 2 --mem=64G -t 60 --pty`

Use the helper script from the **local** repo (so the agent runs this locally):

```bash
bash scripts/run_on_comps0.sh <script>
```

Examples:

- `bash scripts/run_on_comps0.sh scripts/set_run_id.sh`
- `bash scripts/run_on_comps0.sh run_step_1.sh`
- `bash scripts/run_on_comps0.sh run_step_4.sh`

This SSH’s to `sarya@comps0`, then runs `srun --partition gpunodes -c 2 --mem=64G -t 60 --pty bash <script>` in the project directory on comps0. If the repo path on comps0 differs, set `REPO_ON_COMPS0` before calling (e.g. `REPO_ON_COMPS0=/other/path bash scripts/run_on_comps0.sh run_step_1.sh`).

**Run one step at a time** (via the helper above). Do not run all steps in a row. Steps in order:

- Step 1 — QC
- Step 2 — Integration
- Step 3 — Clustering
- Step 4 — Annotation
- Step 5 — Export
- Step 6 — Load samples (Python)
- Step 7 — Preprocess
- Step 8 — LR scoring
- Step 9 — CCI
- Step 10 — Aggregate and rank

After each step, check outputs and decide whether to proceed or re-run (see Execution loop).

**Note**: `run_pipeline.sh` exists only to print this workflow; it exits with an error and must not be used to execute all steps in sequence.

### Run ID and logging

At the **start of a run**, create a run ID so all step logs and run-level outputs are grouped (run on comps0 via the helper):

```bash
bash scripts/run_on_comps0.sh scripts/set_run_id.sh
```

This generates a UUID, sets `run_id` in `config.yml`, creates `runs/<run_id>/` with `logs/` and `run_info.json`, and prints the run ID. Each `run_step_N.sh` then writes that step’s stdout/stderr to `runs/<run_id>/logs/step_N.log` (and still prints to the terminal). Pipeline artifacts (rds, h5ad, csv, pdf) continue to go to `out_dir` from config; use `runs/<run_id>/` for logs and for any run-level notes the agent wants to record.

- **runs/<run_id>/logs/step_1.log … step_10.log** — Per-step execution logs (R/Python stdout and stderr).
- **runs/<run_id>/run_info.json** — Run metadata (`run_id`, `started_at`, `steps_completed`, `notes`). The agent can update this (e.g. append completed steps, add notes) after each step or at the end of the run.

If `config.yml` has no `run_id` (or it is null), step scripts still run but do not write to a log file; you can set `run_id` by running `scripts/set_run_id.sh` at the start of a run.

## Execution Loop

Do not run all steps in a row. For each step:

1. **Start of run (once)**: Run `bash scripts/run_on_comps0.sh scripts/set_run_id.sh` to set `run_id` in config and create `runs/<run_id>/` for logs and run_info.json.
2. Read the pipeline documentation (overview, step_specs, reflection_guidelines) as needed. **Before any reflection**, read `config.yml` for **dataset identity** (`species`, `disease`, `dataset_name`) and use these when evaluating outputs.
3. Run **one** step via comps0: `bash scripts/run_on_comps0.sh run_step_N.sh`. If a run_id is set, that step’s output is also written to `runs/<run_id>/logs/step_N.log`.
4. Inspect that step’s outputs in `out_dir` (from config) and, if logging, in `runs/<run_id>/logs/step_N.log`.
5. **Reflection (mandatory)** at the relevant checkpoint (after 1–2, after 4, after 8–10): Open [reflection_guidelines.md](reflection_guidelines.md) and perform **every** check for that checkpoint. At Checkpoint 2, ensure annotation labels are **consistent with dataset identity** (species and, where applicable, tissue/region implied by disease or sample); if any label contradicts that, remediate and re-run step 4 before proceeding.
6. If all checkpoint checks pass, proceed to the **next** step. If any check fails, apply config or code changes, re-run the relevant step(s), and do not advance until the re-run passes reflection. Update `runs/<run_id>/run_info.json` with `steps_completed`, `notes`, and optionally a short `reflection_outcome` per checkpoint (e.g. `passed`, `action`, `reason`).
7. Repeat: run one step → check outputs → perform full reflection for that checkpoint → next step or re-run. Never batch-run the full pipeline without evaluation between steps.

## File Restrictions

- **Editable**: `separated_steps/*` (pipeline code). Parameter-only changes in `config.yml` are allowed when tuning QC, clustering, annotation, or LR/CCI settings.
- **Read-only**: `.data/*` (raw data). Do not modify. If `evaluation/` exists, treat it as read-only.

## Additional Resources

- Step-by-step specs and required outputs: [step_specs.md](step_specs.md)
- Pipeline flow and key steps: [pipeline_overview.md](pipeline_overview.md)
- Reflection checkpoints: [reflection_guidelines.md](reflection_guidelines.md)
