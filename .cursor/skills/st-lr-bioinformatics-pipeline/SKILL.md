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

**Run one step at a time** (from repo root). Do not run all steps in a row. Use the per-step scripts:

- `bash run_step_1.sh` — QC
- `bash run_step_2.sh` — Integration
- `bash run_step_3.sh` — Clustering
- `bash run_step_4.sh` — Annotation
- `bash run_step_5.sh` — Export
- `bash run_step_6.sh` — Load samples (Python)
- `bash run_step_7.sh` — Preprocess
- `bash run_step_8.sh` — LR scoring
- `bash run_step_9.sh` — CCI
- `bash run_step_10.sh` — Aggregate and rank

Run steps in order (each step expects the previous step’s outputs). After each step, check outputs and decide whether to proceed or re-run (see Execution loop).

**Note**: `run_pipeline.sh` exists only to print this workflow; it exits with an error and must not be used to execute all steps in sequence.

### Run ID and logging

At the **start of a run**, create a run ID so all step logs and run-level outputs are grouped:

```bash
bash scripts/set_run_id.sh
```

This generates a UUID, sets `run_id` in `config.yml`, creates `runs/<run_id>/` with `logs/` and `run_info.json`, and prints the run ID. Each `run_step_N.sh` then writes that step’s stdout/stderr to `runs/<run_id>/logs/step_N.log` (and still prints to the terminal). Pipeline artifacts (rds, h5ad, csv, pdf) continue to go to `out_dir` from config; use `runs/<run_id>/` for logs and for any run-level notes the agent wants to record.

- **runs/<run_id>/logs/step_1.log … step_10.log** — Per-step execution logs (R/Python stdout and stderr).
- **runs/<run_id>/run_info.json** — Run metadata (`run_id`, `started_at`, `steps_completed`, `notes`). The agent can update this (e.g. append completed steps, add notes) after each step or at the end of the run.

If `config.yml` has no `run_id` (or it is null), step scripts still run but do not write to a log file; you can set `run_id` by running `scripts/set_run_id.sh` at the start of a run.

## Execution Loop

Do not run all steps in a row. For each step:

1. **Start of run (once)**: Run `bash scripts/set_run_id.sh` to set `run_id` in config and create `runs/<run_id>/` for logs and run_info.json.
2. Read the pipeline documentation (overview, step_specs, reflection_guidelines) as needed.
3. Run **one** step (e.g. `bash run_step_N.sh`). If a run_id is set, that step’s output is also written to `runs/<run_id>/logs/step_N.log`.
4. Inspect that step’s outputs in `out_dir` (from config) and, if logging, in `runs/<run_id>/logs/step_N.log`.
5. Evaluate whether the results are reasonable (use reflection_guidelines.md at the relevant checkpoints: after 1–2, after 4, after 8–10).
6. If reasonable, proceed to the **next** step. If not, re-run the same step (e.g. after changing `config.yml` or code in `separated_steps/`) or fix and re-run before advancing. Optionally update `runs/<run_id>/run_info.json` (e.g. `steps_completed`, `notes`).
7. Repeat: run one step → check outputs → reason about results → next step or re-run. Never batch-run the full pipeline without evaluation between steps.

## File Restrictions

- **Editable**: `separated_steps/*` (pipeline code). Parameter-only changes in `config.yml` are allowed when tuning QC, clustering, annotation, or LR/CCI settings.
- **Read-only**: `.data/*` (raw data). Do not modify. If `evaluation/` exists, treat it as read-only.

## Additional Resources

- Step-by-step specs and required outputs: [step_specs.md](step_specs.md)
- Pipeline flow and key steps: [pipeline_overview.md](pipeline_overview.md)
- Reflection checkpoints: [reflection_guidelines.md](reflection_guidelines.md)
