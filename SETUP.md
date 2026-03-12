# Setup — Pipeline controller and environments

This project runs a single-cell to ligand–receptor pipeline: **steps 1–5 in R** (Seurat), **steps 6–10 in Python** (scanpy, stlearn). The **controller** (`run_pipeline.py`) provides agent-ready tools (run a step, setup run dir from config) and a minimal CLI: at each step you choose **rerun** or **next**; step configs (`step_1.json` … `step_10.json`) are written by you or an agent.

## Python environment

**Option A — venv + pip**

```bash
python -m venv .venv
.venv\Scripts\activate   # Windows
# source .venv/bin/activate   # Linux/macOS
pip install -r requirements.txt
```

**Option B — Conda**

```bash
conda env create -f environment.yml
conda activate ai-drug-discovery
```

The controller only needs the standard library and PyYAML; pipeline steps 6–10 need scanpy, stlearn, pandas, matplotlib.

## R environment

Install R (e.g. from [CRAN](https://cloud.r-project.org/)) and ensure `Rscript` is on your `PATH`.

Required R packages (install from CRAN or Bioconductor as needed):

- **Required:** `Seurat`, `ggplot2`, `patchwork`, `dplyr`, `yaml`, `jsonlite`, `future`
- **Optional:** `harmony` (integration), `clustree` (clustering plots), `enrichR` (cell-type annotation), `glmGamPoi` (faster SCTransform)

In R:

```r
install.packages(c("Seurat", "ggplot2", "patchwork", "dplyr", "yaml", "jsonlite", "future"))
# Optional:
install.packages("harmony")
install.packages("clustree")
install.packages("enrichR")
```

If you use Bioconductor for Seurat, follow the [Seurat installation guide](https://satijalab.org/seurat/articles/install.html).

## Running the pipeline

The controller is **agent-ready**: it exposes tools (run step, setup run dir) and a minimal CLI. You or an agent write/edit `step_1.json` … `step_10.json`; the controller runs steps and offers **rerun** or **next step** at each pause.

You provide a **run_id** (e.g. `my_run_01`). Configs always go to **`runs/<run_id>/`** (step_1.json … step_10.json, last_step.txt); outputs always go to **`outputs/<run_id>/`**.

**1. Create a run and seed step configs from `config.yml` (once):**

```bash
python run_pipeline.py --run-id my_run_01 --init
```

Optional: merge a custom override when initializing:

```bash
python run_pipeline.py --run-id my_run_01 --init --custom-config my_overrides.yml
```

**2. Interactive loop (rerun / next step):**

```bash
python run_pipeline.py --run-id my_run_01
```

After each step you get: `[r]erun step N | [n]ext step?` — choose **r** to run the same step again (e.g. after editing `step_N.json`) or **n** to continue. Config files are **not** overwritten; edit them manually (or have an agent write them) between runs.

**3. Resume from last completed step:**

```bash
python run_pipeline.py --run-id my_run_01 --resume
```

**4. Single-step run (for agents or scripting):**

```bash
python run_pipeline.py --run-id my_run_01 --step 3
```

Runs only step 3 and exits with that step’s exit code.

### Agent-facing tools (import and call)

From Python (e.g. an agent’s tool layer), use:

- **`get_project_root()`** → `Path` — project root for cwd.
- **`run_step(step: int, run_dir: Path | str, project_root=None)`** → `int` — run one step (1–10); returns 0 on success, else subprocess exit code. Reads `run_dir/step_{step}.json`. Use `run_dir = project_root / "runs" / run_id`.
- **`get_step_config_path(step: int, run_dir: Path)`** → `Path` — path to `step_{step}.json` (for reading/writing params).
- **`setup_run_dir_from_config(config_path, run_dir, custom_config_path=None, project_root=None)`** — create run_dir and write `step_1.json` … `step_10.json`; run_dir must be `project_root / "runs" / run_id`. Each step’s `out_dir` is set to `outputs/<run_id>`.

Example agent flow: call `setup_run_dir_from_config` with `run_dir=project_root / "runs" / run_id`, then in a loop call `run_step(n, run_dir)`, read step outputs, write updated `get_step_config_path(n+1, run_dir)` and call `run_step(n+1, run_dir)`, or re-run step n by calling `run_step(n, run_dir)` again.

---

## Agent controller (OpenAI Responses API)

An **AI agent** can drive the pipeline using the OpenAI Responses API. The controller (`agent_controller.py`) manages the run loop locally: it calls `client.responses.create` with tools, executes any tool calls (init run, get/set step config, run step, read outputs, get metadata) on your machine, and sends the results back until the model returns a final text response. No Assistant ID or thread is stored on the server; state is in the request/response chain.

**Requirements:** Use the project conda environment (it includes `openai`). Activate with `conda activate ai-drug-discovery`. Set your API key:

```bash
export OPENAI_API_KEY=sk-...
```

**Run the agent** (from the project root, with the conda env active):

```bash
conda activate ai-drug-discovery
python agent_controller.py --message "Initialize a run named run_01 from config.yml, then run step 1."
```

Optional arguments:

- **`--run-dir`** — Default run directory for tools (e.g. `runs/run_01`). The agent can still pass a different `run_dir` in each tool call.
- **`--model`** — Model name (default: `gpt-4o-mini`).
- **`--max-tool-rounds`** — Maximum tool-call rounds (default: 50).

**Example:**

```bash
conda activate ai-drug-discovery
python agent_controller.py --run-dir runs/run_01 --message "Run step 1 with default config, then tell me how many spots passed QC."
```

The controller uses the **Responses API** only (`client.responses.create`); it does not use the deprecated Assistants (threads/runs) API.
