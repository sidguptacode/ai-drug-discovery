#!/usr/bin/env bash
# Run pipeline step 4 (Annotation). Uses config.yml unless PIPELINE_STEP_CONFIG is set. Execute from repo root.

set -e
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

# Source conda so merge + prior script use env Python (PyYAML, etc.); set CONDA_SH on cluster if needed
for f in "${CONDA_SH:-}" "/w/20252/cchin/miniconda3/etc/profile.d/conda.sh" "$HOME/miniconda3/etc/profile.d/conda.sh" "$HOME/anaconda3/etc/profile.d/conda.sh" "/opt/conda/etc/profile.d/conda.sh"; do
  [[ -n "$f" && -r "$f" ]] && source "$f" && break
done
if command -v conda &>/dev/null; then
  conda activate ai-drug-discovery
fi
PYTHON=$(command -v python 2>/dev/null || command -v python3)

if [[ -n "${PIPELINE_RUN_ID:-}" ]] && [[ -f "$REPO_ROOT/runs/$PIPELINE_RUN_ID/config_overrides.yml" ]]; then
  MERGED="$REPO_ROOT/runs/$PIPELINE_RUN_ID/config_merged.yml"
  "$PYTHON" "$REPO_ROOT/scripts/merge_config.py" "$REPO_ROOT/config.yml" "$REPO_ROOT/runs/$PIPELINE_RUN_ID/config_overrides.yml" "$MERGED" --run-id "$PIPELINE_RUN_ID"
  CONFIG="$MERGED"
else
  CONFIG="${PIPELINE_STEP_CONFIG:-$REPO_ROOT/config.yml}"
fi
if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: Config not found: $CONFIG (set PIPELINE_STEP_CONFIG if needed)" >&2
  exit 1
fi
export PIPELINE_STEP_CONFIG="$CONFIG"
RUN_ID="${PIPELINE_RUN_ID:-$(bash "$REPO_ROOT/scripts/get_run_id.sh" "$CONFIG")}"
LOG_FILE=""
[[ -n "$RUN_ID" ]] && LOG_DIR="$REPO_ROOT/runs/$RUN_ID/logs" && mkdir -p "$LOG_DIR" && LOG_FILE="$LOG_DIR/step_4.log"

# Run-specific script override: separated_steps/ is static; use runs/<run_id>/steps/ when present
STEP_SCRIPT="separated_steps/step_4.R"
[[ -n "$RUN_ID" ]] && [[ -f "$REPO_ROOT/runs/$RUN_ID/steps/step_4.R" ]] && STEP_SCRIPT="$REPO_ROOT/runs/$RUN_ID/steps/step_4.R"

echo "====== Step 4 — Annotation (intermediate) ======"
# When config `disease` is set, OUT_DIR/step4_disease_marker_prior.json must already exist (written by the pipeline agent; see skill).
# After step 4: reflecting agent writes step4_adjudication_labels.csv + step4_adjudication_report.json; then run_step4b_validate.sh; then step 5.

if [[ -n "$LOG_FILE" ]]; then
  Rscript "$STEP_SCRIPT" 2>&1 | tee "$LOG_FILE"
else
  Rscript "$STEP_SCRIPT"
fi
