#!/usr/bin/env bash
# Run pipeline step 6 (Load samples). Uses config.yml unless PIPELINE_STEP_CONFIG is set. Execute from repo root.

set -e
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

if [[ -n "${PIPELINE_RUN_ID:-}" ]] && [[ -f "$REPO_ROOT/runs/$PIPELINE_RUN_ID/config_overrides.yml" ]]; then
  MERGED="$REPO_ROOT/runs/$PIPELINE_RUN_ID/config_merged.yml"
  python3 "$REPO_ROOT/scripts/merge_config.py" "$REPO_ROOT/config.yml" "$REPO_ROOT/runs/$PIPELINE_RUN_ID/config_overrides.yml" "$MERGED" --run-id "$PIPELINE_RUN_ID"
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
[[ -n "$RUN_ID" ]] && LOG_DIR="$REPO_ROOT/runs/$RUN_ID/logs" && mkdir -p "$LOG_DIR" && LOG_FILE="$LOG_DIR/step_6.log"

# Source conda so it is available in non-interactive/srun environments (set CONDA_SH if needed on cluster)
for f in "${CONDA_SH:-}" "/w/20251/sarya/miniconda3/etc/profile.d/conda.sh" "$HOME/miniconda3/etc/profile.d/conda.sh" "$HOME/anaconda3/etc/profile.d/conda.sh" "/opt/conda/etc/profile.d/conda.sh"; do
  [[ -n "$f" && -r "$f" ]] && source "$f" && break
done
if command -v conda &>/dev/null; then
  conda activate ai-drug-discovery
fi
PYTHON=$(command -v python 2>/dev/null || command -v python3)

# Run-specific script override: separated_steps/ is static; use runs/<run_id>/steps/ when present
STEP_SCRIPT="separated_steps/step_6.py"
[[ -n "$RUN_ID" ]] && [[ -f "$REPO_ROOT/runs/$RUN_ID/steps/step_6.py" ]] && STEP_SCRIPT="$REPO_ROOT/runs/$RUN_ID/steps/step_6.py"

echo "====== Step 6 — Load samples (Python) ======"
if [[ -n "$LOG_FILE" ]]; then
  "$PYTHON" "$STEP_SCRIPT" 2>&1 | tee "$LOG_FILE"
else
  "$PYTHON" "$STEP_SCRIPT"
fi
