#!/usr/bin/env bash
# Run pipeline step 10 (Aggregate and rank). Uses config.yml unless PIPELINE_STEP_CONFIG is set. Execute from repo root.

set -e
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

CONFIG="${PIPELINE_STEP_CONFIG:-$REPO_ROOT/config.yml}"
if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: Config not found: $CONFIG (set PIPELINE_STEP_CONFIG if needed)" >&2
  exit 1
fi
export PIPELINE_STEP_CONFIG="$CONFIG"
RUN_ID="${PIPELINE_RUN_ID:-$(bash "$REPO_ROOT/scripts/get_run_id.sh" "$CONFIG")}"
LOG_FILE=""
[[ -n "$RUN_ID" ]] && LOG_DIR="$REPO_ROOT/runs/$RUN_ID/logs" && mkdir -p "$LOG_DIR" && LOG_FILE="$LOG_DIR/step_10.log"

for f in "${CONDA_SH:-}" "/w/20251/sarya/miniconda3/etc/profile.d/conda.sh" "$HOME/miniconda3/etc/profile.d/conda.sh" "$HOME/anaconda3/etc/profile.d/conda.sh" "/opt/conda/etc/profile.d/conda.sh"; do
  [[ -n "$f" && -r "$f" ]] && source "$f" && break
done
if command -v conda &>/dev/null; then
  conda activate ai-drug-discovery
fi
PYTHON=$(command -v python 2>/dev/null || command -v python3)

echo "====== Step 10 — Aggregate and rank ======"
if [[ -n "$LOG_FILE" ]]; then
  "$PYTHON" separated_steps/step_10.py 2>&1 | tee "$LOG_FILE"
else
  "$PYTHON" separated_steps/step_10.py
fi
