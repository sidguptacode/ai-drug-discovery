#!/usr/bin/env bash
# Run pipeline step 9 (CCI). Uses config.yml unless PIPELINE_STEP_CONFIG is set. Execute from repo root.

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
[[ -n "$RUN_ID" ]] && LOG_DIR="$REPO_ROOT/runs/$RUN_ID/logs" && mkdir -p "$LOG_DIR" && LOG_FILE="$LOG_DIR/step_9.log"

if command -v conda &>/dev/null; then
  eval "$(conda shell.bash hook)"
  conda activate ai-drug-discovery
fi

echo "====== Step 9 — CCI ======"
if [[ -n "$LOG_FILE" ]]; then
  python separated_steps/step_9.py 2>&1 | tee "$LOG_FILE"
else
  python separated_steps/step_9.py
fi
