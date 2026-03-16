#!/usr/bin/env bash
# Run pipeline step 4 (Annotation). Uses config.yml unless PIPELINE_STEP_CONFIG is set. Execute from repo root.

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
[[ -n "$RUN_ID" ]] && LOG_DIR="$REPO_ROOT/runs/$RUN_ID/logs" && mkdir -p "$LOG_DIR" && LOG_FILE="$LOG_DIR/step_4.log"

# Run-specific script override: separated_steps/ is static; use runs/<run_id>/steps/ when present
STEP_SCRIPT="separated_steps/step_4.R"
[[ -n "$RUN_ID" ]] && [[ -f "$REPO_ROOT/runs/$RUN_ID/steps/step_4.R" ]] && STEP_SCRIPT="$REPO_ROOT/runs/$RUN_ID/steps/step_4.R"

echo "====== Step 4 — Annotation ======"
if [[ -n "$LOG_FILE" ]]; then
  Rscript "$STEP_SCRIPT" 2>&1 | tee "$LOG_FILE"
else
  Rscript "$STEP_SCRIPT"
fi
