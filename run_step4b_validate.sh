#!/usr/bin/env bash
# Validate Step 4b adjudication artifacts (CSV + JSON). Does not run adjudication — agent-led.
# Uses merged config when PIPELINE_RUN_ID is set (same pattern as run_step_4.sh).

set -e
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"
PYTHON=$(command -v python 2>/dev/null || command -v python3)

if [[ -n "${PIPELINE_RUN_ID:-}" ]] && [[ -f "$REPO_ROOT/runs/$PIPELINE_RUN_ID/config_overrides.yml" ]]; then
  MERGED="$REPO_ROOT/runs/$PIPELINE_RUN_ID/config_merged.yml"
  "$PYTHON" "$REPO_ROOT/scripts/merge_config.py" "$REPO_ROOT/config.yml" "$REPO_ROOT/runs/$PIPELINE_RUN_ID/config_overrides.yml" "$MERGED" --run-id "$PIPELINE_RUN_ID"
  CONFIG="$MERGED"
else
  CONFIG="${PIPELINE_STEP_CONFIG:-$REPO_ROOT/config.yml}"
fi
if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: Config not found: $CONFIG" >&2
  exit 1
fi

echo "====== Step 4b — Validate adjudication artifacts ======"
exec "$PYTHON" "$REPO_ROOT/scripts/validate_step4b_artifacts.py" --config "$CONFIG" "$@"
