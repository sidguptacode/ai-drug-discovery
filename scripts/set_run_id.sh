#!/usr/bin/env bash
# Start a new pipeline run: generate UUID, create runs/<uuid>/ with config_overrides.yml and run_info.json.
# Does NOT modify config.yml (static). Run from repo root. Uses PIPELINE_STEP_CONFIG for base config path.

set -e
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"
CONFIG="${PIPELINE_STEP_CONFIG:-$REPO_ROOT/config.yml}"
export PIPELINE_STEP_CONFIG="$CONFIG"
RUN_ID=$(python3 scripts/set_run_id.py)
echo "Run ID: $RUN_ID"
echo "Logs: runs/$RUN_ID/logs/"
echo "Run info: runs/$RUN_ID/run_info.json"
echo "Run steps with: PIPELINE_RUN_ID=$RUN_ID bash scripts/run_on_comps0.sh run_step_N.sh"
