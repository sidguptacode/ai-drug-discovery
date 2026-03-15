#!/usr/bin/env bash
# Start a new pipeline run: generate UUID, set run_id in config.yml, create runs/<uuid>/
# Run from repo root. Uses PIPELINE_STEP_CONFIG for config path.

set -e
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"
CONFIG="${PIPELINE_STEP_CONFIG:-$REPO_ROOT/config.yml}"
export PIPELINE_STEP_CONFIG="$CONFIG"
RUN_ID=$(python3 scripts/set_run_id.py)
echo "Run ID: $RUN_ID"
echo "Logs: runs/$RUN_ID/logs/"
echo "Run info: runs/$RUN_ID/run_info.json"
