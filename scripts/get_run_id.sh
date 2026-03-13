#!/usr/bin/env bash
# Print run_id from config. Usage: get_run_id.sh [config_path]
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${1:-${PIPELINE_STEP_CONFIG:-$REPO_ROOT/config.yml}}"
python3 "$(dirname "${BASH_SOURCE[0]}")/get_run_id.py" "$CONFIG" 2>/dev/null || true
