#!/usr/bin/env bash
# Run a pipeline script on comps0 via SSH + srun. Usage: run_on_comps0.sh <script> [script_arg ...]
# Example: run_on_comps0.sh run_step_1.sh
#          run_on_comps0.sh scripts/set_run_id.sh
# Repo path on comps0: set REPO_ON_COMPS0 (default: /w/20251/sarya/ai-drug-discovery)

set -e
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ON_COMPS0="${REPO_ON_COMPS0:-/w/20251/sarya/ai-drug-discovery}"
SCRIPT="$1"
shift || true
if [[ -z "$SCRIPT" ]]; then
  echo "Usage: run_on_comps0.sh <script> [args...]" >&2
  echo "Example: run_on_comps0.sh run_step_1.sh" >&2
  exit 1
fi
# Script path relative to repo root for the remote side
if [[ "$SCRIPT" == /* ]]; then
  SCRIPT_REL="$SCRIPT"
else
  SCRIPT_REL="$SCRIPT"
fi
# Forward run-specific config so steps use merged config when PIPELINE_RUN_ID is set
ENV_EXPORTS=""
[[ -n "${PIPELINE_RUN_ID:-}" ]] && ENV_EXPORTS="export PIPELINE_RUN_ID='$PIPELINE_RUN_ID'; "
[[ -n "${PIPELINE_STEP_CONFIG:-}" ]] && ENV_EXPORTS="${ENV_EXPORTS}export PIPELINE_STEP_CONFIG='$PIPELINE_STEP_CONFIG'; "
ssh sarya@comps0 "${ENV_EXPORTS}cd $REPO_ON_COMPS0 && srun --partition gpunodes -c 2 --mem=64G -t 60 bash $SCRIPT_REL"
