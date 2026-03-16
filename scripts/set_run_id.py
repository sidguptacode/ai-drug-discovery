#!/usr/bin/env python3
"""Start a new pipeline run: generate UUID, create runs/<uuid>/ with config_overrides.yml and run_info.json.
Does NOT modify the base config.yml (static). Run-specific params live in runs/<run_id>/config_overrides.yml."""
import datetime
import json
import os
import sys
import uuid
from pathlib import Path

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML required (pip install pyyaml)", file=sys.stderr)
    sys.exit(1)

REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = Path(os.environ.get("PIPELINE_STEP_CONFIG", REPO_ROOT / "config.yml"))

if not CONFIG_PATH.is_file():
    print(f"ERROR: Config not found: {CONFIG_PATH}", file=sys.stderr)
    sys.exit(1)

run_id = str(uuid.uuid4())
base_cfg = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8")) or {}

# Run-specific overrides: run_id and samples (copy from base so run is valid when merged)
overrides = {
    "run_id": run_id,
    "samples": base_cfg.get("samples") or [],
}

run_dir = REPO_ROOT / "runs" / run_id
run_dir.mkdir(parents=True, exist_ok=True)
(run_dir / "logs").mkdir(exist_ok=True)

overrides_path = run_dir / "config_overrides.yml"
with open(overrides_path, "w", encoding="utf-8") as f:
    yaml.dump(overrides, f, default_flow_style=False, allow_unicode=True, sort_keys=False)

run_info = {
    "run_id": run_id,
    "started_at": datetime.datetime.utcnow().isoformat() + "Z",
    "steps_completed": [],
    "notes": "",
}
info_path = run_dir / "run_info.json"
with open(info_path, "w", encoding="utf-8") as f:
    json.dump(run_info, f, indent=2)

print(run_id)
