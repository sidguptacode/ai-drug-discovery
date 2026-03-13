#!/usr/bin/env python3
"""Generate a run UUID, set it in config.yml, create runs/<uuid>/ and run_info.json."""
import datetime
import json
import os
import re
import sys
import uuid
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = os.environ.get("PIPELINE_STEP_CONFIG", REPO_ROOT / "config.yml")
CONFIG_PATH = Path(CONFIG_PATH)

if not CONFIG_PATH.is_file():
    print(f"ERROR: Config not found: {CONFIG_PATH}", file=sys.stderr)
    sys.exit(1)

run_id = str(uuid.uuid4())
# Update only run_id line so we don't alter comments/formatting
text = CONFIG_PATH.read_text(encoding="utf-8")
text = re.subn(r'^run_id:\s*(?:null|"[^"]*"|[^\s"#]+)?\s*(?:#.*)?$', f'run_id: "{run_id}"', text, count=1, flags=re.MULTILINE)[0]
if not re.search(r'^run_id:\s*', text, re.MULTILINE):
    text += f'\nrun_id: "{run_id}"\n'
CONFIG_PATH.write_text(text, encoding="utf-8")

run_dir = REPO_ROOT / "runs" / run_id
run_dir.mkdir(parents=True, exist_ok=True)
(run_dir / "logs").mkdir(exist_ok=True)
run_info = {
    "run_id": run_id,
    "started_at": None,
    "steps_completed": [],
    "notes": "",
}
# Allow agent to add started_at, steps_completed, notes later
info_path = run_dir / "run_info.json"
run_info["started_at"] = datetime.datetime.utcnow().isoformat() + "Z"
with open(info_path, "w", encoding="utf-8") as f:
    json.dump(run_info, f, indent=2)

print(run_id)
