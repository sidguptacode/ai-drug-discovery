#!/usr/bin/env python3
"""Print run_id from config file. Usage: get_run_id.py <config_path>"""
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit(0)
path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("config.yml")
if path.exists():
    c = yaml.safe_load(path.read_text(encoding="utf-8"))
    rid = c.get("run_id") or ""
    if rid:
        print(rid)
