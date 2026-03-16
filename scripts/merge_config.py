#!/usr/bin/env python3
"""Merge base config YAML with run-specific overrides (deep merge). Overrides win.
Usage: merge_config.py <base.yml> <overrides.yml> [output.yml] [--run-id RUN_ID]
If output.yml is omitted, prints to stdout. If --run-id is set, merged out_dir becomes <repo_root>/outputs/<run_id>."""
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML required (pip install pyyaml)", file=sys.stderr)
    sys.exit(1)


def deep_merge(base: dict, overrides: dict) -> dict:
    """Recursively merge overrides into base. Override values replace base (for non-dicts)."""
    out = dict(base)
    for k, v in overrides.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def main():
    if len(sys.argv) < 3:
        print("Usage: merge_config.py <base.yml> <overrides.yml> [output.yml] [--run-id RUN_ID]", file=sys.stderr)
        sys.exit(1)
    base_path = Path(sys.argv[1]).resolve()
    overrides_path = Path(sys.argv[2])
    out_path = None
    run_id = None
    i = 3
    while i < len(sys.argv):
        if sys.argv[i] == "--run-id" and i + 1 < len(sys.argv):
            run_id = sys.argv[i + 1]
            i += 2
        elif out_path is None and not sys.argv[i].startswith("--"):
            out_path = Path(sys.argv[i])
            i += 1
        else:
            i += 1

    if not base_path.is_file():
        print(f"ERROR: Base config not found: {base_path}", file=sys.stderr)
        sys.exit(1)
    if not overrides_path.is_file():
        print(f"ERROR: Overrides not found: {overrides_path}", file=sys.stderr)
        sys.exit(1)

    base = yaml.safe_load(base_path.read_text(encoding="utf-8")) or {}
    overrides = yaml.safe_load(overrides_path.read_text(encoding="utf-8")) or {}
    merged = deep_merge(base, overrides)

    if run_id:
        repo_root = base_path.parent
        merged["out_dir"] = str((repo_root / "outputs" / run_id).resolve())

    text = yaml.dump(merged, default_flow_style=False, allow_unicode=True, sort_keys=False)
    if out_path is not None:
        out_path.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
