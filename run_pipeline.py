#!/usr/bin/env python3
"""
Pipeline controller: agent-ready "run step" tools + minimal CLI (rerun / next step).

Tools (for agent or programmatic use):
  - get_project_root() -> Path
  - run_step(step, run_dir, project_root=None) -> int   # 0 = success, else exit code; run_dir = runs/<run_id>
  - setup_run_dir_from_config(config_path, run_dir, ..., project_root=None) -> None  # run_dir = runs/<run_id>

CLI: user provides run_id. Configs go to runs/<run_id>/, outputs to outputs/<run_id>/.
Interactive loop at each step: [r]erun | [n]ext. Use --init to create the run and seed JSONs from config.yml once.
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None

# -----------------------------------------------------------------------------
# Per-step config schema (for seeding from full config; agent uses same paths)
# -----------------------------------------------------------------------------
STEP_CONFIG = {
    1: (["data_dir", "out_dir", "samples"], ["qc"]),
    2: (["out_dir", "samples"], ["integration"]),
    3: (["out_dir"], ["clustering"]),
    4: (["out_dir"], ["annotation"]),
    5: (["out_dir", "samples"], []),
    6: (["data_dir", "out_dir", "samples"], []),
    7: (["out_dir", "samples", "dataset_name"], ["preprocessing"]),
    8: (["out_dir", "samples", "dataset_name", "species"], ["lr_scoring"]),
    9: (["out_dir", "samples"], ["cci"]),
    10: (["out_dir", "samples", "dataset_name"], []),
}

FIRST_STEP = 1
LAST_STEP = 10


# -----------------------------------------------------------------------------
# Agent-facing tools (stable API for an agent to call as tools)
# -----------------------------------------------------------------------------

def get_project_root() -> Path:
    """
    Project root (directory containing run_pipeline.py and separated_steps/).
    Agent tool: use this as cwd when invoking run_step.
    """
    return Path(__file__).resolve().parent


def run_step(
    step: int,
    run_dir: Path,
    project_root: Path | None = None,
    *,
    rscript: str = "Rscript",
) -> int:
    """
    Run a single pipeline step (1-10).

    - Reads config from run_dir/step_{step}.json (must exist).
    - Sets env PIPELINE_STEP_CONFIG to that path.
    - Runs R (steps 1-5) or Python (steps 6-10) subprocess; cwd = project_root.

    Returns:
        0 on success, subprocess exit code on failure.

    Agent tool signature:
        run_step(step: int, run_dir: str | Path, project_root: str | Path | None = None) -> int
    """
    if project_root is None:
        project_root = get_project_root()
    run_dir = Path(run_dir)
    project_root = Path(project_root)

    if step < FIRST_STEP or step > LAST_STEP:
        raise ValueError(f"step must be {FIRST_STEP}-{LAST_STEP}, got {step}")

    step_cfg_path = run_dir / f"step_{step}.json"
    if not step_cfg_path.exists():
        print(f"ERROR: Missing {step_cfg_path}", file=sys.stderr)
        return 1

    env = {**os.environ, "PIPELINE_STEP_CONFIG": str(step_cfg_path.resolve())}
    steps_dir = project_root / "separated_steps"

    if step <= 5:
        script_path = steps_dir / f"step_{step}.R"
        cmd = [rscript, str(script_path)]
    else:
        script_path = steps_dir / f"step_{step}.py"
        cmd = [sys.executable, str(script_path)]

    result = subprocess.run(cmd, env=env, cwd=project_root)
    return result.returncode


def get_step_config_path(step: int, run_dir: Path) -> Path:
    """
    Path to the config file for a given step (run_dir/step_{step}.json).
    Agent tool: use this to know where to read/write params for the step.
    """
    return Path(run_dir) / f"step_{step}.json"


# -----------------------------------------------------------------------------
# Setup: seed run_dir with step_1.json ... step_10.json from full config
# -----------------------------------------------------------------------------

def _load_config(config_path: str) -> dict:
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")
    if yaml is None:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    with open(path, encoding="utf-8") as f:
        text = f.read()
    if path.suffix.lower() in (".yml", ".yaml"):
        return yaml.safe_load(text) or {}
    if path.suffix.lower() == ".json":
        return json.loads(text)
    try:
        return yaml.safe_load(text) or {}
    except Exception:
        return json.loads(text)


def _deep_merge(base: dict, override: dict) -> dict:
    out = dict(base)
    for k, v in override.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def _slice_config(full_cfg: dict, step: int) -> dict:
    shared, step_keys = STEP_CONFIG[step]
    out = {}
    for key in shared:
        if key in full_cfg:
            out[key] = full_cfg[key]
    for key in step_keys:
        if key in full_cfg:
            out[key] = full_cfg[key]
    return out


def setup_run_dir_from_config(
    config_path: str,
    run_dir: Path,
    custom_config_path: str | None = None,
    project_root: Path | None = None,
) -> None:
    """
    Create run_dir and write step_1.json ... step_10.json from a full config file.
    Configs go to run_dir (runs/<run_id>); out_dir for each step is set to outputs/<run_id>.
    If custom_config_path is given, merge it over config_path (override wins).

    Agent tool: use when initializing a new run from config.yml (or custom override).
    Pass run_dir = project_root / "runs" / run_id. After this, edit step_*.json before/between run_step calls.
    """
    if project_root is None:
        project_root = get_project_root()
    project_root = Path(project_root)
    full_cfg = _load_config(config_path)
    if custom_config_path and Path(custom_config_path).exists():
        full_cfg = _deep_merge(full_cfg, _load_config(custom_config_path))
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    run_id = run_dir.name
    out_dir = (project_root / "outputs" / run_id).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_dir_str = str(out_dir)
    for s in range(FIRST_STEP, LAST_STEP + 1):
        cfg = _slice_config(full_cfg, s)
        cfg["out_dir"] = out_dir_str
        with open(run_dir / f"step_{s}.json", "w", encoding="utf-8") as f:
            json.dump(cfg, f, indent=2)
    print(f"  Wrote per-step configs to {run_dir}")
    print(f"  Outputs for this run: {out_dir_str}")


# -----------------------------------------------------------------------------
# CLI: interactive loop (rerun / next step)
# -----------------------------------------------------------------------------

def _write_last_step(run_dir: Path, step: int) -> None:
    (Path(run_dir) / "last_step.txt").write_text(str(step))


def _read_last_step(run_dir: Path) -> int | None:
    p = Path(run_dir) / "last_step.txt"
    if not p.exists():
        return None
    try:
        return int(p.read_text().strip())
    except ValueError:
        return None


def _cli_interactive_loop(run_dir: Path, project_root: Path, start_step: int) -> int:
    """Run steps with prompt at each: [r]erun | [n]ext. Returns exit code."""
    current = start_step
    while current <= LAST_STEP:
        print(f"\n====== Running step {current} ======")
        rc = run_step(current, run_dir, project_root)
        if rc != 0:
            print(f"Step {current} failed with exit code {rc}", file=sys.stderr)
            choice = input("  [r]erun step | [n]ext step anyway? [r/n] ").strip().lower() or "r"
            if choice != "n":
                continue
        _write_last_step(run_dir, current)
        if current >= LAST_STEP:
            print("\n====== Pipeline complete (steps 1-10) ======")
            return 0
        while True:
            choice = input("  [r]erun step %d | [n]ext step? [r/n] " % current).strip().lower() or "n"
            if choice == "r":
                break
            if choice == "n":
                current += 1
                break
            print("  Enter r or n.")
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Pipeline controller: run steps 1-10. Interactive loop: rerun or next step (you/agent write step_*.json)."
    )
    parser.add_argument(
        "--run-id",
        metavar="ID",
        required=True,
        help="Run identifier. Configs go to runs/<ID>/, outputs to outputs/<ID>/ (create with --init)",
    )
    parser.add_argument(
        "--init",
        action="store_true",
        help="Create run_dir and seed step_*.json from --config (and optional --custom-config); then run loop",
    )
    parser.add_argument(
        "--config",
        default="config.yml",
        metavar="PATH",
        help="Full config file for --init (default: config.yml)",
    )
    parser.add_argument(
        "--custom-config",
        metavar="PATH",
        help="Override config merged with --config when using --init",
    )
    parser.add_argument(
        "--step",
        type=int,
        metavar="N",
        help="Run only step N and exit (no interactive loop). Agent-friendly single-step run.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Start from last completed step (read runs/<run_id>/last_step.txt); only with interactive loop",
    )
    args = parser.parse_args()

    project_root = get_project_root()
    os.chdir(project_root)

    run_id = args.run_id
    run_dir = project_root / "runs" / run_id

    if args.init:
        config_path = Path(args.config)
        if not config_path.is_absolute():
            config_path = project_root / config_path
        custom = None
        if args.custom_config:
            custom = Path(args.custom_config)
            if not custom.is_absolute():
                custom = project_root / custom
            custom = str(custom)
        if not config_path.exists():
            print(f"ERROR: Config not found: {config_path}", file=sys.stderr)
            sys.exit(1)
        setup_run_dir_from_config(str(config_path), run_dir, custom, project_root)

    if not run_dir.is_dir():
        print(f"ERROR: Run not found: {run_dir}. Use --init to create from config.", file=sys.stderr)
        sys.exit(1)

    # Single-step mode (agent can call: run_pipeline.py --run-id <id> --step 3)
    if args.step is not None:
        if args.step < FIRST_STEP or args.step > LAST_STEP:
            print(f"ERROR: --step must be {FIRST_STEP}-{LAST_STEP}", file=sys.stderr)
            sys.exit(1)
        rc = run_step(args.step, run_dir, project_root)
        sys.exit(rc)

    # Interactive loop
    start_step = FIRST_STEP
    if args.resume:
        last = _read_last_step(run_dir)
        if last is not None:
            start_step = last + 1
            print(f"  Resuming from step {start_step} (last completed: {last})")
    rc = _cli_interactive_loop(run_dir, project_root, start_step)
    sys.exit(rc)


if __name__ == "__main__":
    main()
