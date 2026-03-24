#!/usr/bin/env python3
"""
log_stats.py — Inspect agent run logs produced by AgentLogger.

Usage:
    # Show summary + condensed trace for the most recent log in a run:
    python log_stats.py runs/my_run_01

    # Same, but for a specific log file:
    python log_stats.py logs/my_run_01/2026-03-24_10-00-00.json

    # Show full tool arguments and outputs (verbose):
    python log_stats.py runs/my_run_01 --verbose

    # Print only the summary table (no trace):
    python log_stats.py runs/my_run_01 --summary-only

    # List all log files for a run:
    python log_stats.py runs/my_run_01 --list
"""

import argparse
import json
import os
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_log_files(target: Path, project_root: Path) -> list[Path]:
    """
    Given a target that is either:
      - a .json log file path
      - a run directory (e.g. runs/my_run)
      - a run_id string (e.g. my_run)
    return sorted list of matching log file paths.
    """
    if target.suffix == ".json" and target.exists():
        return [target]

    # Try logs/<run_id>/
    run_id = target.name
    log_dir = project_root / "logs" / run_id
    if not log_dir.exists():
        # Fall back: maybe the user passed the log dir directly
        log_dir = target
    if not log_dir.exists():
        return []
    files = sorted(log_dir.glob("*.json"))
    return files


def _load_log(path: Path) -> dict:
    with open(path, encoding="utf-8") as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def _fmt_tool_call(ev: dict, verbose: bool) -> str:
    name = ev.get("name", "?")
    call_id = ev.get("call_id", "")[:8]
    args_raw = ev.get("arguments", "{}")
    try:
        args = json.loads(args_raw)
    except (json.JSONDecodeError, TypeError):
        args = args_raw

    if verbose:
        args_str = json.dumps(args, indent=4)
    else:
        # Compact single-line, truncated
        args_str = json.dumps(args)
        if len(args_str) > 200:
            args_str = args_str[:200] + "..."
    return f"  TOOL_CALL  [{call_id}] {name}\n    args: {args_str}"


def _fmt_tool_result(ev: dict, verbose: bool) -> str:
    call_id = ev.get("call_id", "")[:8]
    output = ev.get("output", "")
    try:
        parsed = json.loads(output)
        if verbose:
            out_str = json.dumps(parsed, indent=4)
        else:
            # Show status/message/artifacts if present, else truncate
            parts = []
            for key in ("status", "exit_code", "message", "artifacts", "error", "n_clusters", "aggregate_score"):
                if key in parsed:
                    parts.append(f"{key}={json.dumps(parsed[key])}")
            if parts:
                out_str = ", ".join(parts)
            else:
                out_str = json.dumps(parsed)
                if len(out_str) > 300:
                    out_str = out_str[:300] + "..."
    except (json.JSONDecodeError, TypeError):
        out_str = str(output)
        if not verbose and len(out_str) > 300:
            out_str = out_str[:300] + "..."
    return f"  TOOL_RESULT [{call_id}] {out_str}"


def _fmt_api_response(ev: dict) -> str:
    rnd = ev.get("round", "?")
    rid = (ev.get("response_id") or "")[:16]
    items = ev.get("output_summary", [])
    types = []
    for it in items:
        if isinstance(it, dict):
            t = it.get("type", "?")
            if t == "function_call":
                types.append(f"function_call({it.get('name', '?')})")
            elif t == "message":
                types.append("message")
            else:
                types.append(t)
        else:
            types.append(str(it))
    return f"ROUND {rnd}  [{rid}]  output: {', '.join(types) or '(none)'}"


def _fmt_final_response(ev: dict, verbose: bool) -> str:
    text = ev.get("text", "")
    if not verbose and len(text) > 500:
        text = text[:500] + "\n  ... (truncated; use --verbose to see full)"
    return f"\nFINAL RESPONSE:\n{text}"


def _print_trace(log: dict, verbose: bool) -> None:
    events = log.get("events", [])
    print("\n" + "=" * 70)
    print("TRACE")
    print("=" * 70)
    for ev in events:
        ts = ev.get("ts", "")[:19]
        etype = ev.get("type", "")
        if etype == "session_start":
            msg = ev.get("message", "")
            if not verbose and len(msg) > 200:
                msg = msg[:200] + "..."
            print(f"\n[{ts}] SESSION_START  model={ev.get('model','')}  run_id={ev.get('run_id','')}")
            print(f"  message: {msg}")
        elif etype == "api_response":
            print(f"\n[{ts}] {_fmt_api_response(ev)}")
        elif etype == "tool_call":
            print(f"[{ts}] {_fmt_tool_call(ev, verbose)}")
        elif etype == "tool_result":
            print(f"[{ts}] {_fmt_tool_result(ev, verbose)}")
        elif etype == "final_response":
            print(f"[{ts}] {_fmt_final_response(ev, verbose)}")
        elif etype == "session_end":
            print(f"\n[{ts}] SESSION_END  log={ev.get('log_path','')}")


def _print_summary(log: dict, log_path: Path) -> None:
    summary = log.get("summary")
    events = log.get("events", [])

    # Recompute from events if summary missing (older log format)
    if summary is None:
        tool_counts: dict[str, int] = {}
        step_run_counts: dict[str, int] = {}
        total_rounds = 0
        for ev in events:
            if ev["type"] == "tool_call":
                name = ev.get("name", "unknown")
                tool_counts[name] = tool_counts.get(name, 0) + 1
                if name == "pipeline_run_step":
                    try:
                        step_num = json.loads(ev.get("arguments", "{}")).get("step")
                        if step_num is not None:
                            k = str(int(step_num))
                            step_run_counts[k] = step_run_counts.get(k, 0) + 1
                    except (json.JSONDecodeError, TypeError, ValueError):
                        pass
            elif ev["type"] == "api_response":
                total_rounds = max(total_rounds, ev.get("round", 0))
        summary = {
            "total_api_rounds": total_rounds,
            "total_tool_calls": sum(tool_counts.values()),
            "tool_call_counts": tool_counts,
            "step_run_counts": step_run_counts,
            "elapsed_seconds": None,
            "total_input_tokens": None,
            "total_output_tokens": None,
        }

    print("\n" + "=" * 70)
    print(f"SUMMARY  —  {log_path.name}")
    print("=" * 70)
    print(f"  run_id          : {log.get('run_id', '?')}")
    print(f"  api_rounds      : {summary.get('total_api_rounds', '?')}")
    print(f"  total_tool_calls: {summary.get('total_tool_calls', '?')}")
    elapsed = summary.get("elapsed_seconds")
    if elapsed is not None:
        mins, secs = divmod(int(elapsed), 60)
        print(f"  elapsed         : {mins}m {secs}s ({elapsed}s)")
    in_tok = summary.get("total_input_tokens")
    out_tok = summary.get("total_output_tokens")
    if in_tok:
        print(f"  input_tokens    : {in_tok:,}")
    if out_tok:
        print(f"  output_tokens   : {out_tok:,}")

    print("\n  Tool call breakdown:")
    tool_counts = summary.get("tool_call_counts", {})
    if tool_counts:
        max_name = max(len(k) for k in tool_counts)
        for name, cnt in sorted(tool_counts.items(), key=lambda x: -x[1]):
            bar = "█" * cnt
            print(f"    {name:<{max_name}}  {cnt:>3}  {bar}")
    else:
        print("    (none)")

    step_counts = summary.get("step_run_counts", {})
    if step_counts:
        print("\n  pipeline_run_step calls per step:")
        for step, cnt in sorted(step_counts.items(), key=lambda x: int(x[0])):
            bar = "█" * cnt
            print(f"    step {step:>2}  {cnt:>3}  {bar}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Inspect agent run logs (trace + statistics)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "target",
        help="Run directory (e.g. runs/my_run_01), run_id, or direct .json log path",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Show full tool arguments and outputs",
    )
    parser.add_argument(
        "--summary-only", "-s", action="store_true",
        help="Print only the summary table, skip the trace",
    )
    parser.add_argument(
        "--list", "-l", action="store_true",
        help="List all available log files for the run and exit",
    )
    parser.add_argument(
        "--all", "-a", action="store_true",
        help="Process all log files for the run (default: most recent only)",
    )
    parser.add_argument(
        "--project-root", default=None,
        help="Project root directory (default: directory of this script)",
    )
    args = parser.parse_args()

    project_root = Path(args.project_root) if args.project_root else Path(__file__).parent.resolve()
    target = Path(args.target)
    if not target.is_absolute():
        target = (project_root / target).resolve()

    log_files = _find_log_files(target, project_root)

    if not log_files:
        print(f"No log files found for: {args.target}", file=sys.stderr)
        print(f"  Looked in: {project_root / 'logs' / Path(args.target).name}", file=sys.stderr)
        sys.exit(1)

    if args.list:
        print(f"Log files for '{Path(args.target).name}':")
        for p in log_files:
            size_kb = p.stat().st_size / 1024
            print(f"  {p.name}  ({size_kb:.1f} KB)")
        return

    files_to_process = log_files if args.all else [log_files[-1]]

    for log_path in files_to_process:
        log = _load_log(log_path)
        _print_summary(log, log_path)
        if not args.summary_only:
            _print_trace(log, args.verbose)
        print()


if __name__ == "__main__":
    main()
