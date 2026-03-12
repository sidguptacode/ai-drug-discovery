#!/usr/bin/env python3
"""
OpenAI Agent Pipeline Controller.
Uses client.responses.create with tools; run loop managed locally (execute tool calls, append
function_call_output, continue until no more function_call items).
"""

import argparse
import json
import os
import subprocess
import sys
import threading
from datetime import datetime
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()

from openai import OpenAI

# Local pipeline modules
import run_pipeline
from pipeline_outputs import read_step_outputs, read_all_step_outputs, step_outputs_exist
from pipeline_metadata import get_pipeline_metadata


# -----------------------------------------------------------------------------
# Tool definitions (OpenAI function-calling format for Responses API)
# -----------------------------------------------------------------------------

PIPELINE_TOOLS = [
    {
        "type": "function",
        "name": "pipeline_init_run",
        "description": "Create run directory and seed step_1.json through step_10.json from config. Run dir can be relative to project root (e.g. runs/run_01).",
        "parameters": {
            "type": "object",
            "properties": {
                "run_dir": {"type": "string", "description": "Path to run directory (e.g. runs/run_01)"},
                "config_path": {"type": "string", "description": "Path to full config file", "default": "config.yml"},
                "custom_config_path": {"type": "string", "description": "Optional override config to merge"},
            },
            "required": ["run_dir"],
        },
    },
    {
        "type": "function",
        "name": "pipeline_get_step_config",
        "description": "Return path and full JSON contents of step_{step}.json for a run.",
        "parameters": {
            "type": "object",
            "properties": {
                "run_dir": {"type": "string", "description": "Path to run directory"},
                "step": {"type": "integer", "description": "Step number 1-10"},
            },
            "required": ["run_dir", "step"],
        },
    },
    {
        "type": "function",
        "name": "pipeline_set_step_config",
        "description": "Update hyperparameters for a step. Only provided keys are changed (merge semantics).",
        "parameters": {
            "type": "object",
            "properties": {
                "run_dir": {"type": "string", "description": "Path to run directory"},
                "step": {"type": "integer", "description": "Step number 1-10"},
                "updates": {"type": "object", "description": "Partial JSON to merge on top of existing step config"},
            },
            "required": ["run_dir", "step", "updates"],
        },
    },
    {
        "type": "function",
        "name": "pipeline_run_step",
        "description": "Run pipeline step N. If force is false, skip if output files already exist and return a message; if force is true, always run. Returns exit code and last N lines of stdout/stderr when executed.",
        "parameters": {
            "type": "object",
            "properties": {
                "run_dir": {"type": "string", "description": "Path to run directory"},
                "step": {"type": "integer", "description": "Step number 1-10"},
                "force": {"type": "boolean", "description": "If false, skip when outputs exist; if true, always run", "default": False},
                "capture_output_lines": {"type": "integer", "description": "Last N lines of stdout/stderr to return", "default": 30},
            },
            "required": ["run_dir", "step"],
        },
    },
    {
        "type": "function",
        "name": "pipeline_read_outputs",
        "description": "Get reflection summary for step(s). When use_stub is true, only verifies that the step's output files exist and returns a short message that the step is complete and need not be reflected on (no file content is inspected). When false, returns detailed summaries from CSV/JSON output files (no RDS/h5ad content).",
        "parameters": {
            "type": "object",
            "properties": {
                "run_dir": {"type": "string", "description": "Path to run directory"},
                "step": {"type": "integer", "description": "Step number 1-10, or omit for all steps"},
                "use_stub": {"type": "boolean", "description": "If true, only check that outputs exist and tell the model the step is complete with no reflection needed; if false, return full summaries from output files", "default": True},
            },
            "required": ["run_dir"],
        },
    },
    {
        "type": "function",
        "name": "pipeline_get_metadata",
        "description": "Get dataset and folder layout and per-step parameter documentation. Use when setting or correcting hyperparameters. If run_dir provided, include current run's out_dir and samples.",
        "parameters": {
            "type": "object",
            "properties": {
                "run_dir": {"type": "string", "description": "Optional path to run directory to include current run info"},
            },
        },
    },
]


# -----------------------------------------------------------------------------
# Agent interaction logging (stdout + logs/<run_id>/<timestamp>.json)
# -----------------------------------------------------------------------------

def _run_id_from_run_dir(run_dir_default: str | Path | None) -> str:
    """Derive a short run_id for log folder naming. Multiple agent runs can share the same run_id."""
    if not run_dir_default:
        return "default"
    name = Path(run_dir_default).name
    return name or "default"


def _unique_log_path(log_dir: Path) -> Path:
    """Return path to a new log file; use timestamp and suffix if collision."""
    log_dir.mkdir(parents=True, exist_ok=True)
    base = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    path = log_dir / f"{base}.json"
    n = 0
    while path.exists():
        n += 1
        path = log_dir / f"{base}_{n}.json"
    return path


class AgentLogger:
    """Logs every agent interaction to stdout and to a JSON log file under logs/<run_id>/."""

    def __init__(self, project_root: Path, run_id: str):
        self.project_root = Path(project_root)
        self.run_id = run_id
        self.log_dir = self.project_root / "logs" / run_id
        self.log_path = _unique_log_path(self.log_dir)
        self.events = []

    def _emit(self, event_type: str, payload: dict) -> None:
        event = {"ts": datetime.now().isoformat(), "type": event_type, **payload}
        self.events.append(event)
        # Print to stdout
        if event_type == "session_start":
            print(f"\n[agent log] session_start run_id={self.run_id} log_file={self.log_path}", flush=True)
            print(f"  message: {payload.get('message', '')[:200]}{'...' if len(payload.get('message', '') or '') > 200 else ''}", flush=True)
        elif event_type == "api_response":
            print(f"\n[agent log] api_response round={payload.get('round')} response_id={payload.get('response_id', '')}", flush=True)
            for item in payload.get("output_summary", []):
                print(f"  output_item: {item}", flush=True)
        elif event_type == "tool_call":
            print(f"\n[agent log] tool_call name={payload.get('name')} call_id={payload.get('call_id', '')}", flush=True)
            print(f"  arguments: {payload.get('arguments', '')[:500]}{'...' if len(payload.get('arguments', '') or '') > 500 else ''}", flush=True)
        elif event_type == "tool_result":
            print(f"[agent log] tool_result call_id={payload.get('call_id', '')} output_len={len(payload.get('output', '') or '')}", flush=True)
            out_preview = (payload.get("output") or "")[:300]
            if len(payload.get("output") or "") > 300:
                out_preview += "..."
            print(f"  output: {out_preview}", flush=True)
        elif event_type == "final_response":
            print(f"\n[agent log] final_response", flush=True)
            print(f"  {payload.get('text', '')[:500]}{'...' if len(payload.get('text', '') or '') > 500 else ''}", flush=True)
        elif event_type == "session_end":
            print(f"[agent log] session_end log_file={self.log_path}", flush=True)

    def session_start(self, message: str, model: str) -> None:
        self._emit("session_start", {"message": message, "model": model, "run_id": self.run_id})

    def api_response(self, round_num: int, response_id: str, output_items: list) -> None:
        output_summary = []
        output_items_serialized = []
        for item in output_items:
            if isinstance(item, dict):
                t = item.get("type", "unknown")
                output_items_serialized.append(dict(item))
                if t == "message":
                    output_summary.append({"type": t, "role": item.get("role"), "content_preview": str(item.get("content", ""))[:200]})
                elif t == "function_call":
                    output_summary.append({"type": t, "name": item.get("name"), "call_id": item.get("call_id")})
                else:
                    output_summary.append({"type": t})
            else:
                try:
                    output_items_serialized.append(item.model_dump() if hasattr(item, "model_dump") else {"_raw": str(item)})
                except Exception:
                    output_items_serialized.append({"_raw": str(item)})
                output_summary.append(str(type(item).__name__))
        self._emit("api_response", {
            "round": round_num, "response_id": response_id,
            "output_summary": output_summary,
            "output_items": output_items_serialized,
        })

    def tool_call(self, name: str, call_id: str, arguments: str) -> None:
        self._emit("tool_call", {"name": name, "call_id": call_id, "arguments": arguments})

    def tool_result(self, call_id: str, output: str) -> None:
        self._emit("tool_result", {"call_id": call_id, "output": output})

    def final_response(self, text: str) -> None:
        self._emit("final_response", {"text": text})

    def session_end(self) -> None:
        self._emit("session_end", {"log_path": str(self.log_path)})
        with open(self.log_path, "w", encoding="utf-8") as f:
            json.dump({"run_id": self.run_id, "events": self.events}, f, indent=2)
        print(f"[agent log] wrote {len(self.events)} events to {self.log_path}", flush=True)


def _resolve_run_dir(run_dir: str | Path, project_root: Path) -> Path:
    """Resolve run_dir to absolute; if relative, interpret relative to project root."""
    p = Path(run_dir)
    if not p.is_absolute():
        p = project_root / p
    return p.resolve()


def _stream_reader(pipe, stream, lines: list) -> None:
    """Read pipe line-by-line; write to stream and append to lines."""
    for line in iter(pipe.readline, ""):
        stream.write(line)
        stream.flush()
        lines.append(line.rstrip("\n"))
    pipe.close()


def _run_step_with_capture(step: int, run_dir: Path, project_root: Path, last_n: int = 30) -> tuple[int, str, str]:
    """Run pipeline step; stream stdout/stderr to terminal with boundaries; return (exit_code, last_n lines of stdout, last_n lines of stderr)."""
    run_dir = Path(run_dir)
    project_root = Path(project_root)
    step_cfg_path = run_dir / f"step_{step}.json"
    if not step_cfg_path.exists():
        return 1, "", f"ERROR: Missing {step_cfg_path}"
    env = {**os.environ, "PIPELINE_STEP_CONFIG": str(step_cfg_path.resolve())}
    steps_dir = project_root / "separated_steps"
    if step <= 5:
        cmd = ["Rscript", str(steps_dir / f"step_{step}.R")]
    else:
        cmd = [sys.executable, str(steps_dir / f"step_{step}.py")]

    print("\n" + "=" * 60, flush=True)
    print(f"  Step {step} subprocess output (stdout/stderr)", flush=True)
    print("=" * 60, flush=True)

    proc = subprocess.Popen(
        cmd,
        env=env,
        cwd=project_root,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )
    out_lines: list[str] = []
    err_lines: list[str] = []
    t_out = threading.Thread(target=_stream_reader, args=(proc.stdout, sys.stdout, out_lines))
    t_err = threading.Thread(target=_stream_reader, args=(proc.stderr, sys.stderr, err_lines))
    t_out.daemon = True
    t_err.daemon = True
    t_out.start()
    t_err.start()
    t_out.join()
    t_err.join()
    returncode = proc.wait()

    print("=" * 60, flush=True)
    print(f"  End step {step} subprocess output", flush=True)
    print("=" * 60 + "\n", flush=True)

    out_tail = "\n".join(out_lines[-last_n:]) if out_lines else ""
    err_tail = "\n".join(err_lines[-last_n:]) if err_lines else ""
    return returncode, out_tail, err_tail


def execute_tool(
    name: str,
    arguments: str,
    project_root: Path,
    run_dir_default: str | Path | None = None,
) -> str:
    """
    Execute one pipeline tool by name with JSON arguments.
    Returns a string (JSON or text) to send as function_call_output.
    """
    project_root = Path(project_root)
    try:
        args = json.loads(arguments) if isinstance(arguments, str) else arguments
    except json.JSONDecodeError as e:
        return json.dumps({"error": f"Invalid JSON arguments: {e}"})

    run_dir = args.get("run_dir", run_dir_default)
    if run_dir is None and name != "pipeline_get_metadata":
        return json.dumps({"error": "run_dir is required"})
    if run_dir is not None:
        run_dir = _resolve_run_dir(run_dir, project_root)

    if name == "pipeline_init_run":
        config_path = args.get("config_path", "config.yml")
        if not Path(config_path).is_absolute():
            config_path = str(project_root / config_path)
        custom = args.get("custom_config_path")
        if custom and not Path(custom).is_absolute():
            custom = str(project_root / custom)
        run_pipeline.setup_run_dir_from_config(config_path, run_dir, custom, project_root)
        return json.dumps({"status": "ok", "run_dir": str(run_dir), "message": "Run directory initialized; step_1.json ... step_10.json written."})

    if name == "pipeline_get_step_config":
        step = int(args["step"])
        path = run_pipeline.get_step_config_path(step, run_dir)
        if not path.exists():
            return json.dumps({"error": f"Missing {path}"})
        with open(path, encoding="utf-8") as f:
            cfg = json.load(f)
        return json.dumps({"path": str(path), "config": cfg})

    if name == "pipeline_set_step_config":
        step = int(args["step"])
        updates = args["updates"]
        path = run_pipeline.get_step_config_path(step, run_dir)
        if not path.exists():
            return json.dumps({"error": f"Missing {path}"})
        with open(path, encoding="utf-8") as f:
            cfg = json.load(f)
        for k, v in updates.items():
            if isinstance(v, dict) and k in cfg and isinstance(cfg[k], dict):
                cfg[k] = {**cfg[k], **v}
            else:
                cfg[k] = v
        with open(path, "w", encoding="utf-8") as f:
            json.dump(cfg, f, indent=2)
        return json.dumps({"status": "ok", "path": str(path), "message": "Step config updated."})

    if name == "pipeline_run_step":
        step = int(args["step"])
        force = args.get("force", False)
        capture_output_lines = int(args.get("capture_output_lines", 30))
        if not force and step_outputs_exist(step, run_dir):
            return json.dumps({
                "status": "skipped",
                "message": "Step already completed; output files present. Use force=true to re-run.",
                "step": step,
            })
        exit_code, out_tail, err_tail = _run_step_with_capture(step, run_dir, project_root, capture_output_lines)
        if exit_code == 0:
            (run_dir / "last_step.txt").write_text(str(step))
        return json.dumps({
            "step": step,
            "exit_code": exit_code,
            "stdout_last_lines": out_tail,
            "stderr_last_lines": err_tail,
            "message": "Step completed successfully." if exit_code == 0 else f"Step failed with exit code {exit_code}.",
        })

    if name == "pipeline_read_outputs":
        step_arg = args.get("step")
        use_stub = args.get("use_stub", True)
        if step_arg is None:
            result = read_all_step_outputs(run_dir, use_stub=use_stub)
        else:
            result = read_step_outputs(int(step_arg), run_dir, use_stub=use_stub)
        return json.dumps(result)

    if name == "pipeline_get_metadata":
        run_dir_arg = args.get("run_dir")
        if run_dir_arg is not None:
            run_dir_arg = _resolve_run_dir(run_dir_arg, project_root)
        meta = get_pipeline_metadata(run_dir_arg)
        return json.dumps(meta)

    return json.dumps({"error": f"Unknown tool: {name}"})


def _default_instructions() -> str:
    return """You are an agent that drives a 10-step spatial transcriptomics pipeline. The goal is to take raw spatial transcriptomics data (multiple samples), run QC and single-cell analysis, annotate cell types, then score ligand–receptor (LR) interactions and produce a ranked list of LR pairs as the final output.

**Overall flow:** Steps 1–5 run in R (Seurat): load samples, QC filter, integrate across samples, cluster, and annotate cell types with marker genes and EnrichR. Step 5 exports cell-type metadata. Steps 6–10 run in Python (scanpy, stlearn): load per-sample data with that metadata, normalize and preprocess, run LR scoring in spatial neighbourhoods, run cell–cell interaction (CCI) testing, then aggregate and rank LR pairs across samples into a single ground-truth table. Steps must be run in order (each step depends on the previous).

**Your tools:** Initialize a run (pipeline_init_run), read or update step configs (pipeline_get_step_config, pipeline_set_step_config), run a step (pipeline_run_step; use force=true to re-run if outputs already exist), read step outputs for reflection (pipeline_read_outputs), and get data/folder metadata and per-step descriptions (pipeline_get_metadata). Call pipeline_get_metadata when you need to set or correct hyperparameters or understand what a step does.

**Paths:** Under data_dir each sample has a folder named exactly <sample_id>. Required files per sample: {sample_id}_filtered_feature_bc_matrix.h5, {sample_id}_tissue_positions_list.csv, {sample_id}_scalefactors_json.json. Outputs go to out_dir (e.g. outputs/<run_id>). Step configs live in run_dir/step_1.json ... step_10.json."""


def run_agent_loop(
    client: OpenAI,
    *,
    model: str = "gpt-4o-mini",
    message: str,
    instructions: str | None = None,
    project_root: Path | None = None,
    run_dir_default: str | None = None,
    max_tool_rounds: int = 50,
) -> str:
    """
    Run the local tool loop: responses.create -> on function_call execute tools and append
    function_call_output -> create again until no more function_calls. Returns final response text.
    """
    if project_root is None:
        project_root = run_pipeline.get_project_root()
    project_root = Path(project_root)
    instructions = instructions or _default_instructions()

    run_id = _run_id_from_run_dir(run_dir_default)
    logger = AgentLogger(project_root, run_id)
    logger.session_start(message, model)

    if run_dir_default:
        user_content = f"The current user task relates to run_id {run_id}.\n\n{message}"
    else:
        user_content = message
    input_list = [{"role": "user", "content": user_content}]
    for round_num in range(1, max_tool_rounds + 1):
        response = client.responses.create(
            model=model,
            tools=PIPELINE_TOOLS,
            instructions=instructions,
            input=input_list,
        )

        output_items = getattr(response, "output", [])
        if hasattr(output_items, "__iter__") and not isinstance(output_items, list):
            output_items = list(output_items)

        response_id = getattr(response, "id", "") or ""
        logger.api_response(round_num, response_id, output_items)

        # Append all response output items to input_list (for context on next request)
        for item in output_items:
            if isinstance(item, dict):
                input_list.append(dict(item))
            elif hasattr(item, "model_dump"):
                input_list.append(item.model_dump())
            else:
                input_list.append({"role": getattr(item, "role", None), "content": getattr(item, "content", item)})

        has_function_call = False
        for item in output_items:
            it = item.model_dump() if hasattr(item, "model_dump") else (item if isinstance(item, dict) else {})
            if it.get("type") != "function_call":
                continue
            has_function_call = True
            name = it.get("name") or (getattr(item, "name", None) if hasattr(item, "name") else None)
            arguments = it.get("arguments") or (getattr(item, "arguments", "{}") if hasattr(item, "arguments") else "{}")
            call_id = it.get("call_id") or it.get("id") or (getattr(item, "call_id", None) or getattr(item, "id", None))
            logger.tool_call(name, call_id, arguments)
            result = execute_tool(name, arguments, project_root, run_dir_default)
            logger.tool_result(call_id, result)
            input_list.append({"type": "function_call_output", "call_id": call_id, "output": result})

        if not has_function_call:
            output_text = getattr(response, "output_text", None)
            if output_text:
                logger.final_response(output_text)
                logger.session_end()
                return output_text
            if output_items:
                for item in output_items:
                    c = getattr(item, "content", None)
                    if c and isinstance(c, list):
                        for block in c:
                            if getattr(block, "type", None) == "output_text":
                                text = getattr(block, "text", "") or ""
                                logger.final_response(text)
                                logger.session_end()
                                return text
                    if isinstance(item, dict) and "content" in item:
                        for block in item.get("content", []):
                            if block.get("type") == "output_text":
                                text = block.get("text", "") or ""
                                logger.final_response(text)
                                logger.session_end()
                                return text
            logger.final_response("(No text output)")
            logger.session_end()
            return "(No text output)"
    logger.final_response("(Max tool rounds reached)")
    logger.session_end()
    return "(Max tool rounds reached)"


def main() -> None:
    parser = argparse.ArgumentParser(description="OpenAI Agent Pipeline Controller (Responses API, local run loop)")
    parser.add_argument("--message", "-m", required=True, help="User message for the agent")
    parser.add_argument("--run-dir", default=None, help="Default run directory for tools (e.g. runs/run_01)")
    parser.add_argument("--model", default="gpt-4.1", help="Model name (default: gpt-4.1)")
    parser.add_argument("--max-tool-rounds", type=int, default=20, help="Max tool-call rounds")
    args = parser.parse_args()

    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        print("ERROR: Set OPENAI_API_KEY", file=sys.stderr)
        sys.exit(1)

    client = OpenAI()
    project_root = run_pipeline.get_project_root()
    final = run_agent_loop(
        client,
        model=args.model,
        message=args.message,
        project_root=project_root,
        run_dir_default=args.run_dir,
        max_tool_rounds=args.max_tool_rounds,
    )
    print("--------------------------------")
    print("Agent loop output text:")
    print(final)
    print("--------------------------------")


if __name__ == "__main__":
    main()
