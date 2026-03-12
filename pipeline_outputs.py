#!/usr/bin/env python3
"""
Reflection layer for the pipeline. use_stub=True: only verify step output files exist,
then tell the model the step is complete and need not be reflected on (no file content
inspected). use_stub=False: full summaries from CSV/JSON only; RDS/h5ad: existence and
non-empty only—never send content to the LLM.
"""

import json
from pathlib import Path


def _load_step_config(run_dir: Path, step: int) -> dict:
    """Load step config JSON; returns dict with out_dir, samples (if present), etc."""
    path = run_dir / f"step_{step}.json"
    if not path.exists():
        return {}
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _out_dir(run_dir: Path, step: int) -> Path | None:
    cfg = _load_step_config(run_dir, step)
    out = cfg.get("out_dir")
    return Path(out) if out else None


def _samples(run_dir: Path, step: int) -> list[str]:
    cfg = _load_step_config(run_dir, step)
    s = cfg.get("samples")
    if s is None:
        return []
    return list(s) if isinstance(s, list) else [s]


def _file_exists_and_nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def read_step_outputs(
    step: int,
    run_dir: Path,
    *,
    use_stub: bool = True,
) -> dict:
    """
    Return reflection summary for a pipeline step.

    When use_stub=True (default): do not inspect any file content. Only check that
    the step's key output files exist; if they do, return a short message that the
    step is complete and need not be reflected on. If outputs are missing, return
    status indicating the step is not complete.
    When use_stub=False: summarize from CSV/JSON only; for RDS/h5ad only check
    existence and non-empty (do not send content to LLM).
    """
    run_dir = Path(run_dir)
    if step < 1 or step > 10:
        return {"step": step, "status": "error", "message": "step must be 1-10"}

    if use_stub:
        if step_outputs_exist(step, run_dir):
            return {
                "step": step,
                "status": "completed",
                "message": "Step completed. Outputs exist; no need to reflect on this step.",
                "artifacts": [],
            }
        return {
            "step": step,
            "status": "incomplete",
            "message": "Step outputs are missing or incomplete. Run the step before reflecting.",
            "artifacts": [],
        }

    out_dir = _out_dir(run_dir, step)
    if not out_dir or not out_dir.is_dir():
        return {
            "step": step,
            "status": "unknown",
            "message": "out_dir not found or not a directory",
            "artifacts": [],
        }

    # Full implementation per step (CSV/JSON only; RDS/h5ad existence only)
    if step == 1:
        rds = out_dir / "step1_seurat_list.rds"
        pdf = out_dir / "step1_qc.pdf"
        artifacts = []
        if _file_exists_and_nonempty(rds):
            artifacts.append("step1_seurat_list.rds")
        if _file_exists_and_nonempty(pdf):
            artifacts.append("step1_qc.pdf")
        return {"step": 1, "status": "completed", "artifacts": artifacts, "message": "See step1_qc.pdf for QC plots."}

    if step == 2:
        rds = out_dir / "step2_seurat_integrated.rds"
        pdf = out_dir / "step2_integration.pdf"
        artifacts = [n for n in ["step2_seurat_integrated.rds", "step2_integration.pdf"] if _file_exists_and_nonempty(out_dir / n)]
        return {"step": 2, "status": "completed", "artifacts": artifacts, "message": "See step2_integration.pdf."}

    if step == 3:
        artifacts = [n for n in ["step3_seurat_clustered.rds", "step3_clustree.pdf", "step3_clusters.pdf"] if _file_exists_and_nonempty(out_dir / n)]
        return {"step": 3, "status": "completed", "artifacts": artifacts, "message": "See step3_clusters.pdf."}

    if step == 4:
        markers_path = out_dir / "step4_markers.csv"
        scores_path = out_dir / "step4_annotation_scores.csv"
        out = {"step": 4, "status": "completed", "artifacts": [], "n_markers": 0, "clusters": [], "annotation_scores": []}


        # Slight bloat here in the relfection later, not necessary to send these markers to LLM
        # if markers_path.exists() and markers_path.stat().st_size > 0:
        #     import csv
        #     with open(markers_path, encoding="utf-8") as f:
        #         rows = list(csv.DictReader(f))
        #     out["n_markers"] = len(rows)
        #     cluster_col = "cluster" if rows and "cluster" in rows[0] else None
        #     if cluster_col:
        #         clusters = sorted({r[cluster_col] for r in rows})
        #         out["clusters"] = clusters
        #         top3 = {}
        #         for r in rows:
        #             c = r[cluster_col]
        #             if c not in top3:
        #                 top3[c] = []
        #             if len(top3[c]) < 3:
        #                 top3[c].append(r.get("gene", r.get("Gene", "")))
        #         out["top_markers_per_cluster"] = top3
        #     out["artifacts"].append("step4_markers.csv")

        if scores_path.exists() and scores_path.stat().st_size > 0:
            import csv
            with open(scores_path, encoding="utf-8") as f:
                out["annotation_scores"] = list(csv.DictReader(f))
            out["artifacts"].append("step4_annotation_scores.csv")
        out["message"] = "See step4_annotation.pdf for plots."
        return out

    if step == 5:
        meta_path = out_dir / "cell_type_metadata.csv"
        out = {"step": 5, "status": "completed", "artifacts": [], "n_cells": 0, "cell_types": []}
        if meta_path.exists() and meta_path.stat().st_size > 0:
            import csv
            with open(meta_path, encoding="utf-8") as f:
                rows = list(csv.DictReader(f))
            out["n_cells"] = len(rows)
            ct_col = "cell_type" if rows and "cell_type" in rows[0] else "cell_type_label"
            if rows and ct_col in rows[0]:
                out["cell_types"] = sorted({r[ct_col] for r in rows})
            out["artifacts"].append("cell_type_metadata.csv")
        return out

    if step == 6:
        samples = _samples(run_dir, 6)
        artifacts = [f"step6_{s}.h5ad" for s in samples if _file_exists_and_nonempty(out_dir / f"step6_{s}.h5ad")]
        return {"step": 6, "status": "completed", "artifacts": artifacts, "message": f"{len(artifacts)} sample h5ad files (content not sent)."}

    if step == 7:
        samples = _samples(run_dir, 7)
        artifacts = [f"step7_{s}.h5ad" for s in samples if _file_exists_and_nonempty(out_dir / f"step7_{s}.h5ad")]
        if _file_exists_and_nonempty(out_dir / "step7_count_distributions.pdf"):
            artifacts.append("step7_count_distributions.pdf")
        return {"step": 7, "status": "completed", "artifacts": artifacts, "message": "See step7_count_distributions.pdf."}

    if step == 8:
        samples = _samples(run_dir, 8)
        out = {"step": 8, "status": "completed", "samples": {}, "artifacts": []}
        for s in samples:
            csv_path = out_dir / f"{s}_lr_scores.csv"
            if csv_path.exists() and csv_path.stat().st_size > 0:
                import csv
                import statistics
                with open(csv_path, encoding="utf-8") as f:
                    rows = list(csv.DictReader(f))
                n_spots = len(rows)
                cols = [k for k in (rows[0].keys() if rows else []) if k and not str(k).startswith("Unnamed")]
                n_pairs = len(cols)
                means = []
                for c in cols:
                    vals = []
                    for row in rows:
                        v = row.get(c)
                        try:
                            vals.append(float(v))
                        except (TypeError, ValueError):
                            pass
                    means.append((c, statistics.mean(vals)) if vals else (c, 0.0))
                means.sort(key=lambda x: -x[1])
                top5 = [m[0] for m in means[:5]]
                out["samples"][s] = {"n_spots": n_spots, "n_pairs": n_pairs, "top5_lr_pairs": top5}
                out["artifacts"].append(f"{s}_lr_scores.csv")
        return out

    if step == 9:
        samples = _samples(run_dir, 9)
        out = {"step": 9, "status": "completed", "samples": {}, "artifacts": []}
        for s in samples:
            csv_path = out_dir / f"{s}_cci_results.csv"
            if csv_path.exists() and csv_path.stat().st_size > 0:
                import csv
                with open(csv_path, encoding="utf-8") as f:
                    rows = list(csv.DictReader(f))
                n = len(rows)
                top5 = rows[:5] if rows else []
                out["samples"][s] = {"n_cci_rows": n, "top5": [{"lr_pair": r.get("lr_pair"), "sender": r.get("sender"), "receiver": r.get("receiver"), "cci_score": r.get("cci_score")} for r in top5]}
                out["artifacts"].append(f"{s}_cci_results.csv")
        return out

    if step == 10:
        csv_path = out_dir / "GROUND_TRUTH_lr_pairs_ranked.csv"
        out = {"step": 10, "status": "completed", "artifacts": [], "total_ranked": 0, "top20": []}
        if csv_path.exists() and csv_path.stat().st_size > 0:
            import csv
            with open(csv_path, encoding="utf-8") as f:
                rows = list(csv.DictReader(f))
            out["total_ranked"] = len(rows)
            out["top20"] = rows[:20]
            out["artifacts"].append("GROUND_TRUTH_lr_pairs_ranked.csv")
        out["message"] = "See step10_ranked_lr_pairs.pdf."
        return out

    return {"step": step, "status": "completed", "artifacts": [], "message": "No summary defined."}


def read_all_step_outputs(run_dir: Path, *, use_stub: bool = True) -> dict:
    """Return reflection summary for all steps 1-10. Keyed by step number."""
    run_dir = Path(run_dir)
    return {s: read_step_outputs(s, run_dir, use_stub=use_stub) for s in range(1, 11)}


# Key output files per step for force/skip check (all must exist and be non-empty)
_STEP_KEY_OUTPUTS = {
    1: ["step1_seurat_list.rds", "step1_qc.pdf"],
    2: ["step2_seurat_integrated.rds"],
    3: ["step3_seurat_clustered.rds"],
    4: ["step4_seurat_annotated.rds", "step4_markers.csv", "step4_annotation_scores.csv"],
    5: ["cell_type_metadata.csv"],
    # 6, 7, 8, 9: depend on samples
    10: ["GROUND_TRUTH_lr_pairs_ranked.csv"],
}


def step_outputs_exist(step: int, run_dir: Path) -> bool:
    """
    Return True if the step's key output files exist and are non-empty.
    Used when force=False to skip re-running the step.
    """
    run_dir = Path(run_dir)
    if step < 1 or step > 10:
        return False
    out_dir = _out_dir(run_dir, step)
    if not out_dir or not out_dir.is_dir():
        return False
    samples = _samples(run_dir, step)

    if step == 6:
        if not samples:
            return False
        return all(_file_exists_and_nonempty(out_dir / f"step6_{s}.h5ad") for s in samples)
    if step == 7:
        pdf = _file_exists_and_nonempty(out_dir / "step7_count_distributions.pdf")
        if not samples:
            return pdf
        return pdf and any(_file_exists_and_nonempty(out_dir / f"step7_{s}.h5ad") for s in samples)
    if step == 8:
        if not samples:
            return False
        return all(_file_exists_and_nonempty(out_dir / f"{s}_lr_scores.csv") for s in samples)
    if step == 9:
        if not samples:
            return False
        return all(_file_exists_and_nonempty(out_dir / f"{s}_cci_results.csv") for s in samples)
    keys = _STEP_KEY_OUTPUTS.get(step, [])
    return all(_file_exists_and_nonempty(out_dir / k) for k in keys)
