#!/usr/bin/env python3
"""
Pipeline metadata for the agent: data folder layout, output layout, per-step parameters.
"""

import json
from pathlib import Path


# Output files per step (for metadata description only)
OUTPUT_LAYOUT = {
    1: ["step1_seurat_list.rds", "step1_qc.pdf"],
    2: ["step2_seurat_integrated.rds", "step2_integration.pdf"],
    3: ["step3_seurat_clustered.rds", "step3_clustree.pdf", "step3_clusters.pdf"],
    4: ["step4_seurat_annotated.rds", "step4_markers.csv", "step4_enrichr_<community>_<db>.csv", "step4_annotation.pdf", "step4_annotation_scores.csv"],
    5: ["cell_type_metadata.csv", "seurat_integrated.rds"],
    6: ["step6_<sample>.h5ad"],
    7: ["step7_<sample>.h5ad", "step7_count_distributions.pdf"],
    8: ["step8_<sample>.h5ad", "<sample>_lr_scores.csv", "<sample>_lr_pvals.csv", "step8_<sample>_lr_spatial.pdf"],
    9: ["<sample>_cci_results.csv"],
    10: ["GROUND_TRUTH_lr_pairs_ranked.csv", "step10_ranked_lr_pairs.pdf"],
}

# Per-step: short description (1-2 lines) + config keys and hyperparameters
STEP_PARAMS = {
    1: {
        "description": "QC: Load each sample from data_dir, apply nCount and optional mitochondrial filters. Produces per-sample Seurat objects and a QC plot.",
        "shared": ["data_dir", "out_dir", "samples"],
        "hyperparameters": {
            "qc.ncount_min": "integer, minimum total UMI counts per spot",
            "qc.ncount_max": "integer, maximum total UMI counts per spot",
            "qc.mt_cutoff": "number or null, percent mitochondrial filter; null = no filter",
        },
    },
    2: {
        "description": "Integration: SCTransform per sample, then RPCA anchor integration across samples. Produces one integrated Seurat object and UMAP.",
        "shared": ["out_dir", "samples"],
        "hyperparameters": {
            "integration.n_features": "integer, variable features for SCTransform",
            "integration.n_pcs": "integer, PCs computed",
            "integration.n_dims": "integer, PCs for neighbour graph, UMAP",
        },
    },
    3: {
        "description": "Clustering: Build neighbour graph and run Louvain at given resolution(s). Choose one resolution for downstream; produces cluster labels and plots.",
        "shared": ["out_dir"],
        "hyperparameters": {
            "clustering.resolutions": "list of floats, resolution values to try",
            "clustering.chosen_resolution": "float, resolution used for downstream",
        },
    },
    4: {
        "description": "Annotation: Find marker genes per cluster, run EnrichR for cell-type labels, assign labels. Produces annotated object and annotation scores.",
        "shared": ["out_dir"],
        "hyperparameters": {
            "annotation.n_marker_genes": "integer, top N marker genes for EnrichR",
            "annotation.min_pct": "float, min fraction of cells expressing gene",
            "annotation.logfc_threshold": "float, minimum log2 fold-change",
            "annotation.enrichr_dbs": "list of strings, EnrichR database names",
            "annotation.primary_db": "string, primary annotation database",
            "annotation.label_disqualify_patterns": "optional list of strings (regex); EnrichR terms matching any are excluded. When provided, supersedes the default list for this run.",
            "annotation.label_prefer_patterns": "optional list of strings (regex); terms matching any are preferred; best by p-value, else best overall with fallback prefix. When provided, supersedes the default list for this run.",
            "annotation.label_fallback_prefix": "optional string; prefix when chosen label did not match prefer patterns. When provided, supersedes the default for this run.",
        },
    },
    5: {
        "description": "Export: Write cell_type_metadata.csv (barcode, cluster, cell_type_label) and RDS for the Python pipeline. Step 6 needs this CSV.",
        "shared": ["out_dir", "samples"],
        "hyperparameters": {},
    },
    6: {
        "description": "Load Python: Read each sample’s counts and coords from data_dir, attach cell-type labels from step 5. Produces one h5ad per sample.",
        "shared": ["data_dir", "out_dir", "samples"],
        "hyperparameters": {},
    },
    7: {
        "description": "Preprocessing: Normalize, log-transform, filter by min_genes/min_cells. Produces normalized h5ad per sample and count distribution plot.",
        "shared": ["out_dir", "samples", "dataset_name"],
        "hyperparameters": {
            "preprocessing.min_genes": "integer, minimum genes per spot",
            "preprocessing.min_cells": "integer, minimum cells expressing a gene",
            "preprocessing.normalize_target": "integer, library size normalization target",
        },
    },
    8: {
        "description": "LR scoring: Load LR database (e.g. CellPhoneDB), score ligand–receptor pairs in spatial neighbourhoods per sample. Produces LR scores and p-values per spot.",
        "shared": ["out_dir", "samples", "dataset_name", "species"],
        "hyperparameters": {
            "lr_scoring.databases": "list of strings, LR databases (e.g. cellphonedb)",
            "lr_scoring.n_pairs": "integer, permutation pairs for significance",
            "lr_scoring.min_spots": "integer, minimum spots for LR pair expression",
            "lr_scoring.n_neighbors": "integer, spatial neighbourhood size",
        },
    },
    9: {
        "description": "CCI: Test cell–cell interactions per sender–receiver and LR pair; keep significant ones. Produces one CCI results CSV per sample.",
        "shared": ["out_dir", "samples"],
        "hyperparameters": {
            "cci.min_spots": "integer, minimum spots for CCI testing",
            "cci.significance_cutoff": "float, p-value threshold (e.g. 0.05)",
        },
    },
    10: {
        "description": "Aggregate and rank: Combine CCI results across samples, rank LR pairs by recurrence and mean score. Final output is GROUND_TRUTH_lr_pairs_ranked.csv.",
        "shared": ["out_dir", "samples", "dataset_name"],
        "hyperparameters": {},
    },
}


def _load_step_config(run_dir: Path, step: int) -> dict:
    path = run_dir / f"step_{step}.json"
    if not path.exists():
        return {}
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _read_last_step(run_dir: Path) -> int | None:
    """Read last completed step from run_dir/last_step.txt (same as run_pipeline controller)."""
    p = Path(run_dir) / "last_step.txt"
    if not p.exists():
        return None
    try:
        return int(p.read_text().strip())
    except (ValueError, OSError):
        return None


def _run_dir_established(run_dir: Path) -> bool:
    """True if run_dir exists and has step_1.json through step_10.json (initialized)."""
    if not run_dir.is_dir():
        return False
    for step in range(1, 11):
        if not (run_dir / f"step_{step}.json").exists():
            return False
    return True


def get_pipeline_metadata(run_dir: Path | str | None = None) -> dict:
    """
    Return dataset and folder layout and per-step parameter documentation.
    If run_dir is provided, include: which run dir it is, whether it is established,
    last completed step (from last_step.txt), and current out_dir/samples when established.
    """
    out = {
        "data_folder_layout": {
            "description": "data_dir is the root; under it each sample has its own folder named exactly <sample_id>.",
            "required_files_per_sample": [
                "{sample_id}_filtered_feature_bc_matrix.h5",
                "{sample_id}_tissue_positions_list.csv",
                "{sample_id}_scalefactors_json.json",
            ],
            "note": "Sample names in config 'samples' must match folder names under data_dir.",
        },
        "output_folder_layout": {
            "description": "out_dir is set per run (e.g. outputs/<run_id>). Steps write the following files:",
            "by_step": OUTPUT_LAYOUT,
        },
        "per_step_parameters": STEP_PARAMS,
        "step_config_location": "Each step reads config from run_dir/step_{step}.json. Edit these JSON files to change hyperparameters.",
    }
    if run_dir is not None:
        run_dir = Path(run_dir).resolve()
        run_id = run_dir.name or "default"
        established = _run_dir_established(run_dir)
        last_completed_step = _read_last_step(run_dir) if run_dir.is_dir() else None

        current_run = {
            "run_dir": str(run_dir),
            "run_id": run_id,
            "established": established,
            "last_completed_step": last_completed_step,
            "note": "Run dir is established when it contains step_1.json ... step_10.json (e.g. after pipeline_init_run). last_completed_step is the highest step number already completed (from last_step.txt); next step to run is last_completed_step + 1.",
        }
        if last_completed_step is not None and last_completed_step < 10:
            current_run["next_step_to_run"] = last_completed_step + 1
        elif last_completed_step == 10:
            current_run["next_step_to_run"] = None
            current_run["pipeline_complete"] = True
        elif not established:
            current_run["next_step_to_run"] = 1
            current_run["note"] = "Run dir not yet established. Call pipeline_init_run first to create step_1.json ... step_10.json, then run steps from 1."

        if established:
            cfg = _load_step_config(run_dir, 1)
            current_run["out_dir"] = cfg.get("out_dir")
            current_run["samples"] = cfg.get("samples", [])
        out["current_run"] = current_run
    return out
