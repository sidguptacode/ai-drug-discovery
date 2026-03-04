#!/usr/bin/env python3
"""
step_7.py — Preprocess (normalize, log-transform, filter)
Reads:  OUT_DIR/step6_<sample>.h5ad      (from step_6.py)
Writes: OUT_DIR/step7_<sample>.h5ad      (one per sample)
        OUT_DIR/step7_count_distributions.pdf
"""

import os, sys, warnings
warnings.filterwarnings("ignore")

import numpy  as np
import pandas as pd
import scanpy as sc
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CONFIG_PATH = "/scratch/baderlab/sgupta/ai-drug-discovery/config.yml"
with open(CONFIG_PATH) as f:
    cfg = yaml.safe_load(f)

OUT_DIR = cfg["out_dir"]
SAMPLES = cfg["samples"]
DATASET = cfg["dataset_name"]
PRE     = cfg["preprocessing"]

sc.settings.verbosity = 1

print(f"====== step_7.py | Preprocess | scanpy {sc.__version__} ======")
print(f"  min_genes={PRE['min_genes']} | min_cells={PRE['min_cells']} "
      f"| normalize_target={PRE['normalize_target']}")

fig, axes = plt.subplots(2, len(SAMPLES), figsize=(3 * len(SAMPLES), 6))

for i, samp in enumerate(SAMPLES):
    in_path = os.path.join(OUT_DIR, f"step6_{samp}.h5ad")
    if not os.path.exists(in_path):
        print(f"ERROR: {in_path} not found. Run step_6.py first.", file=sys.stderr)
        sys.exit(1)

    adata    = sc.read_h5ad(in_path)
    n_before = adata.n_obs
    g_before = adata.n_vars

    counts_before = np.array(adata.X.sum(axis=1)).flatten()
    axes[0, i].hist(counts_before, bins=30, color="steelblue", alpha=0.8, edgecolor="none")
    axes[0, i].set_title(samp, fontsize=6, rotation=20)
    if i == 0:
        axes[0, i].set_ylabel("spots (before)", fontsize=7)
    axes[0, i].tick_params(labelsize=5)

    sc.pp.filter_cells(adata, min_genes=PRE["min_genes"])
    sc.pp.filter_genes(adata, min_cells=PRE["min_cells"])
    sc.pp.normalize_total(adata, target_sum=PRE["normalize_target"])
    sc.pp.log1p(adata)

    print(f"  {samp}: spots {n_before}->{adata.n_obs}, genes {g_before}->{adata.n_vars}")

    counts_after = np.array(adata.X.sum(axis=1)).flatten()
    axes[1, i].hist(counts_after, bins=30, color="coral", alpha=0.8, edgecolor="none")
    if i == 0:
        axes[1, i].set_ylabel("spots (after)", fontsize=7)
    axes[1, i].tick_params(labelsize=5)

    out_path = os.path.join(OUT_DIR, f"step7_{samp}.h5ad")
    adata.write_h5ad(out_path)
    print(f"    Saved: step7_{samp}.h5ad")

fig.suptitle(f"Step 7: Count distributions - {DATASET}", fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "step7_count_distributions.pdf"),
            dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: step7_count_distributions.pdf")

print(f"\n====== step_7.py COMPLETE | {len(SAMPLES)} samples ======")
