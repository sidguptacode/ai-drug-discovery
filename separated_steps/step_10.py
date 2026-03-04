#!/usr/bin/env python3
"""
step_10.py — Aggregate and rank LR pairs across samples
Reads:  OUT_DIR/<sample>_cci_results.csv (from step_9.py)
Writes: OUT_DIR/GROUND_TRUTH_lr_pairs_ranked.csv
        OUT_DIR/step10_ranked_lr_pairs.pdf
"""

import os, sys, warnings
warnings.filterwarnings("ignore")

import numpy  as np
import pandas as pd
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

print(f"====== step_10.py | Aggregate & rank LR pairs ======")

EMPTY_CCI_COLS = ["lr_pair","sender","receiver","cci_score","p_val","sample"]

cci_results_all = {}
for samp in SAMPLES:
    cci_path = os.path.join(OUT_DIR, f"{samp}_cci_results.csv")
    if not os.path.exists(cci_path):
        print(f"ERROR: {cci_path} not found. Run step_9.py first.", file=sys.stderr)
        sys.exit(1)
    cci_results_all[samp] = pd.read_csv(cci_path)
    print(f"  [{samp}] {len(cci_results_all[samp])} CCI rows loaded")

all_cci = pd.concat(
    [df for df in cci_results_all.values() if not df.empty],
    ignore_index=True
)

if all_cci.empty:
    print("  [WARN] No CCI data - saving empty ground truth file.")
    pd.DataFrame(columns=["lr_pair","sender","receiver",
                           "n_samples","mean_cci_score","rank"]
                 ).to_csv(os.path.join(OUT_DIR, "GROUND_TRUTH_lr_pairs_ranked.csv"), index=False)
    sys.exit(0)

print(f"\n  Total CCI rows         : {len(all_cci)}")
print(f"  Unique LR pairs        : {all_cci['lr_pair'].nunique()}")
print(f"  Unique cell-type pairs : {all_cci.groupby(['sender','receiver']).ngroups}")
print(f"  Samples with data      : {all_cci['sample'].nunique()}")

ranked = (all_cci
    .groupby(["lr_pair","sender","receiver"])
    .agg(n_samples      = ("sample",    "nunique"),
         mean_cci_score = ("cci_score", "mean"))
    .reset_index()
    .sort_values(["n_samples","mean_cci_score"], ascending=[False,False])
    .reset_index(drop=True))
ranked["rank"] = ranked.index + 1

ranked.to_csv(os.path.join(OUT_DIR, "GROUND_TRUTH_lr_pairs_ranked.csv"), index=False)
print(f"\n  Saved: GROUND_TRUTH_lr_pairs_ranked.csv ({len(ranked)} entries)")

# ── Visualise top 5 senders by interaction breadth ────────────────────────────
top_senders = (all_cci.groupby("sender")["sample"].nunique()
               .nlargest(5).index.tolist())
print(f"  Top 5 senders: {top_senders}")

fig, axes_list = plt.subplots(1, len(top_senders),
                               figsize=(8 * len(top_senders), 8))
if len(top_senders) == 1:
    axes_list = [axes_list]
colors = plt.cm.tab10(np.linspace(0, 0.9, len(top_senders)))

for ax, sender, color in zip(axes_list, top_senders, colors):
    subset = ranked[ranked["sender"] == sender].head(20)
    if subset.empty:
        ax.text(0.5, 0.5, f"No data\n{sender}", ha="center", va="center",
                transform=ax.transAxes); ax.axis("off"); continue
    ax.barh(range(len(subset)), subset["n_samples"],
            color=color, alpha=0.85, edgecolor="white")
    ax.set_yticks(range(len(subset)))
    ax.set_yticklabels([f"{r['lr_pair']}  ->  {r['receiver']}"
                        for _, r in subset.iterrows()], fontsize=7, fontstyle="italic")
    ax.invert_yaxis()
    ax.set_xlabel("Number of samples (recurrence)", fontsize=10)
    ax.set_xlim(0, len(SAMPLES) + 1)
    ax.set_title(f"Sender: {sender}", fontsize=11, fontweight="bold", color=color)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    for j, (_, row) in enumerate(subset.iterrows()):
        ax.text(row["n_samples"] + 0.1, j, str(int(row["n_samples"])),
                va="center", fontsize=7)

fig.suptitle(f"Ranked LR pairs per sender - {DATASET}", fontsize=12, fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "step10_ranked_lr_pairs.pdf"),
            dpi=200, bbox_inches="tight")
plt.close()
print("  Saved: step10_ranked_lr_pairs.pdf")

print(f"\n  Top 20 ranked LR pairs:")
print(f"  {'Rank':<6} {'LR pair':<35} {'Sender':<30} {'Receiver':<30} n_samples")
print("  " + "-"*110)
for _, row in ranked.head(20).iterrows():
    print(f"  {int(row['rank']):<6} {row['lr_pair']:<35} "
          f"{row['sender']:<30} {row['receiver']:<30} {int(row['n_samples'])}")

print(f"\n====== step_10.py COMPLETE - {DATASET} ======")
print("Ground truth: GROUND_TRUTH_lr_pairs_ranked.csv")
