#!/usr/bin/env python3
"""
pipeline.py — ST LR pipeline (dataset-agnostic)
All parameters read from config.yaml

Requirements:
    pip install stlearn==0.4.12 scanpy pandas matplotlib pyyaml
"""

import os
import sys
import json
import warnings
warnings.filterwarnings("ignore")

import numpy  as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy  as sc

try:
    import stlearn as st
    STLEARN_VERSION = st.__version__
except ImportError:
    print("ERROR: stlearn not installed. Run: pip install stlearn==0.4.12", file=sys.stderr)
    sys.exit(1)

sc.settings.verbosity = 1

# ── Load config ────────────────────────────────────────────────────────────────
with open("config.yaml") as f:
    cfg = yaml.safe_load(f)

DATASET     = cfg["dataset_name"]
DATA_DIR    = cfg["data_dir"]
OUT_DIR     = cfg["out_dir"]
SAMPLES     = cfg["samples"]
SPECIES     = cfg["species"]

PRE         = cfg["preprocessing"]
LR          = cfg["lr_scoring"]
CCI_CFG     = cfg["cci"]

os.makedirs(OUT_DIR, exist_ok=True)

print(f"====== pipeline.py | {DATASET} | scanpy {sc.__version__} | stlearn {STLEARN_VERSION} ======")
print(f"  Samples: {len(SAMPLES)} | Species: {SPECIES}")


# =============================================================================
# STEP 6: Load samples + attach cell-type metadata from pipeline.R
# =============================================================================
print("\n================ STEP 6: Load samples ================")

meta_path = os.path.join(OUT_DIR, "cell_type_metadata.csv")
if not os.path.exists(meta_path):
    print(f"ERROR: {meta_path} not found. Run pipeline.R first.", file=sys.stderr)
    sys.exit(1)

meta_df = pd.read_csv(meta_path)
print(f"  Loaded metadata: {len(meta_df)} spots across {meta_df['sample'].nunique()} samples")
print(f"  Cell types: {sorted(meta_df['cell_type_label'].dropna().unique().tolist())}")


def load_sample(sample_name: str, meta_df: pd.DataFrame) -> sc.AnnData:
    h5_path    = os.path.join(DATA_DIR, sample_name,
                               f"{sample_name}_filtered_feature_bc_matrix.h5")
    coord_path = os.path.join(DATA_DIR, sample_name,
                               f"{sample_name}_tissue_positions_list.csv")
    sf_path    = os.path.join(DATA_DIR, sample_name,
                               f"{sample_name}_scalefactors_json.json")

    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()

    try:
        coords = pd.read_csv(coord_path, header=None,
                             names=["barcode","in_tissue","array_row","array_col",
                                    "pxl_row","pxl_col"])
        if coords["in_tissue"].dtype == object:
            raise ValueError("header row detected")
    except Exception:
        coords = pd.read_csv(coord_path)
        coords.columns = [c.lower() for c in coords.columns]

    coords = coords[coords["in_tissue"] == 1].set_index("barcode")
    shared = adata.obs_names.intersection(coords.index)
    adata  = adata[shared].copy()
    adata.obsm["spatial"] = coords.loc[shared, ["pxl_col","pxl_row"]].values.astype(float)

    with open(sf_path) as f:
        adata.uns["spatial"] = {sample_name: {"scalefactors": json.load(f), "images": {}}}

    adata.obs["sample"] = sample_name

    sample_meta = meta_df[meta_df["sample"] == sample_name].set_index("barcode")
    common_bc   = adata.obs_names.intersection(sample_meta.index)
    if len(common_bc) < 10:
        strip   = lambda bc: bc.rsplit("_", 1)[0]
        a_map   = {strip(b): b for b in adata.obs_names}
        m_map   = {strip(b): b for b in sample_meta.index}
        common_stripped = set(a_map) & set(m_map)
        common_bc   = pd.Index([a_map[b] for b in common_stripped])
        meta_idx    = pd.Index([m_map[b] for b in common_stripped])
        sample_meta = sample_meta.loc[meta_idx]
        sample_meta.index = common_bc

    adata = adata[common_bc].copy()
    for col in ["community", "cell_type_label"]:
        if col in sample_meta.columns:
            adata.obs[col] = sample_meta.loc[common_bc, col].values

    print(f"    {sample_name}: {adata.n_obs} spots, {adata.n_vars} genes")
    return adata


adatas_raw = {samp: load_sample(samp, meta_df) for samp in SAMPLES}


# =============================================================================
# STEP 7: Preprocess
# =============================================================================
print("\n================ STEP 7: Preprocess ================")
print(f"  min_genes={PRE['min_genes']} | min_cells={PRE['min_cells']} "
      f"| normalize_target={PRE['normalize_target']}")

fig, axes = plt.subplots(2, len(SAMPLES), figsize=(3 * len(SAMPLES), 6))
adatas = {}

for i, samp in enumerate(SAMPLES):
    adata    = adatas_raw[samp].copy()
    n_before = adata.n_obs
    g_before = adata.n_vars

    counts_before = np.array(adata.X.sum(axis=1)).flatten()
    axes[0, i].hist(counts_before, bins=30, color="steelblue", alpha=0.8, edgecolor="none")
    axes[0, i].set_title(samp, fontsize=6, rotation=20)
    if i == 0: axes[0, i].set_ylabel("spots (before)", fontsize=7)
    axes[0, i].tick_params(labelsize=5)

    sc.pp.filter_cells(adata, min_genes=PRE["min_genes"])
    sc.pp.filter_genes(adata, min_cells=PRE["min_cells"])
    sc.pp.normalize_total(adata, target_sum=PRE["normalize_target"])
    sc.pp.log1p(adata)

    print(f"  {samp}: spots {n_before}→{adata.n_obs}, genes {g_before}→{adata.n_vars}")

    counts_after = np.array(adata.X.sum(axis=1)).flatten()
    axes[1, i].hist(counts_after, bins=30, color="coral", alpha=0.8, edgecolor="none")
    if i == 0: axes[1, i].set_ylabel("spots (after)", fontsize=7)
    axes[1, i].tick_params(labelsize=5)

    adatas[samp] = adata

fig.suptitle(f"Step 7: Count distributions — {DATASET}", fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "step7_count_distributions.pdf"), dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: outputs/step7_count_distributions.pdf")


# =============================================================================
# STEP 8: LR scoring
# =============================================================================
print("\n================ STEP 8: LR scoring ================")
print(f"  DBs: {LR['databases']} | n_pairs={LR['n_pairs']} | "
      f"min_spots={LR['min_spots']} | n_neighbors={LR['n_neighbors']}")

print("  Loading LR pairs via stlearn...")
try:
    lrs = st.tl.cci.load_lrs(LR["databases"], species=SPECIES)
except Exception:
    lrs = st.tl.cci.load_lrs(LR["databases"])
print(f"  Loaded {len(lrs)} LR pairs")

lr_scores_all = {}
lr_pvals_all  = {}


def extract_lr_matrices(adata, lrs):
    for score_key, pval_key in [("lr_scores","lr_pvals"), ("lrscore","lrpval")]:
        if score_key in adata.obsm:
            n      = adata.obsm[score_key].shape[1]
            scores = pd.DataFrame(adata.obsm[score_key],
                                  index=adata.obs_names, columns=lrs[:n])
            pvals  = pd.DataFrame(
                adata.obsm[pval_key] if pval_key in adata.obsm
                else np.ones((adata.n_obs, n)),
                index=adata.obs_names, columns=lrs[:n])
            return scores, pvals
    sk = next((k for k in adata.obsm if "score" in k.lower()), None)
    pk = next((k for k in adata.obsm if "pval"  in k.lower()), None)
    if sk:
        n      = adata.obsm[sk].shape[1]
        scores = pd.DataFrame(adata.obsm[sk], index=adata.obs_names, columns=lrs[:n])
        pvals  = pd.DataFrame(adata.obsm[pk] if pk else np.ones((adata.n_obs, n)),
                              index=adata.obs_names, columns=lrs[:n])
        return scores, pvals
    return pd.DataFrame(index=adata.obs_names), pd.DataFrame(index=adata.obs_names)


for samp in SAMPLES:
    print(f"\n  [{samp}] Running LR scoring...")
    adata = adatas[samp].copy()

    try:
        st.tl.cci.run(adata, lrs,
                      min_spots   = LR["min_spots"],
                      distance    = None,
                      n_pairs     = LR["n_pairs"],
                      n_cpus      = 1,
                      use_label   = "cell_type_label",
                      no_neighbors= LR["n_neighbors"],
                      verbose     = False)
    except TypeError:
        st.tl.cci.run(adata, lrs,
                      min_spots = LR["min_spots"],
                      distance  = None,
                      n_pairs   = LR["n_pairs"],
                      n_cpus    = 1,
                      use_label = "cell_type_label",
                      verbose   = False)

    scores_mat, pvals_mat = extract_lr_matrices(adata, lrs)
    scores_mat.to_csv(os.path.join(OUT_DIR, f"{samp}_lr_scores.csv"))
    pvals_mat.to_csv( os.path.join(OUT_DIR, f"{samp}_lr_pvals.csv"))
    lr_scores_all[samp] = scores_mat
    lr_pvals_all[samp]  = pvals_mat
    adatas[samp] = adata
    print(f"    Saved: {samp}_lr_scores.csv ({scores_mat.shape[1]} pairs)")

    if scores_mat.shape[1] > 0:
        top5 = scores_mat.mean(axis=0).nlargest(5).index.tolist()
        fig, ax_arr = plt.subplots(1, 5, figsize=(20, 4))
        for j, lr_pair in enumerate(top5):
            vals = scores_mat[lr_pair].values
            im   = ax_arr[j].scatter(adata.obsm["spatial"][:,0],
                                     -adata.obsm["spatial"][:,1],
                                     c=vals, cmap="Reds", s=4,
                                     vmin=0, vmax=np.percentile(vals, 95))
            plt.colorbar(im, ax=ax_arr[j], fraction=0.04, pad=0.01)
            ax_arr[j].set_title(lr_pair, fontsize=8)
            ax_arr[j].axis("off")
        fig.suptitle(f"Step 8: Top 5 LR scores — {samp}", fontsize=10)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, f"step8_{samp}_lr_spatial.pdf"),
                    dpi=150, bbox_inches="tight")
        plt.close()


# =============================================================================
# STEP 9: CCI analysis
# =============================================================================
print("\n================ STEP 9: CCI analysis ================")
print(f"  min_spots={CCI_CFG['min_spots']} | significance_cutoff={CCI_CFG['significance_cutoff']}")

cci_results_all = {}
EMPTY_CCI_COLS  = ["lr_pair","sender","receiver","cci_score","p_val","sample"]


def extract_cci_rows(adata, samp):
    rows = []
    sig  = CCI_CFG["significance_cutoff"]

    if "per_lr_cci_pvals" in adata.uns:
        pvals_dict = adata.uns["per_lr_cci_pvals"]
        score_dict = adata.uns.get("per_lr_results", {})
        for lr_pair, pval_df in pvals_dict.items():
            if not isinstance(pval_df, pd.DataFrame):
                try:    pval_df = pd.DataFrame(pval_df)
                except: continue
            score_df = score_dict.get(lr_pair)
            for sender in pval_df.index:
                for receiver in pval_df.columns:
                    pval  = pval_df.loc[sender, receiver]
                    score = (score_df.loc[sender, receiver]
                             if score_df is not None
                             and sender in score_df.index
                             and receiver in score_df.columns
                             else np.nan)
                    if pval <= sig:
                        rows.append({"lr_pair": lr_pair, "sender": sender,
                                     "receiver": receiver, "cci_score": score,
                                     "p_val": pval, "sample": samp})

    elif "ccis" in adata.uns:
        ccis = adata.uns["ccis"]
        if isinstance(ccis, pd.DataFrame):
            ccis = ccis.copy(); ccis["sample"] = samp
            rows = ccis.to_dict("records")
        elif isinstance(ccis, dict):
            for lr_pair, lr_df in ccis.items():
                if isinstance(lr_df, pd.DataFrame):
                    lr_df = lr_df.copy()
                    lr_df["lr_pair"] = lr_pair; lr_df["sample"] = samp
                    rows.extend(lr_df.to_dict("records"))
    else:
        cci_keys = [k for k in adata.uns if any(x in k.lower()
                    for x in ["cci","lr_cci","cell_cell"])]
        print(f"    Available CCI keys: {cci_keys if cci_keys else 'none found'}")
    return rows


for samp in SAMPLES:
    print(f"  [{samp}] Running CCI analysis...")
    adata = adatas[samp]
    try:
        st.tl.cci.run_cci(adata, use_label="cell_type_label",
                           min_spots=CCI_CFG["min_spots"],
                           sig_spots=True, verbose=False)
    except Exception as e:
        print(f"    [WARN] run_cci failed: {e}")
        cci_results_all[samp] = pd.DataFrame(columns=EMPTY_CCI_COLS)
        continue

    rows   = extract_cci_rows(adata, samp)
    cci_df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=EMPTY_CCI_COLS)
    cci_df.to_csv(os.path.join(OUT_DIR, f"{samp}_cci_results.csv"), index=False)
    cci_results_all[samp] = cci_df
    adatas[samp] = adata
    print(f"    {len(cci_df)} significant interactions → {samp}_cci_results.csv")


# =============================================================================
# STEP 10: Aggregate and rank
# =============================================================================
print("\n================ STEP 10: Aggregate & rank LR pairs ================")

all_cci = pd.concat(
    [df for df in cci_results_all.values() if not df.empty],
    ignore_index=True
)

if all_cci.empty:
    print("  [WARN] No CCI data — saving empty ground truth file.")
    pd.DataFrame(columns=["lr_pair","sender","receiver",
                           "n_samples","mean_cci_score","rank"]
                 ).to_csv(os.path.join(OUT_DIR, "GROUND_TRUTH_lr_pairs_ranked.csv"), index=False)
    sys.exit(0)

print(f"  Total CCI rows         : {len(all_cci)}")
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
print(f"\n  Saved: outputs/GROUND_TRUTH_lr_pairs_ranked.csv ({len(ranked)} entries)")

# Visualize top 5 senders by interaction breadth
top_senders = (all_cci.groupby("sender")["sample"].nunique()
               .nlargest(5).index.tolist())
print(f"  Top 5 senders: {top_senders}")

fig, axes_list = plt.subplots(1, len(top_senders),
                               figsize=(8 * len(top_senders), 8))
if len(top_senders) == 1: axes_list = [axes_list]
colors = plt.cm.tab10(np.linspace(0, 0.9, len(top_senders)))

for ax, sender, color in zip(axes_list, top_senders, colors):
    subset = ranked[ranked["sender"] == sender].head(20)
    if subset.empty:
        ax.text(0.5, 0.5, f"No data\n{sender}", ha="center", va="center",
                transform=ax.transAxes); ax.axis("off"); continue
    ax.barh(range(len(subset)), subset["n_samples"],
            color=color, alpha=0.85, edgecolor="white")
    ax.set_yticks(range(len(subset)))
    ax.set_yticklabels([f"{r['lr_pair']}  →  {r['receiver']}"
                        for _, r in subset.iterrows()], fontsize=7, fontstyle="italic")
    ax.invert_yaxis()
    ax.set_xlabel("Number of samples (recurrence)", fontsize=10)
    ax.set_xlim(0, len(SAMPLES) + 1)
    ax.set_title(f"Sender: {sender}", fontsize=11, fontweight="bold", color=color)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    for j, (_, row) in enumerate(subset.iterrows()):
        ax.text(row["n_samples"] + 0.1, j, str(int(row["n_samples"])), va="center", fontsize=7)

fig.suptitle(f"Ranked LR pairs per sender — {DATASET}", fontsize=12, fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "step10_ranked_lr_pairs.pdf"), dpi=200, bbox_inches="tight")
plt.close()
print("  Saved: outputs/step10_ranked_lr_pairs.pdf")

print("\n  ── Top 20 ranked LR pairs ──")
print(f"  {'Rank':<6} {'LR pair':<35} {'Sender':<30} {'Receiver':<30} {'n_samples'}")
print("  " + "─"*110)
for _, row in ranked.head(20).iterrows():
    print(f"  {int(row['rank']):<6} {row['lr_pair']:<35} "
          f"{row['sender']:<30} {row['receiver']:<30} {int(row['n_samples'])}")

print(f"\n====== pipeline.py COMPLETE — {DATASET} ======")
print("Ground truth: outputs/GROUND_TRUTH_lr_pairs_ranked.csv")