#!/usr/bin/env python3
"""
step_8.py — LR scoring with stlearn
Reads:  OUT_DIR/step7_<sample>.h5ad      (from step_7.py)
Writes: OUT_DIR/step8_<sample>.h5ad      (one per sample, with LR scores in obsm)
        OUT_DIR/<sample>_lr_scores.csv
        OUT_DIR/<sample>_lr_pvals.csv
        OUT_DIR/step8_<sample>_lr_spatial.pdf
"""

import os, sys, json, warnings
warnings.filterwarnings("ignore")

import numpy  as np
import pandas as pd
import scanpy as sc
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import stlearn as st
    STLEARN_VERSION = st.__version__
except ImportError:
    print("ERROR: stlearn not installed. Run: pip install stlearn==0.4.12", file=sys.stderr)
    sys.exit(1)

CONFIG_PATH = os.environ.get("PIPELINE_STEP_CONFIG", "config.yml")
with open(CONFIG_PATH, encoding="utf-8") as f:
    if CONFIG_PATH.lower().endswith(".json"):
        cfg = json.load(f)
    else:
        cfg = yaml.safe_load(f)

OUT_DIR = cfg["out_dir"]
SAMPLES = cfg["samples"]
DATASET = cfg["dataset_name"]
SPECIES = cfg["species"]
LR      = cfg["lr_scoring"]

sc.settings.verbosity = 1

print(f"====== step_8.py | LR scoring | stlearn {STLEARN_VERSION} ======")
print(f"  DBs: {LR['databases']} | n_pairs={LR['n_pairs']} | min_spots={LR['min_spots']}")

print("  Loading LR pairs via stlearn...")
try:
    lrs = st.tl.cci.load_lrs(LR["databases"], species=SPECIES)
except Exception:
    lrs = st.tl.cci.load_lrs(LR["databases"])
print(f"  Loaded {len(lrs)} LR pairs")


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
    in_path = os.path.join(OUT_DIR, f"step7_{samp}.h5ad")
    if not os.path.exists(in_path):
        print(f"ERROR: {in_path} not found. Run step_7.py first.", file=sys.stderr)
        sys.exit(1)

    print(f"\n  [{samp}] Running LR scoring...")
    adata = sc.read_h5ad(in_path)

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
        fig.suptitle(f"Step 8: Top 5 LR scores - {samp}", fontsize=10)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, f"step8_{samp}_lr_spatial.pdf"),
                    dpi=150, bbox_inches="tight")
        plt.close()
        print(f"    Saved: step8_{samp}_lr_spatial.pdf")

    out_path = os.path.join(OUT_DIR, f"step8_{samp}.h5ad")
    adata.write_h5ad(out_path)
    print(f"    Saved: step8_{samp}.h5ad")

print(f"\n====== step_8.py COMPLETE | {len(SAMPLES)} samples ======")
