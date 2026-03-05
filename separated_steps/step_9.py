#!/usr/bin/env python3
"""
step_9.py — Cell-cell interaction (CCI) analysis with stlearn
Reads:  OUT_DIR/step8_<sample>.h5ad      (from step_8.py)
Writes: OUT_DIR/<sample>_cci_results.csv (one per sample)
"""

import os, sys, json, warnings
warnings.filterwarnings("ignore")

import numpy  as np
import pandas as pd
import scanpy as sc
import yaml

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
CCI_CFG = cfg["cci"]

sc.settings.verbosity = 1

print(f"====== step_9.py | CCI analysis | stlearn {STLEARN_VERSION} ======")
print(f"  min_spots={CCI_CFG['min_spots']} | significance_cutoff={CCI_CFG['significance_cutoff']}")

EMPTY_CCI_COLS = ["lr_pair","sender","receiver","cci_score","p_val","sample"]


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
    in_path = os.path.join(OUT_DIR, f"step8_{samp}.h5ad")
    if not os.path.exists(in_path):
        print(f"ERROR: {in_path} not found. Run step_8.py first.", file=sys.stderr)
        sys.exit(1)

    print(f"  [{samp}] Running CCI analysis...")
    adata = sc.read_h5ad(in_path)

    try:
        st.tl.cci.run_cci(adata, use_label="cell_type_label",
                          min_spots=CCI_CFG["min_spots"],
                          sig_spots=True, verbose=False)
    except Exception as e:
        print(f"    [WARN] run_cci failed: {e}")
        pd.DataFrame(columns=EMPTY_CCI_COLS).to_csv(
            os.path.join(OUT_DIR, f"{samp}_cci_results.csv"), index=False)
        continue

    rows   = extract_cci_rows(adata, samp)
    cci_df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=EMPTY_CCI_COLS)
    cci_df.to_csv(os.path.join(OUT_DIR, f"{samp}_cci_results.csv"), index=False)
    print(f"    {len(cci_df)} significant interactions -> {samp}_cci_results.csv")

print(f"\n====== step_9.py COMPLETE | {len(SAMPLES)} samples ======")
