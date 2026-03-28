#!/usr/bin/env python3
"""
step_9.py — Cell-cell interaction (CCI) analysis with stlearn
Reads:  OUT_DIR/step8_<sample>.h5ad      (from step_8.py)
Writes: OUT_DIR/<sample>_cci_results.csv (one per sample)

Fixes vs original:
  - n_pairs passed to run_cci (was defaulting to 100 instead of config value)
  - BH FDR correction applied across all tests per sample before sig cutoff
  - adj_p_val column added to output
"""

import os, sys, warnings
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

_args = sys.argv[1:]
CONFIG_PATH = _args[0] if _args else \
    "/scratch/baderlab/sgupta/workflows_march/mar9_ai_drug_discovery/config_dipg.yml"
with open(CONFIG_PATH) as f:
    cfg = yaml.safe_load(f)

OUT_DIR    = cfg["out_dir"]
SAMPLES    = cfg["samples"]
LR_CFG     = cfg["lr_scoring"]
CCI_CFG    = cfg["cci"]

N_PAIRS    = LR_CFG["n_pairs"]
SIG_CUTOFF = CCI_CFG["significance_cutoff"]
FDR_CORR   = CCI_CFG.get("fdr_correction", True)

sc.settings.verbosity = 1

print(f"====== step_9.py | CCI analysis | stlearn {STLEARN_VERSION} ======")
print(f"  n_pairs={N_PAIRS} | min_spots={CCI_CFG['min_spots']} | "
      f"significance_cutoff={SIG_CUTOFF} | fdr_correction={FDR_CORR}")

EMPTY_CCI_COLS = ["lr_pair", "sender", "receiver", "cci_score",
                  "p_val", "adj_p_val", "sample"]


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR correction. Returns adjusted p-values."""
    n = len(pvals)
    if n == 0:
        return np.array([])
    pvals = np.asarray(pvals, dtype=float)
    order = np.argsort(pvals)
    adjusted = pvals[order] * n / np.arange(1, n + 1)
    # Enforce monotonicity (cumulative min from largest to smallest rank)
    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])
    result = np.empty(n)
    result[order] = adjusted
    return np.minimum(result, 1.0)


USE_LABEL = "cell_type_label"


def extract_cci_rows(adata, samp):
    """Extract ALL CCI rows without p-value filtering (FDR applied separately).

    Handles both stlearn 1.1.x (keys suffixed with _{use_label}) and
    older stlearn builds (keys named per_lr_cci_pvals / per_lr_results).
    """
    rows = []

    # stlearn 1.1.x: keys are per_lr_cci_pvals_{use_label}
    pvals_key_new = f"per_lr_cci_pvals_{USE_LABEL}"
    score_key_new = f"per_lr_cci_{USE_LABEL}"
    # older stlearn: keys are per_lr_cci_pvals / per_lr_results
    pvals_key_old = "per_lr_cci_pvals"
    score_key_old = "per_lr_results"

    if pvals_key_new in adata.uns:
        pvals_dict = adata.uns[pvals_key_new]
        score_dict = adata.uns.get(score_key_new, {})
    elif pvals_key_old in adata.uns:
        pvals_dict = adata.uns[pvals_key_old]
        score_dict = adata.uns.get(score_key_old, {})
    else:
        pvals_dict = None

    if pvals_dict is not None:
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
                    for x in ["cci", "lr_cci", "cell_cell"])]
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
        try:
            st.tl.cci.run_cci(adata, use_label="cell_type_label",
                              min_spots=CCI_CFG["min_spots"],
                              n_pairs=N_PAIRS,
                              sig_spots=True, verbose=False)
        except TypeError:
            # Older stlearn builds may not accept n_pairs in run_cci
            print(f"    [WARN] run_cci does not accept n_pairs in this stlearn build; "
                  f"permutation count may differ from n_pairs={N_PAIRS}")
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

    if not cci_df.empty and "p_val" in cci_df.columns:
        if FDR_CORR:
            cci_df["adj_p_val"] = bh_fdr(cci_df["p_val"].values)
            n_before = len(cci_df)
            cci_df = cci_df[cci_df["adj_p_val"] <= SIG_CUTOFF].copy()
            print(f"    FDR correction: {n_before} tests -> "
                  f"{len(cci_df)} significant (adj p <= {SIG_CUTOFF})")
        else:
            cci_df["adj_p_val"] = cci_df["p_val"]
            cci_df = cci_df[cci_df["p_val"] <= SIG_CUTOFF].copy()
    else:
        if "adj_p_val" not in cci_df.columns:
            cci_df["adj_p_val"] = np.nan

    cci_df.to_csv(os.path.join(OUT_DIR, f"{samp}_cci_results.csv"), index=False)
    print(f"    {len(cci_df)} significant interactions -> {samp}_cci_results.csv")

print(f"\n====== step_9.py COMPLETE | {len(SAMPLES)} samples ======")
