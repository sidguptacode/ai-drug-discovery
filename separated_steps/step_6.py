#!/usr/bin/env python3
"""
step_6.py — Load samples and attach cell-type metadata
Reads:  OUT_DIR/cell_type_metadata.csv   (from step_5.R)
        DATA_DIR/<sample>/<sample>_filtered_feature_bc_matrix.h5
        DATA_DIR/<sample>/<sample>_tissue_positions_list.csv
        DATA_DIR/<sample>/<sample>_scalefactors_json.json
Writes: OUT_DIR/step6_<sample>.h5ad      (one per sample)
"""

import os, sys, json, warnings
warnings.filterwarnings("ignore")

import numpy  as np
import pandas as pd
import scanpy as sc
import yaml

CONFIG_PATH = os.environ.get("PIPELINE_STEP_CONFIG", "config.yml")
with open(CONFIG_PATH, encoding="utf-8") as f:
    if CONFIG_PATH.lower().endswith(".json"):
        cfg = json.load(f)
    else:
        cfg = yaml.safe_load(f)

DATA_DIR = cfg["data_dir"]
OUT_DIR  = cfg["out_dir"]
SAMPLES  = cfg["samples"]

os.makedirs(OUT_DIR, exist_ok=True)
sc.settings.verbosity = 1

print(f"====== step_6.py | Load samples | scanpy {sc.__version__} ======")
print(f"  Samples: {len(SAMPLES)}")

meta_path = os.path.join(OUT_DIR, "cell_type_metadata.csv")
if not os.path.exists(meta_path):
    print(f"ERROR: {meta_path} not found. Run step_5.R first.", file=sys.stderr)
    sys.exit(1)

meta_df = pd.read_csv(meta_path)
print(f"  Loaded metadata: {len(meta_df)} spots | "
      f"{meta_df['sample'].nunique()} samples | "
      f"cell types: {sorted(meta_df['cell_type_label'].dropna().unique().tolist())}")


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

    with open(sf_path) as f:
        scalefactors = json.load(f)
    hires_sf = scalefactors["tissue_hires_scalef"]

    adata.obsm["spatial"] = coords.loc[shared, ["pxl_col_in_fullres",
                                                 "pxl_row_in_fullres"]].values.astype(float)

    # stlearn's calc_neighbours reads imagerow/imagecol (hires-scaled pixels)
    adata.obs["imagerow"] = (coords.loc[shared, "pxl_row_in_fullres"].values * hires_sf).astype(float)
    adata.obs["imagecol"] = (coords.loc[shared, "pxl_col_in_fullres"].values * hires_sf).astype(float)

    adata.uns["spatial"] = {sample_name: {
        "scalefactors": scalefactors,
        "use_quality":  "hires",
        "images":       {},
    }}

    adata.obs["sample"] = sample_name

    # Attach community + cell_type_label from Step 5 metadata
    sample_meta = meta_df[meta_df["sample"] == sample_name].set_index("barcode")
    common_bc   = adata.obs_names.intersection(sample_meta.index)
    if len(common_bc) < 10:
        # Try stripping sample suffix from barcodes
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


for samp in SAMPLES:
    print(f"\n  [{samp}] Loading...")
    adata    = load_sample(samp, meta_df)
    out_path = os.path.join(OUT_DIR, f"step6_{samp}.h5ad")
    adata.write_h5ad(out_path)
    print(f"    Saved: step6_{samp}.h5ad")

print(f"\n====== step_6.py COMPLETE | {len(SAMPLES)} samples ======")
