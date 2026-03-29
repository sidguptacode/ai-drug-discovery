#!/usr/bin/env python3
"""
step_6.py — Load samples and attach cell-type metadata
Reads:  OUT_DIR/cell_type_metadata.csv   (from step_5.R)
        (h5)  DATA_DIR/<sample>/<sample>_filtered_feature_bc_matrix.h5
              DATA_DIR/<sample>/<sample>_tissue_positions_list.csv
              DATA_DIR/<sample>/<sample>_scalefactors_json.json
        (mtx, flat)   DATA_DIR/<sample>_{matrix.mtx,barcodes.tsv,features.tsv}
                      DATA_DIR/<sample>_tissue_positions_list.csv
                      DATA_DIR/<sample>_scalefactors_json.json
        (mtx, subdir) DATA_DIR/<sample>_{matrix.mtx,barcodes.tsv,features.tsv}
                      DATA_DIR/<sample>_spatial/spatial/tissue_positions_list.csv
                      DATA_DIR/<sample>_spatial/spatial/scalefactors_json.json
Writes: OUT_DIR/step6_<sample>.h5ad      (one per sample)

Usage:  python step_6.py [config.yml]
        Falls back to config_dipg.yml if no argument given.
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

DATA_DIR       = cfg["data_dir"]
OUT_DIR        = cfg["out_dir"]
SAMPLES        = cfg["samples"]
DATA_FORMAT    = cfg.get("data_format",    "h5")
SPATIAL_LAYOUT = cfg.get("spatial_layout", "subdir")

os.makedirs(OUT_DIR, exist_ok=True)
sc.settings.verbosity = 1

print(f"====== step_6.py | Load samples | scanpy {sc.__version__} ======")
print(f"  format={DATA_FORMAT} | spatial={SPATIAL_LAYOUT} | samples={len(SAMPLES)}")

meta_path = os.path.join(OUT_DIR, "cell_type_metadata.csv")
if not os.path.exists(meta_path):
    print(f"ERROR: {meta_path} not found. Run step_5.R first.", file=sys.stderr)
    sys.exit(1)

meta_df = pd.read_csv(meta_path)
print(f"  Loaded metadata: {len(meta_df)} spots | "
      f"{meta_df['sample'].nunique()} samples | "
      f"cell types: {sorted(meta_df['cell_type_label'].dropna().unique().tolist())}")


def get_spatial_paths(samp):
    """Return (coord_path, sf_path) for a sample, mirroring step_1.R logic."""
    if DATA_FORMAT == "h5":
        return (
            os.path.join(DATA_DIR, samp, f"{samp}_tissue_positions_list.csv"),
            os.path.join(DATA_DIR, samp, f"{samp}_scalefactors_json.json"),
        )
    elif SPATIAL_LAYOUT == "subdir":
        # ovarian (GSE211956): spatial in {sample}_spatial/spatial/
        return (
            os.path.join(DATA_DIR, f"{samp}_spatial", "spatial",
                         "tissue_positions_list.csv"),
            os.path.join(DATA_DIR, f"{samp}_spatial", "spatial",
                         "scalefactors_json.json"),
        )
    else:
        # flat: all files alongside MTX with {sample}_ prefix
        return (
            os.path.join(DATA_DIR, f"{samp}_tissue_positions_list.csv"),
            os.path.join(DATA_DIR, f"{samp}_scalefactors_json.json"),
        )


def read_tissue_positions(coord_path):
    """Read tissue_positions_list robustly (with or without header row).
    Always returns a DataFrame indexed by barcode with columns
    pxl_row_in_fullres and pxl_col_in_fullres, filtered to in_tissue==1."""
    try:
        df = pd.read_csv(coord_path, header=None,
                         names=["barcode", "in_tissue", "array_row",
                                "array_col", "pxl_row", "pxl_col"])
        if df["in_tissue"].dtype == object:
            raise ValueError("header row detected")
    except Exception:
        df = pd.read_csv(coord_path)
        df.columns = [c.lower() for c in df.columns]

    df = df[df["in_tissue"] == 1].set_index("barcode")

    # Normalise pixel column names to the fullres variant used downstream
    if "pxl_row_in_fullres" not in df.columns:
        df = df.rename(columns={"pxl_row": "pxl_row_in_fullres",
                                 "pxl_col": "pxl_col_in_fullres"})
    return df


def load_counts(samp):
    """Load count matrix from h5 or MTX depending on DATA_FORMAT."""
    if DATA_FORMAT == "h5":
        adata = sc.read_10x_h5(
            os.path.join(DATA_DIR, samp, f"{samp}_filtered_feature_bc_matrix.h5"))
    else:
        mtx_path  = os.path.join(DATA_DIR, f"{samp}_matrix.mtx")
        bc_path   = os.path.join(DATA_DIR, f"{samp}_barcodes.tsv")
        feat_path = os.path.join(DATA_DIR, f"{samp}_features.tsv")
        adata = sc.read_mtx(mtx_path).T   # MTX is genes × cells; transpose to cells × genes
        adata.obs_names = pd.read_csv(bc_path, header=None, sep="\t")[0].values
        features        = pd.read_csv(feat_path, header=None, sep="\t")
        adata.var_names = features[1].values   # column 1 = gene symbol
    adata.var_names_make_unique()
    return adata


def load_sample(sample_name: str, meta_df: pd.DataFrame) -> sc.AnnData:
    coord_path, sf_path = get_spatial_paths(sample_name)

    adata  = load_counts(sample_name)
    coords = read_tissue_positions(coord_path)

    # Intersect: keeps only barcodes present in both matrix and on-tissue positions
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
        strip           = lambda bc: bc.rsplit("_", 1)[0]
        a_map           = {strip(b): b for b in adata.obs_names}
        m_map           = {strip(b): b for b in sample_meta.index}
        common_stripped = set(a_map) & set(m_map)
        common_bc       = pd.Index([a_map[b] for b in common_stripped])
        meta_idx        = pd.Index([m_map[b] for b in common_stripped])
        sample_meta     = sample_meta.loc[meta_idx]
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
