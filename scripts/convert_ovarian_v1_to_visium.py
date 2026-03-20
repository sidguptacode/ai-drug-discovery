#!/usr/bin/env python3
"""
Convert ovarian_v1 flat MTX/counts in .data/ovarian_v1 to Visium-style layout:
one folder per sample with <sample>_filtered_feature_bc_matrix.h5,
<sample>_tissue_positions_list.csv, <sample>_scalefactors_json.json.
Writes to an output directory (default: outputs/ovarian_v1_visium_format).
Does not modify .data (read-only).
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / ".data" / "ovarian_v1"
OUT_BASE = REPO_ROOT / "outputs" / "ovarian_v1_visium_format"

# MTX samples: (file prefix, sample_id for folder name)
MTX_SAMPLES = [
    ("GSM6506110_SP1", "GSM6506110_SP1"),
    ("GSM6506111_SP2", "GSM6506111_SP2"),
    ("GSM6506112_SP3", "GSM6506112_SP3"),
    ("GSM6506113_SP4", "GSM6506113_SP4"),
    ("GSM6506114_SP5", "GSM6506114_SP5"),
    ("GSM6506115_SP6", "GSM6506115_SP6"),
    ("GSM6506116_SP7", "GSM6506116_SP7"),
    ("GSM6506117_SP8", "GSM6506117_SP8"),
]


def write_10x_h5(out_path: Path, matrix: sp.spmatrix, barcodes: list[str], feature_ids: list[str], feature_names: list[str], feature_type: str = "Gene Expression") -> None:
    """Write matrix in 10x H5 filtered_feature_bc_matrix format. Shape (n_features, n_barcodes)."""
    # 10x stores (features x barcodes) in CSC
    if matrix.shape[0] != len(feature_ids) or matrix.shape[1] != len(barcodes):
        matrix = sp.csc_matrix(matrix)
    else:
        matrix = sp.csc_matrix(matrix)
    n_f, n_b = matrix.shape
    with h5py.File(out_path, "w") as f:
        m = f.create_group("matrix")
        m.create_dataset("shape", data=np.array([n_f, n_b], dtype=np.int64))
        m.create_dataset("data", data=matrix.data.astype(np.float32))
        m.create_dataset("indices", data=matrix.indices.astype(np.int32))
        m.create_dataset("indptr", data=matrix.indptr.astype(np.int64))
        dt = h5py.special_dtype(vlen=str)
        m.create_dataset("barcodes", data=np.array(barcodes, dtype=object), dtype=dt)
        feats = m.create_group("features")
        feats.create_dataset("id", data=np.array(feature_ids, dtype=object), dtype=dt)
        feats.create_dataset("name", data=np.array(feature_names, dtype=object), dtype=dt)
        feats.create_dataset("feature_type", data=np.array([feature_type] * n_f, dtype=object), dtype=dt)
        feats.create_dataset("genome", data=np.array(["GRCh38"] * n_f, dtype=object), dtype=dt)


def convert_sample(prefix: str, sample_id: str, out_dir: Path) -> None:
    data_dir = DATA_DIR
    mtx_path = data_dir / f"{prefix}_matrix.mtx"
    barcodes_path = data_dir / f"{prefix}_barcodes.tsv"
    features_path = data_dir / f"{prefix}_features.tsv"
    if not mtx_path.exists():
        print(f"  Skip {sample_id}: {mtx_path} not found", file=sys.stderr)
        return
    matrix = sio.mmread(mtx_path)
    if matrix.shape[0] != matrix.shape[1]:
        # MTX from Space Ranger is (genes x barcodes)
        pass
    matrix = sp.csc_matrix(matrix)
    n_genes, n_cells = matrix.shape
    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")[0].astype(str).tolist()
    if len(barcodes) != n_cells:
        # sometimes MTX is (cells x genes)
        if matrix.shape[0] == len(barcodes):
            matrix = matrix.T
            n_genes, n_cells = matrix.shape
        else:
            raise ValueError(f"{sample_id}: barcodes {len(barcodes)} vs matrix {matrix.shape}")
    feats = pd.read_csv(features_path, header=None, sep="\t")
    if feats.shape[0] != n_genes:
        if matrix.shape[1] == feats.shape[0]:
            matrix = matrix.T
            n_genes, n_cells = matrix.shape
        else:
            raise ValueError(f"{sample_id}: features {feats.shape[0]} vs matrix {matrix.shape}")
    feature_ids = feats[0].astype(str).tolist()
    feature_names = (feats[1] if feats.shape[1] > 1 else feats[0]).astype(str).tolist()
    sample_out = out_dir / sample_id
    sample_out.mkdir(parents=True, exist_ok=True)
    h5_path = sample_out / f"{sample_id}_filtered_feature_bc_matrix.h5"
    write_10x_h5(h5_path, matrix, barcodes, feature_ids, feature_names)
    # Dummy spatial: grid layout so stlearn has coordinates
    n_spots = len(barcodes)
    n_cols = max(1, int(np.ceil(np.sqrt(n_spots))))
    n_rows = max(1, (n_spots + n_cols - 1) // n_cols)
    positions = []
    for i, bc in enumerate(barcodes):
        row, col = i // n_cols, i % n_cols
        positions.append({
            "barcode": bc,
            "in_tissue": 1,
            "array_row": row,
            "array_col": col,
            "pxl_row": row,
            "pxl_col": col,
            "pxl_row_in_fullres": float(row),
            "pxl_col_in_fullres": float(col),
        })
    pos_df = pd.DataFrame(positions)
    pos_path = sample_out / f"{sample_id}_tissue_positions_list.csv"
    pos_df.to_csv(pos_path, index=False)
    sf_path = sample_out / f"{sample_id}_scalefactors_json.json"
    sf_path.write_text(json.dumps({
        "tissue_hires_scalef": 1.0,
        "spot_diameter_fullres": 55
    }))
    print(f"  {sample_id}: {n_cells} spots, {n_genes} genes -> {sample_out}")


def main():
    out_dir = Path(os.environ.get("OVARIAN_V1_VISIUM_OUT", str(OUT_BASE)))
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Converting ovarian_v1 MTX -> {out_dir}")
    for prefix, sample_id in MTX_SAMPLES:
        convert_sample(prefix, sample_id, out_dir)
    print("Done. Use in config: data_dir:", out_dir, "samples:", [s[1] for s in MTX_SAMPLES])


if __name__ == "__main__":
    main()
