#!/usr/bin/env python
"""cluster_correspondence.py -- cross-tool contingency table: how the Seurat 9-cluster
solution and the scanpy leiden solution partition the same barcodes.
usage: cluster_correspondence.py <pbmc3k.h5ad> <seurat_meta.csv> <out_prefix>
"""
import sys
import pandas as pd

h5ad, meta_csv, prefix = sys.argv[1], sys.argv[2], sys.argv[3]
import scanpy as sc
adata = sc.read_h5ad(h5ad)
seu = pd.read_csv(meta_csv)
obs = adata.obs[["leiden"]].copy()
obs["barcode"] = obs.index
m = obs.merge(seu, on="barcode", how="inner")
ct = pd.crosstab(m["seurat_label"], m["leiden"])
ct.to_csv(f"{prefix}_counts.csv")
share = ct.div(ct.sum(axis=0), axis=1).round(3)          # share of each scanpy cluster
share.to_csv(f"{prefix}_share_by_scanpy_cluster.csv")
share_r = ct.div(ct.sum(axis=1), axis=0).round(3)        # share of each Seurat label
share_r.to_csv(f"{prefix}_share_by_seurat_label.csv")
print(f"barcodes matched: {len(m)} of {adata.n_obs} (scanpy) / {len(seu)} (Seurat)")
print("\ncounts (rows = Seurat label, cols = scanpy cluster):")
print(ct.to_string())
print("\nshare of each scanpy cluster (rows sum to 1 over the column):")
print(share.to_string())
