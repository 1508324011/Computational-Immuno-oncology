#!/usr/bin/env python
"""marker_panel_audit.py -- mean log-normalised expression of a canonical PBMC marker
panel per cluster, computed from the delivered pbmc3k.h5ad (scanpy side).

This is the table that decides whether the tutorial's positional renaming
(`rename_categories('leiden', new_cluster_names)`) actually put the right name on
the right cluster.

usage: python marker_panel_audit.py <h5ad> <cluster_col> <out_csv>
"""
import sys

import pandas as pd
import scanpy as sc

h5ad, cluster_col, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]
adata = sc.read_h5ad(h5ad)
raw = adata.raw.to_adata() if adata.raw is not None else adata

panel = {
    "CD4 naive": ["IL7R", "CCR7", "LEF1", "SELL", "TCF7", "MAL"],
    "CD4 memory": ["IL32", "LTB", "CD2", "S100A4"],
    "CD8 T": ["CD8A", "CD8B", "GZMK", "CCL5", "CTSW"],
    "NK": ["GNLY", "GZMB", "NKG7", "KLRF1", "PRF1", "FGFBP2", "SPON2"],
    "B": ["MS4A1", "CD79A", "CD79B", "CD19", "TCL1A"],
    "CD14+ Mono": ["CD14", "LYZ", "S100A8", "S100A9", "FCN1"],
    "FCGR3A+ Mono": ["FCGR3A", "MS4A7", "LST1"],
    "DC": ["FCER1A", "CST3", "CLEC10A", "SERPINF1"],
    "Platelet": ["PPBP", "PF4", "GP9"],
    "technical": ["RPS12", "RPL32", "RPS27", "PTPRC", "ACTB"],
}
genes = [g for v in panel.values() for g in v]
genes = [g for g in genes if g in set(raw.var_names)]
missing = sorted({g for v in panel.values() for g in v} - set(genes))

groups = list(adata.obs[cluster_col].cat.categories)
rows = []
for g in genes:
    vec = raw[:, g].X
    vec = vec.toarray().ravel() if hasattr(vec, "toarray") else vec.ravel()
    rec = {"gene": g,
           "programme": next(k for k, v in panel.items() if g in v)}
    for grp in groups:
        mask = (adata.obs[cluster_col] == grp).to_numpy()
        rec[grp] = round(float(vec[mask].mean()), 4)
    rows.append(rec)
out = pd.DataFrame(rows)
out.to_csv(out_csv, index=False)

print(f"clusters: {groups}")
print(f"cells per cluster: {adata.obs[cluster_col].value_counts().sort_index().to_dict()}")
if missing:
    print(f"markers absent from this object: {missing}")
pd.set_option("display.width", 200, "display.max_columns", 30)
print(out.set_index(["programme", "gene"]).drop(columns=[]).to_string())
