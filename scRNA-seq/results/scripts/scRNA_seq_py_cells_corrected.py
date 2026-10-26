# @@CORRECTED version of the tutorial's Python cells (Part 8).
# @@Every change is tagged  # --- FIX Pn (SC-xx): reason --- ; everything else is verbatim.
# @@Extras the tutorial implies but never writes are tagged  # --- ADD ... --- .
#@@CELL 0
import os
# --- FIX P1 (SC-10): the tutorial filters n_genes < 2000 here while the Seurat half uses
#     nFeature < 2500, so one dataset is described by two thresholds. This switch keeps the
#     tutorial value by default (so the 8 hard-coded names remain valid and the figures stay
#     comparable with the reference set) and can be set to 2500 to reproduce the Seurat half.
N_GENES_MAX = int(os.environ.get("N_GENES_MAX", "2000"))

#@@CELL 54
import pandas as pd
import scanpy as sc

#@@CELL 55
results_file = "pbmc3k.h5ad"  # the file that will store the analysis results
adata = sc.read_10x_mtx(
    "/lustre1/user/bjx131_pkuhpc/scrna_run_20261026/hg19",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file for faster subsequent reading
)
adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
# --- FIX P0 (SC-14): gene symbols repeat in genes.tsv; var_names_make_unique() renames
#     them (ACE -> ACE-1, ...), so gene lookups differ from the Seurat side (which uses
#     gene-1 style suffixes too, but applies them to a different symbol set). Count them. ---
_sym = pd.read_csv("/lustre1/user/bjx131_pkuhpc/scrna_run_20261026/hg19/genes.tsv",
                   sep="\t", header=None)[1]
print(f"loaded: {adata.n_obs} cells x {adata.n_vars} genes; duplicated symbols in genes.tsv: "
      f"{int(_sym.duplicated().sum())}; '_'-containing symbols: {int(_sym.str.contains('_').sum())}")

#@@CELL 57
sc.pl.highest_expr_genes(adata, n_top=20,save="scanpy_highest_expr_genes.png")

#@@CELL 59
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print(f"after filter_cells(min_genes=200)+filter_genes(min_cells=3): "
      f"{adata.n_obs} cells x {adata.n_vars} genes")

# annotate the group of mitochondrial genes as "mt"
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    save="scanpy_violin_plot_QC.png"  # 保存为 violin_plot.png 文件
)

#@@CELL 61
sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", save="scanpy_scatter_pct_counts_mt.png")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", save="scanpy_scatter_n_genes_by_counts.png")
n_before = adata.n_obs
adata = adata[adata.obs.n_genes_by_counts < N_GENES_MAX, :]
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()
print(f"QC (n_genes<{N_GENES_MAX}, pct_mt<5): {n_before} -> {adata.n_obs} cells")
sc.pp.normalize_total(adata, target_sum=1e4)
#Logarithmize the data:
sc.pp.log1p(adata)

#@@CELL 64
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, save="scanpy_HVG.png")
adata.raw = adata.copy()
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
print(f"HVG: {adata.n_vars} genes kept (raw keeps {adata.raw.n_vars})")

#@@CELL 67
# --- FIX P2 (SC-11): `sc.pp.scale` was already applied at the end of the previous cell;
#     running it twice is a no-op on already-scaled data but silently re-clips, so drop it. ---
sc.tl.pca(adata, svd_solver="arpack")
sc.pl.pca(adata, color="CST3",save="scanpy_PCA.png")
# --- FIX P3 (SC-09): literal `log=True,,save=...` in the tutorial is a SyntaxError, which
#     killed the whole cell (PCA + neighbors + UMAP never ran and the ratio plot is missing
#     from the 28 reference figures). One comma removed. ---
sc.pl.pca_variance_ratio(adata, log=True, save="scanpy_PCA_ratio.png")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

#@@CELL 70
sc.tl.leiden(
    adata,
    resolution=0.9,
    random_state=0,
    flavor="igraph",
    n_iterations=2,
    directed=False,
)
print(f"leiden clusters: {adata.obs['leiden'].nunique()} -> {list(adata.obs['leiden'].cat.categories)}")
sc.pl.umap(adata, color=["leiden", "MS4A1", "GNLY","NKG7", "CD3E","IL7R","CD14", "FCER1A","CST3", "FCGR3A", "LYZ", "PPBP","CD8A"],save="scanpy_markers.png")
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False,save="scanpy_rank_genes_groups.png")
sc.pl.rank_genes_groups_violin(adata, groups="0", n_genes=10,save="scanpy_rank_genes_groups_violin.png")

#@@CELL 72
marker_genes = [
    *["IL7R", "CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "CD14"],
    *["LGALS3", "S100A8", "GNLY", "NKG7", "KLRB1"],
    *["FCGR3A", "MS4A7", "FCER1A", "CST3", "PPBP"],
]
sc.pl.dotplot(adata, marker_genes, groupby="leiden",save="scanpy_rank_genes_groups_dotplot.png")
# --- ADD (SC-13): the tutorial keeps the marker statistics inside the AnnData only; export. ---
for i, cl in enumerate(adata.obs["leiden"].cat.categories):
    df_cl = sc.get.rank_genes_groups_df(adata, group=cl)
    df_cl.head(200).to_csv(f"cluster{i}_markers.csv", index=False)
print(f"cluster sizes: {adata.obs['leiden'].value_counts().sort_index().to_dict()}")

#@@CELL 75
new_cluster_names = [
    "CD4 T",
    "CD8 T",
    "B",
    "CD14+ Monocytes",
    "NK",
    "FCGR3A+ Monocytes",
    "Dendritic",
    "Platelet"
]
# --- FIX P4 (SC-12): rename_categories() demands exactly one new name per category and
#     fails with an opaque error otherwise. Check the count and print the audit table first. ---
cats = list(adata.obs["leiden"].cat.categories)
audit = pd.DataFrame({
    "cluster": cats,
    "proposed_name": new_cluster_names if len(cats) == len(new_cluster_names) else None,
    "top5_markers": [", ".join(sc.get.rank_genes_groups_df(adata, group=c)["names"].head(5))
                     for c in cats],
    "n_cells": [int((adata.obs["leiden"] == c).sum()) for c in cats],
})
print(audit.to_string(index=False))
audit.to_csv("cluster_label_audit_py.csv", index=False)
if len(cats) == len(new_cluster_names):
    adata.rename_categories("leiden", new_cluster_names)
    print(f"renamed {len(cats)} clusters with the tutorial's names")
else:
    print(f"WARNING: {len(cats)} leiden clusters but {len(new_cluster_names)} hard-coded names "
          f"-> keeping numeric labels (see cluster_label_audit_py.csv)")
sc.pl.umap(
    adata, color="leiden", legend_loc="on data", title="", frameon=False, save="scanpy_final_umap.png"
)
adata.write("pbmc3k.h5ad", compression="gzip")
adata.obs.to_csv("pbmc3k_obs.csv")
print(f"written: pbmc3k.h5ad ({adata.n_obs} cells), pbmc3k_obs.csv, cluster_label_audit_py.csv")
