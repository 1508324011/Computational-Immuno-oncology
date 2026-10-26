# @@SOURCE class/scRNA-seq/scRNA_seq_analyze_totorial_2026.html (Python cells, verbatim)
# @@NOTE only change: tutorial placeholder "yourfolder/hg19" -> absolute run path
# @@NOTE cell 67 is kept VERBATIM including the `log=True,,save=` double comma (expected SyntaxError)
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

#@@CELL 57
sc.pl.highest_expr_genes(adata, n_top=20,save="scanpy_highest_expr_genes.png")

#@@CELL 59
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

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
adata = adata[adata.obs.n_genes_by_counts < 2000, :]
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
#Logarithmize the data:
sc.pp.log1p(adata)

#@@CELL 64
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, save="scanpy_HVG.png")
adata.raw = adata.copy()
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)

#@@CELL 67
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pl.pca(adata, color="CST3",save="scanpy_PCA.png")
sc.pl.pca_variance_ratio(adata, log=True,,save="scanpy_PCA_ratio.png")
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

adata.rename_categories("leiden", new_cluster_names)
sc.pl.umap(
    adata, color="leiden", legend_loc="on data", title="", frameon=False, save="scanpy_final_umap.png"
)
adata.write("pbmc3k.h5ad", compression="gzip")
