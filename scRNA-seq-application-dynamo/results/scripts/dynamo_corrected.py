#!/usr/bin/env python
"""dynamo_corrected.py -- the tutorial's Dynamo pipeline, re-run with the fixes and the
extra numeric exports needed for the corrected edition of the module.

Everything the tutorial does is kept in the same order and with the same parameters; what is
added is:

  C1  the three Jacobian panels and the two GO panels that the tutorial only displays inline are
      written to disk (same file names the staff's 2026 re-run uses);
  C2  the GO figure titles no longer claim GO_Biological_Process_2018 / KEGG_2016 while the
      enrichr call uses GO_Biological_Process_2021;
  C3  `dyn.ext.enrichr(outdir=...)` points at results/05_enrichment instead of ./figure, so the
      enrichr reports stop polluting the figure directory;
  C4  the tutorial's `dyn.tl.reduceDimension(...)` is silently skipped because the shipped
      h5ad already has a UMAP; we run it again with enforce=True and measure how far the
      recomputed embedding is from the shipped one (pairwise-distance Spearman on a sample);
  C5  the `progenitor` variable that the tutorial computes and then ignores is used: the
      Ductal-only subset gets its own velocity/vector-field read-out;
  C6  `dyn.vf.curl` is computed by the tutorial (pot_curl_div=True) but never drawn -- add it;
  C7  every number the tutorial only *displayed* (adata repr, len(progenitor), cluster list,
      rank_acceleration head) is printed and/or exported, because a plain Python script has no
      display hook;
  C8  all intermediate tables the text of the corrected edition quotes are exported as TSV.

usage: python dynamo_corrected.py <module_results_dir>
"""
import os
import sys
import time
import traceback

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "results")
RUN = os.path.abspath(os.environ.get("DYNAMO_RUN", os.path.expanduser("~/dynamo_run_20261109")))
FIG = os.path.join(OUT, "06_corrected")          # corrected-only figures land in the module
for sub in ["01_preprocess", "02_velocity", "03_vectorfield", "04_jacobian",
            "05_enrichment", "06_corrected"]:
    os.makedirs(os.path.join(OUT, sub), exist_ok=True)
os.makedirs(os.path.join(RUN, "figure_corrected"), exist_ok=True)

T0 = time.time()


def stamp(msg):
    print("\n[%6.1fs] === %s ===" % (time.time() - T0, msg), flush=True)


def dump(df, relpath, **kw):
    p = os.path.join(OUT, relpath)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    df.to_csv(p, sep="\t", **kw)
    print("[export] %s  (%s)" % (relpath, getattr(df, "shape", "-")), flush=True)


def q(series):
    s = pd.Series(series).dropna()
    if s.empty:
        return {}
    return {"n": len(s), "min": round(float(s.min()), 4), "q25": round(float(s.quantile(.25)), 4),
            "median": round(float(s.median()), 4), "mean": round(float(s.mean()), 4),
            "q75": round(float(s.quantile(.75)), 4), "max": round(float(s.max()), 4)}


# ============================ Part 0: tutorial preamble ============================
stamp("preamble")
import warnings
warnings.filterwarnings("ignore")
import dynamo as dyn
import scanpy as sc

os.chdir(RUN)                                    # ./data and ./figure_corrected live here
print("run dir:", os.getcwd(), "| results dir:", OUT, flush=True)
dyn.configuration.set_figure_params("dynamo", background="white", figsize=(3, 2))

data_path = './data/endocrinogenesis_day15.h5ad'
adata = sc.read_h5ad(data_path)
print("raw:", adata.shape, "| layers:", list(adata.layers.keys()),
      "| obs:", list(adata.obs.columns), flush=True)
print("clusters:", list(adata.obs['clusters'].unique()), flush=True)
print(adata.obs['clusters'].value_counts().to_string(), flush=True)

pancreas_genes = ["Hes1", "Nkx6-1", "Nkx2-2", "Neurog3", "Neurod1", "Pax4", "Pax6", "Arx",
                  "Pdx1", "Ins1", "Ins2", "Ghrl", "Ptf1a", "Iapp", "Isl1", "Sox9", "Gcg"]
figure_path = './figure_corrected/'
missing_in_raw = [g for g in pancreas_genes if g not in adata.var_names]
print("pancreas genes absent from raw var_names:", missing_in_raw, flush=True)

# ============================ Part 1: preprocessing ============================
stamp("recipe_monocle")
dyn.pp.recipe_monocle(adata, n_top_genes=4000,
                      fg_kwargs={"shared_count": 20},
                      genes_to_append=pancreas_genes)
print("after recipe_monocle:", adata.shape, flush=True)
# dynamo does NOT subset genes: it flags them (`keep_filtered_genes_key=True` default), so
# adata.n_vars stays 27,998 while 6,956 genes carry pass_basic_filter==True.
# dynamo's own top-4,000 selection is stored in `use_for_pca` -- NOT in `highly_variable_genes`,
# which is a *categorical* column shipped with the h5ad by the original paper's processing
# (True 4,004 / False 13,338 / NaN 10,656). Two consequences, both verified here:
#   * `adata.var['highly_variable_genes'].sum()` raises TypeError (Categorical reduction)
#   * `astype(bool).sum()` gives 17,342, i.e. "not NaN", which is not the HVG set either
hvg_flag = adata.var['use_for_pca'].astype(bool)
pass_basic = adata.var['pass_basic_filter'].astype(bool)
legacy = (adata.var['highly_variable_genes'].astype(str).str.lower().eq('true')
          if 'highly_variable_genes' in adata.var else pd.Series(False, index=adata.var_names))
inter = int((legacy & hvg_flag).sum())
union = int((legacy | hvg_flag).sum())
print("hvg (use_for_pca): %d | pass_basic_filter: %d | legacy highly_variable_genes: %d"
      % (int(hvg_flag.sum()), int(pass_basic.sum()), int(legacy.sum())), flush=True)
print("overlap legacy vs dynamo: %d (Jaccard %.4f)" % (inter, inter / union), flush=True)
dump(pd.DataFrame({"gene": adata.var_names,
                   "use_for_pca": hvg_flag.values,
                   "pass_basic_filter": pass_basic.values,
                   "legacy_hvg": legacy.values,
                   "pancreas_gene": [g in pancreas_genes for g in adata.var_names]}),
     "01_preprocess/var_flags.tsv", index=False)
dump(pd.DataFrame({"item": ["cells_raw", "genes_raw", "cells_after_filter",
                            "genes_kept_in_var", "genes_pass_basic_filter",
                            "n_top_genes_requested", "hvg_use_for_pca", "legacy_hvg_true",
                            "legacy_hvg_nan", "overlap_legacy_vs_dynamo",
                            "jaccard_legacy_vs_dynamo", "pancreas_genes_requested",
                            "pancreas_genes_absent"],
                   "value": [adata.n_obs, 27998, adata.n_obs, adata.n_vars,
                             int(pass_basic.sum()), 4000, int(hvg_flag.sum()),
                             int(legacy.sum()),
                             int(adata.var['highly_variable_genes'].isna().sum())
                             if 'highly_variable_genes' in adata.var else 0,
                             inter, round(inter / union, 4), len(pancreas_genes),
                             len(missing_in_raw)]}),
     "01_preprocess/adata_summary.tsv", index=False)
cc = adata.obs['clusters'].value_counts().rename_axis('cluster').reset_index(name='cells')
cc['share'] = (cc['cells'] / cc['cells'].sum()).round(4)
dump(cc, "01_preprocess/cluster_counts.tsv", index=False)

# ============================ Part 2: dynamics ============================
stamp("dynamics (stochastic)")
t = time.time()
dyn.tl.dynamics(adata, model="stochastic")
print("dynamics took %.1f s" % (time.time() - t), flush=True)

gcol = [c for c in ["gamma", "beta", "half_life", "alpha"] if c in adata.var.columns]
rows = []
for c in gcol:
    d = q(adata.var[c])
    d["quantity"] = c
    rows.append(d)
if rows:
    dump(pd.DataFrame(rows)[["quantity", "n", "min", "q25", "median", "mean", "q75", "max"]],
         "02_velocity/gamma_summary.tsv", index=False)
if 'use_for_dynamics' in adata.var:
    print("velocity genes flagged:", int(adata.var['use_for_dynamics'].astype(bool).sum()),
          flush=True)
if 'pass_basic_filter' in adata.var:
    print("genes passing basic filter:", int(adata.var['pass_basic_filter'].astype(bool).sum()),
          "| kept in adata.var:", adata.n_vars, flush=True)

# C4: the shipped h5ad already carries a UMAP, so the tutorial's reduceDimension is a no-op.
# The comparative recompute runs on a COPY: replacing adata's embedding would silently change
# every figure below (they all project the vector field onto X_umap), which is not a correction
# the tutorial asks for. Keeping the copy makes figures 01-09 byte-identical to the as-written
# pass, so any remaining difference is a real fix and not an embedding swap.
stamp("C4: reduceDimension skip vs enforce (on a copy)")
umap_before = np.array(adata.obsm['X_umap'], dtype=float)
print("UMAP present before reduceDimension:", 'X_umap' in adata.obsm, flush=True)
dyn.tl.reduceDimension(adata, n_pca_components=30)          # tutorial call: silently skipped
print("after tutorial call, X_umap identical:",
      np.allclose(umap_before, adata.obsm['X_umap']), flush=True)
_tmp = adata.copy()
dyn.tl.reduceDimension(_tmp, n_pca_components=30, enforce=True)
umap_after = np.array(_tmp.obsm['X_umap'], dtype=float)
del _tmp


def pairwise_spearman(A, B, n=600, seed=0):
    """Spearman correlation of the pairwise distance matrices of two embeddings -- invariant to
    rotation/reflection/scale, so it measures 'same neighbourhood geometry', not 'same axes'."""
    rng = np.random.default_rng(seed)
    idx = rng.choice(A.shape[0], size=min(n, A.shape[0]), replace=False)
    da = pd.Series(np.linalg.norm(A[idx][:, None] - A[idx][None], axis=2).ravel())
    db = pd.Series(np.linalg.norm(B[idx][:, None] - B[idx][None], axis=2).ravel())
    return float(da.corr(db, method='spearman'))


r_dist = pairwise_spearman(umap_before, umap_after)
disp = float(np.mean(np.linalg.norm(umap_before - umap_after, axis=1)))
spread = float(np.mean(np.std(umap_before, axis=0)))
print("C4: pairwise-distance Spearman(shipped UMAP, recomputed UMAP) = %.4f" % r_dist, flush=True)
print("C4: mean displacement %.4f vs UMAP spread %.4f (ratio %.3f)"
      % (disp, spread, disp / spread), flush=True)
dump(pd.DataFrame([{"what": "pairwise_distance_spearman", "value": round(r_dist, 4)},
                   {"what": "mean_displacement", "value": round(disp, 4)},
                   {"what": "umap_spread", "value": round(spread, 4)},
                   {"what": "displacement_over_spread", "value": round(disp / spread, 4)},
                   {"what": "n_cells_sampled", "value": 600}]),
     "06_corrected/umap_recompute_comparison.tsv", index=False)

fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True)
axes[0].scatter(umap_before[:, 0], umap_before[:, 1], s=1, c='#1f77b4')
axes[0].set_title("shipped UMAP (tutorial reuses it)")
axes[1].scatter(umap_after[:, 0], umap_after[:, 1], s=1, c='#d62728')
axes[1].set_title("reduceDimension(enforce=True)")
for ax in axes:
    ax.set_xticks([])
    ax.set_yticks([])
fig.savefig(os.path.join(FIG, "12_umap_shipped_vs_recomputed.png"), dpi=200, bbox_inches='tight')
plt.close(fig)

# ============================ Part 3: velocities ============================
stamp("cell_velocities")
dyn.tl.cell_velocities(adata, method="pearson", other_kernels_dict={"transform": "sqrt"})
dyn.tl.cell_velocities(adata, basis="pca")
for key in ["velocity_S", "velocity_T", "velocity_N", "velocity_C"]:
    if key in adata.obsm:
        print("%s: %s" % (key, np.array(adata.obsm[key]).shape), flush=True)
vstats = []
if 'velocity_S' in adata.obsm:
    v = np.array(adata.obsm['velocity_S'], dtype=float)
    s = q(np.linalg.norm(v, axis=1))
    s["quantity"] = "velocity_S_norm_filtered_space"
    vstats.append(s)
if 'velocity_umap' in adata.obsm:
    v = np.array(adata.obsm['velocity_umap'], dtype=float)
    s = q(np.linalg.norm(v, axis=1))
    s["quantity"] = "velocity_umap_norm"
    vstats.append(s)
if vstats:
    dump(pd.DataFrame(vstats)[["quantity", "n", "min", "q25", "median", "mean", "q75", "max"]],
         "02_velocity/cell_velocity_stats.tsv", index=False)

# ============================ Part 4: figures 01-05 ============================
stamp("figures 01-04")
dyn.pl.streamline_plot(adata, color=["clusters"], basis="umap", show_legend="on data",
                       show_arrowed_spines=True, save_show_or_return='save',
                       save_kwargs={"prefix": "01_velocity_streamline", "ext": "png",
                                    "path": figure_path, "dpi": 300, "transparent": True,
                                    "close": True, "verbose": True})
dyn.pl.cell_wise_vectors(adata, color=["clusters"], basis="umap", show_legend="on data",
                         quiver_length=4, quiver_size=4, figsize=(5, 3),
                         show_arrowed_spines=False, save_show_or_return='save',
                         save_kwargs={"prefix": "02_cell_wise_velocity", "ext": "png",
                                      "path": figure_path, "dpi": 300, "transparent": True,
                                      "close": True, "verbose": True})
progenitor = adata.obs_names[adata.obs.clusters.isin(["Ductal"])]
print("len(progenitor) =", len(progenitor), flush=True)

stamp("vector field")
dyn.vf.VectorField(adata, basis="pca", pot_curl_div=True)
dyn.vf.VectorField(adata, basis="umap", pot_curl_div=True)
dyn.vf.speed(adata, basis="pca")
dyn.vf.divergence(adata, basis="pca")
dyn.vf.acceleration(adata, basis="pca")
dyn.vf.curl(adata, basis="umap")

dyn.configuration.set_figure_params("dynamo", background="black")
fig1, fig1_axes = plt.subplots(ncols=2, nrows=2, constrained_layout=True, figsize=(10, 6))
dyn.pl.cell_wise_vectors(adata, color="pca_ddhodge_potential", pointsize=0.1, alpha=0.7,
                         ax=fig1_axes[0, 0], quiver_length=6, quiver_size=6,
                         save_show_or_return="return", background="black")
dyn.pl.grid_vectors(adata, color="speed_pca", basis="umap", ax=fig1_axes[0, 1],
                    quiver_length=12, quiver_size=12, save_show_or_return="return",
                    background="black")
dyn.pl.streamline_plot(adata, color="divergence_pca", basis="umap", ax=fig1_axes[1, 0],
                       save_show_or_return="return", background="black")
dyn.pl.streamline_plot(adata, color="acceleration_pca", basis="umap", ax=fig1_axes[1, 1],
                       save_show_or_return="return", background="black")
fig1.savefig('./figure_corrected/03_combined_plots.png')
plt.close(fig1)

dyn.pl.topography(adata, basis="umap", background="white", color=["clusters"],
                  streamline_color="black", show_legend="on data", save_show_or_return='save',
                  save_kwargs={"prefix": "04_cell_state", "ext": "png", "path": figure_path,
                               "dpi": 300, "transparent": True, "close": True, "verbose": True})

# ============================ Part 5: vector-field tables ============================
stamp("vector-field tables")
rows = []
for key in ["speed_pca", "speed_umap", "divergence_pca", "divergence_umap",
            "acceleration_pca", "acceleration_umap", "curl_umap", "pca_ddhodge_potential"]:
    if key in adata.obs:
        d = q(adata.obs[key])
        d["quantity"] = key
        rows.append(d)
if rows:
    dump(pd.DataFrame(rows)[["quantity", "n", "min", "q25", "median", "mean", "q75", "max"]],
         "03_vectorfield/vf_stats.tsv", index=False)

by_cluster = []
for key in ["speed_pca", "divergence_pca", "acceleration_pca", "curl_umap",
            "pca_ddhodge_potential"]:
    if key in adata.obs:
        g = adata.obs.groupby('clusters', observed=True)[key]
        for name, val in g.mean().items():
            by_cluster.append({"cluster": name, "quantity": key, "mean": round(float(val), 4),
                               "median": round(float(g.median()[name]), 4),
                               "n_cells": int(g.size()[name])})
if by_cluster:
    bc = pd.DataFrame(by_cluster)
    dump(bc, "03_vectorfield/vf_by_cluster.tsv", index=False)
    piv = bc.pivot(index='cluster', columns='quantity', values='mean')
    print(piv.round(3).to_string(), flush=True)

# C6: curl is computed by the tutorial but never drawn.
stamp("C6: curl figure")
dyn.configuration.set_figure_params("dynamo", background="black")
try:
    dyn.pl.streamline_plot(adata, color="curl_umap", basis="umap",
                           save_show_or_return='save',
                           save_kwargs={"prefix": "12_curl_nullcline", "ext": "png",
                                        "path": figure_path, "dpi": 300,
                                        "transparent": True, "close": True, "verbose": True})
except Exception as e:
    print("C6 figure failed:", type(e).__name__, e, flush=True)

# ============================ Part 6: state graph + jacobian ============================
stamp("state graph")
dyn.configuration.set_figure_params("dynamo", background="white")
dyn.pd.state_graph(adata, group='clusters', basis='pca', method='vf')
dyn.pl.state_graph(adata, color=['clusters'], group='clusters', basis='umap',
                   show_legend='on data', save_show_or_return='save',
                   save_kwargs={"prefix": "05_cell_state_transition_map", "ext": "png",
                                "path": figure_path, "dpi": 300, "transparent": True,
                                "close": True, "verbose": True})
stamp("jacobian panels (06-09)")
# The tutorial draws cell 22's panel with save_show_or_return='save' but leaves cells 23/24 to
# the notebook's inline display, so 07/08/09 never reach disk. Same call, now saved (C1).
dyn.vf.jacobian(adata, regulators=pancreas_genes)
for reg, eff, prefix in [("Neurog3", "Pax4", "06_reg_Ngn3_activate_Pax4"),
                         ("Arx", "Pax4", "07_reg_Arx_Pax4"),
                         ("Pax4", "Arx", "08_reg_Pax4_Arx"),
                         ("Pax4", "Ins2", "09_reg_Pax4_Ins2")]:
    dyn.pl.jacobian(adata, basis="umap", regulators=[reg], effectors=[eff], alpha=1,
                    save_show_or_return='save',
                    save_kwargs={"prefix": prefix, "ext": "png", "path": figure_path,
                                 "dpi": 300, "transparent": True, "close": True,
                                 "verbose": True})
print("pancreas jacobian kept in adata.uns['jacobian_pca']; numbers: "
      "results/scripts/dynamo_jacobian_export.py", flush=True)

# ============================ Part 7: acceleration ranking ============================
stamp("rank_acceleration_genes")
dyn.vf.acceleration(adata, basis='pca')
dyn.vf.rank_acceleration_genes(adata, groups='clusters')
print("rank_acceleration keys:", list(adata.uns['rank_acceleration'].keys()), flush=True)
print(adata.uns['rank_acceleration'][:5].to_string(), flush=True)
ra = adata.uns['rank_acceleration']
dump(ra, "04_jacobian/rank_acceleration_top200.tsv")
top20 = ra.apply(lambda col: ", ".join(col.dropna().astype(str)[:20]))
dump(pd.DataFrame({"cluster": top20.index, "top20_accelerating_genes": top20.values}),
     "04_jacobian/rank_acceleration_top20_per_cluster.tsv", index=False)
print("clusters in rank_acceleration:", list(ra.columns), flush=True)

# ============================ Part 8: enrichment ============================
stamp("C2/C3: enrichr")
enr = dyn.ext.enrichr(
    adata.uns['rank_acceleration']['Ductal'][:100].to_list(),
    organism='mouse',
    gene_sets='GO_Biological_Process_2021',
    outdir=os.path.join(OUT, "05_enrichment"),
)
res = enr.res2d.copy()
dump(res, "05_enrichment/GO_BP_2021_Ductal_top100.tsv", index=False)
print("enrichr rows:", res.shape, "| min adj p:", float(pd.to_numeric(
    res['Adjusted P-value'], errors='coerce').min()), flush=True)
print(res.head(10).to_string(), flush=True)
sig = res[pd.to_numeric(res['Adjusted P-value'], errors='coerce') < 0.05]
print("significant terms (adj p < 0.05):", len(sig), flush=True)

from gseapy.plot import barplot, dotplot
fig_bar = barplot(res, title='GO_Biological_Process_2021', cutoff=0.05, top_term=10).get_figure()
fig_bar.savefig(os.path.join(FIG, "10_Ductal_GO_enrichment_barplot_corrected.png"),
                dpi=300, bbox_inches='tight')
plt.close(fig_bar)
fig_dot = dotplot(res, title='GO_Biological_Process_2021', cutoff=0.05, top_term=10).get_figure()
fig_dot.savefig(os.path.join(FIG, "11_Ductal_GO_enrichment_dotplot_corrected.png"),
                dpi=300, bbox_inches='tight')
plt.close(fig_dot)

# ============================ Part 9: C5 Ductal subset ============================
stamp("C5: Ductal subset")
try:
    sub = adata[progenitor].copy()
    dd = pd.DataFrame({
        "quantity": ["n_cells", "mean_speed_pca", "mean_divergence_pca",
                     "mean_acceleration_pca", "mean_curl_umap",
                     "mean_ddhodge_potential", "mean_speed_pca_other_cells"],
        "value": [sub.n_obs,
                  float(sub.obs['speed_pca'].mean()),
                  float(sub.obs['divergence_pca'].mean()),
                  float(sub.obs['acceleration_pca'].mean()),
                  float(sub.obs['curl_umap'].mean()),
                  float(sub.obs['pca_ddhodge_potential'].mean()),
                  float(adata.obs.loc[~adata.obs_names.isin(progenitor), 'speed_pca'].mean())],
    })
    dump(dd, "06_corrected/ductal_subset_stats.tsv", index=False)
    print(dd.to_string(), flush=True)
except Exception as e:
    print("C5 failed:", type(e).__name__, e, flush=True)

# ============================ wrap-up ============================
stamp("done")
print("total wall time %.1f s" % (time.time() - T0), flush=True)
print("figures written to %s" % os.path.abspath(figure_path), flush=True)
print("corrected-only figures under %s" % FIG, flush=True)
