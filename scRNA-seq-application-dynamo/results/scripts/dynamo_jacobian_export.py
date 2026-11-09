#!/usr/bin/env python
"""dynamo_jacobian_export.py -- export the numbers behind the tutorial's Jacobian panels.

Why a separate script: `dyn.pl.jacobian` only *draws*; the per-cell Jacobian values it colours by
are never written to disk by the tutorial (cells 21-24), and `dyn.pd.state_graph` stores its
transition data in a nested dict that pandas cannot take directly. This script re-runs the
analysis up to the vector field and dumps both as plain TSVs, so the corrected edition can quote
real numbers instead of eyeballing colours.

Shape semantics (verified against dynamo 1.3.2 source, vectorfield/utils.py:497):
    adata.uns['jacobian_pca']['jacobian']       -> (n_PCs, n_PCs, n_cells), PCA-space
    adata.uns['jacobian_pca']['jacobian_gene']  -> (n_effector, n_regulator, n_cells), because
        subset_jacobian_transformation computes  Qi @ J @ Qj.T  with Qi = effector loadings and
        Qj = regulator loadings, i.e. entry [i, j] = d(effector i)/d(regulator j).
So the gene-level tensor is read as Jg[eff_idx, reg_idx, cell].

Also covered: the pair the tutorial's own comment names but its code never plots,
Pdx1 -> Ins2 ("Pdx1 activates Ins2 in beta cells", cell 24 comment).

usage: DYNAMO_RUN=~/dynamo_run_20261109 python dynamo_jacobian_export.py <module_results_dir>
"""
import json
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
import dynamo as dyn
import scanpy as sc

OUT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "results")
RUN = os.path.abspath(os.environ.get("DYNAMO_RUN", os.path.expanduser("~/dynamo_run_20261109")))
os.chdir(RUN)
T0 = time.time()

PANCREAS = ["Hes1", "Nkx6-1", "Nkx2-2", "Neurog3", "Neurod1", "Pax4", "Pax6", "Arx", "Pdx1",
            "Ins1", "Ins2", "Ghrl", "Ptf1a", "Iapp", "Isl1", "Sox9", "Gcg"]


def stamp(msg):
    print("\n[%6.1fs] === %s ===" % (time.time() - T0, msg), flush=True)


def dump(df, relpath, **kw):
    p = os.path.join(OUT, relpath)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    df.to_csv(p, sep="\t", **kw)
    print("[export] %s  (%s)" % (relpath, getattr(df, "shape", "-")), flush=True)


stamp("pipeline up to the vector field")
dyn.configuration.set_figure_params("dynamo", background="white", figsize=(3, 2))
adata = sc.read_h5ad('./data/endocrinogenesis_day15.h5ad')
dyn.pp.recipe_monocle(adata, n_top_genes=4000, fg_kwargs={"shared_count": 20},
                      genes_to_append=PANCREAS)
dyn.tl.dynamics(adata, model="stochastic")
dyn.tl.reduceDimension(adata, n_pca_components=30)
dyn.tl.cell_velocities(adata, method="pearson", other_kernels_dict={"transform": "sqrt"})
dyn.tl.cell_velocities(adata, basis="pca")
dyn.vf.VectorField(adata, basis="pca", pot_curl_div=True)
dyn.vf.VectorField(adata, basis="umap", pot_curl_div=True)
dyn.vf.speed(adata, basis="pca")
dyn.vf.divergence(adata, basis="pca")
dyn.vf.acceleration(adata, basis="pca")
dyn.vf.curl(adata, basis="umap")

# ---------------------------------------------------------------- state graph
stamp("state graph transition matrix")
dyn.pd.state_graph(adata, group='clusters', basis='pca', method='vf')
graph = adata.uns.get('clusters_graph')
print("type(clusters_graph):", type(graph), flush=True)
labels = list(adata.obs['clusters'].cat.categories) if hasattr(
    adata.obs['clusters'], 'cat') else sorted(adata.obs['clusters'].unique())
if isinstance(graph, dict):
    print("clusters_graph keys:", [(k, np.asarray(v).shape) for k, v in graph.items()],
          flush=True)
    for k, v in graph.items():
        try:
            arr = np.asarray(v, dtype=float)
        except (ValueError, TypeError):
            print("skip clusters_graph[%r]: non-numeric (%s)" % (k, type(v).__name__), flush=True)
            continue
        if arr.ndim != 2 or arr.shape[0] != len(labels):
            continue
        df = pd.DataFrame(arr, index=labels, columns=labels).round(4)
        dump(df, "03_vectorfield/clusters_graph_%s.tsv" % k)
        print(df.to_string(), flush=True)
        if arr.shape[1] == len(labels):
            rowsum = df.sum(axis=1)
            dump(pd.DataFrame({"cluster": labels, "row_sum": rowsum.round(4).values}),
                 "03_vectorfield/clusters_graph_rowsum_%s.tsv" % k, index=False)

# ---------------------------------------------------------------- jacobian
stamp("jacobian (all pancreas genes at once)")
dyn.vf.jacobian(adata, regulators=PANCREAS, effectors=PANCREAS)
ret = adata.uns['jacobian_pca']
regs = list(ret['regulators'])
effs = list(ret['effectors'])
J_pc = np.asarray(ret['jacobian'], dtype=float)
Jg = np.asarray(ret['jacobian_gene'], dtype=float)
cell_idx = np.asarray(ret['cell_idx'])
print("jacobian (PCA space):", J_pc.shape, "| jacobian_gene:", Jg.shape,
      "| n_reg:", len(regs), "| n_eff:", len(effs), "| n_cells:", len(cell_idx), flush=True)
with open(os.path.join(OUT, "04_jacobian", "jacobian_metadata.json"), "w") as f:
    json.dump({"jacobian_pca_shape": list(J_pc.shape),
               "jacobian_gene_shape": list(Jg.shape),
               "regulators": regs, "effectors": effs,
               "n_cells": int(len(cell_idx)),
               "axis_note": "jacobian_gene[i, j, c] = d(effector i)/d(regulator j) in cell c "
                            "(dynamo vectorfield/utils.py: subset_jacobian_transformation, "
                            "Qi=effector loadings, Qj=regulator loadings)"}, f, indent=2)

clusters = adata.obs['clusters'].values[cell_idx] if len(cell_idx) else adata.obs['clusters'].values

# gene-level matrix averaged over all cells: rows = regulators, cols = effectors (reporting
# orientation -- the tensor itself is (effector, regulator, cell))
mean_all = Jg.mean(axis=2)                                  # (n_eff, n_reg)
dump(pd.DataFrame(np.round(mean_all.T, 5), index=regs, columns=effs),
     "04_jacobian/jacobian_mean_all_cells.tsv")

is_ductal = np.asarray(clusters) == 'Ductal'
if is_ductal.any():
    mean_ductal = Jg[:, :, is_ductal].mean(axis=2)
    dump(pd.DataFrame(np.round(mean_ductal.T, 5), index=regs, columns=effs),
         "04_jacobian/jacobian_mean_ductal_cells.tsv")

PAIRS = [("Neurog3", "Pax4"), ("Arx", "Pax4"), ("Pax4", "Arx"), ("Pax4", "Ins2"),
         ("Pdx1", "Ins2"), ("Neurog3", "Neurod1"), ("Pax4", "Pax6"), ("Pdx1", "Pax4"),
         ("Nkx6-1", "Ins2"), ("Gcg", "Arx")]
rows, per_cluster = [], []
for r_, e_ in PAIRS:
    if r_ not in regs or e_ not in effs:
        print("skip %s->%s (not in jacobian)" % (r_, e_), flush=True)
        continue
    v = Jg[effs.index(e_), regs.index(r_), :]
    s = pd.Series(v)
    rows.append({"regulator": r_, "effector": e_, "n_cells": int(s.size),
                 "min": round(float(s.min()), 5), "q25": round(float(s.quantile(.25)), 5),
                 "median": round(float(s.median()), 5), "mean": round(float(s.mean()), 5),
                 "q75": round(float(s.quantile(.75)), 5), "max": round(float(s.max()), 5),
                 "frac_positive": round(float((s > 0).mean()), 4)})
    dump(pd.DataFrame({"cell": adata.obs_names.values[cell_idx], "clusters": clusters,
                       "jacobian": v}), "04_jacobian/jac_%s_%s.tsv" % (r_, e_), index=False)
    for name, vals in pd.Series(v).groupby(pd.Series(clusters)):
        per_cluster.append({"regulator": r_, "effector": e_, "cluster": name,
                            "n_cells": int(vals.size), "mean": round(float(vals.mean()), 5),
                            "median": round(float(vals.median()), 5),
                            "frac_positive": round(float((vals > 0).mean()), 4)})
if rows:
    rt = pd.DataFrame(rows)
    dump(rt, "04_jacobian/jacobian_pairs.tsv", index=False)
    print(rt.to_string(), flush=True)
if per_cluster:
    pc = pd.DataFrame(per_cluster)
    dump(pc, "04_jacobian/jacobian_pairs_by_cluster.tsv", index=False)
    print(pc.pivot_table(index='cluster', columns=['regulator', 'effector'],
                         values='mean').round(4).to_string(), flush=True)

stamp("done")
print("total %.1f s" % (time.time() - T0), flush=True)
