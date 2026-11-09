#!/usr/bin/env python
"""dynamo_vf_percell_export.py -- per-cell vector-field read-outs plus the cluster-level tests
the corrected edition needs when it checks the tutorial's own reading of the panels.

The tutorial's prose (cell 20 and the markdown after it) makes four claims that are checkable
against the numbers it never prints:

  "a saddle point in ductal cells (negative divergence)"
  "exit from this state to early endocrine progenitors (high positive acceleration)"
  "the bifurcation point for late progenitors ... (high acceleration and positive divergence)"
  "stable cell types (negative divergence)"

So: dump speed / divergence / acceleration / curl / ddhodge potential per cell, then run
Kruskal-Wallis across the eight clusters and a per-cluster sign summary.

usage: DYNAMO_RUN=~/dynamo_run_20261109 python dynamo_vf_percell_export.py <module_results_dir>
"""
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
import dynamo as dyn
import scanpy as sc

OUT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "results")
RUN = os.path.abspath(os.environ.get("DYNAMO_RUN", os.path.expanduser("~/dynamo_run_20261109")))
os.chdir(RUN)
T0 = time.time()

PANCREAS = ["Hes1", "Nkx6-1", "Nkx2-2", "Neurog3", "Neurod1", "Pax4", "Pax6", "Arx", "Pdx1",
            "Ins1", "Ins2", "Ghrl", "Ptf1a", "Iapp", "Isl1", "Sox9", "Gcg"]
QUANTS = ["speed_pca", "divergence_pca", "divergence_umap", "acceleration_pca", "curl_umap",
          "pca_ddhodge_potential"]


def stamp(msg):
    print("\n[%6.1fs] === %s ===" % (time.time() - T0, msg), flush=True)


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

stamp("per-cell table")
cols = [c for c in QUANTS if c in adata.obs]
tab = adata.obs[["clusters"] + cols].copy()
if 'velocity_umap' in adata.obsm:
    tab['velocity_umap_norm'] = np.linalg.norm(np.asarray(adata.obsm['velocity_umap'], float), axis=1)
tab.index.name = 'cell'
p = os.path.join(OUT, "03_vectorfield", "vf_per_cell.tsv")
tab.to_csv(p, sep="\t")
print("[export] 03_vectorfield/vf_per_cell.tsv  (%s)" % (tab.shape,), flush=True)

stamp("cluster-level tests")
rows = []
for c in cols:
    groups = [g[c].dropna().values for _, g in tab.groupby('clusters', observed=True)]
    if len(groups) < 2:
        continue
    H, p_kw = stats.kruskal(*groups)
    ductal = tab.loc[tab.clusters == 'Ductal', c].dropna()
    rest = tab.loc[tab.clusters != 'Ductal', c].dropna()
    U, p_mw = stats.mannwhitneyu(ductal, rest, alternative='two-sided')
    rows.append({
        "quantity": c,
        "kruskal_H": round(float(H), 3),
        "kruskal_p": float("%.3e" % p_kw),
        "ductal_mean": round(float(ductal.mean()), 5),
        "rest_mean": round(float(rest.mean()), 5),
        "mannwhitney_p_ductal_vs_rest": float("%.3e" % p_mw),
        "n_ductal": int(ductal.size), "n_rest": int(rest.size),
        "frac_negative_ductal": round(float((ductal < 0).mean()), 4),
        "frac_negative_rest": round(float((rest < 0).mean()), 4),
    })
tests = pd.DataFrame(rows)
p = os.path.join(OUT, "03_vectorfield", "vf_cluster_tests.tsv")
tests.to_csv(p, sep="\t", index=False)
print("[export] 03_vectorfield/vf_cluster_tests.tsv  (%s)" % (tests.shape,), flush=True)
print(tests.to_string(), flush=True)

# per-cluster sign summary: is divergence negative everywhere (the tutorial's "stable" claim)?
rows = []
for name, g in tab.groupby('clusters', observed=True):
    r = {"cluster": name, "n_cells": int(len(g))}
    for c in ["divergence_pca", "acceleration_pca", "curl_umap"]:
        r[c + "_mean"] = round(float(g[c].mean()), 5)
        r[c + "_frac_pos"] = round(float((g[c] > 0).mean()), 4)
    rows.append(r)
sign = pd.DataFrame(rows)
p = os.path.join(OUT, "03_vectorfield/vf_cluster_sign_summary.tsv")
sign.to_csv(p, sep="\t", index=False)
print("[export] 03_vectorfield/vf_cluster_sign_summary.tsv  (%s)" % (sign.shape,), flush=True)
print(sign.to_string(), flush=True)

stamp("done")
print("total %.1f s" % (time.time() - T0), flush=True)
