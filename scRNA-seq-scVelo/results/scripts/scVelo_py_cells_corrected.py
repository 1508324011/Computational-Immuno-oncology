# @@PURPOSE Corrected / annotated edition of the scVelo practical (module: class/scRNA-seq-scVelo).
#            Same data, same 7-step storyline, but every problem logged in
#            scVelo_analysis_corrected.md (VELO-01 ... VELO-16) is addressed, and the dynamical
#            (EM) model -- the model the theory slides actually teach -- is finally run.
# @@USAGE   cd <run dir with data/Pancreas/endocrinogenesis_day15.h5ad> && python <this file>
# @@OUTPUT  figures_corrected/*.png, tables_corrected/*.csv|json, corrected_adata.h5ad
# @@NOTE    The regression run that reproduces the teacher's notebook line-for-line lives in
#            scVelo_py_cells_aswritten.py. This file is deliberately NOT a drop-in replacement
#            for the notebook: corrections are tagged with the VELO ids they close.
# @@TIME    ~10-20 min on one node (recover_dynamics on 3,696 x 2,000 genes is the long pole).

import atexit
import json
import os
import random
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

# ================================= harness (not tutorial material) ======================
# The tutorial never saves a figure, which is issue VELO-11. Instead of monkeypatching every
# plotting call, the harness intercepts matplotlib's show() and writes each rendered figure
# under a stable name. (A production script would use scVelo's own `save=` argument, which
# honours settings.figdir / plot_prefix / file_format_figs; the hook is used here only so that
# the as-written and the corrected run produce identically named panel files.)
FIGNAMES = [
    "corrected_01_proportions.png",
    "corrected_02_velocity_embedding_stream.png",
    "corrected_03_velocity_embedding.png",
    "corrected_04_phaseportrait_markers_stochastic.png",
    "corrected_05_phaseportrait_markers_dynamical.png",
    "corrected_06_cpe_groups_velocity.png",
    "corrected_07_top5_ngn3_high_ep.png",
    "corrected_08_top5_pre_endocrine.png",
    "corrected_09_cellcycle_gradients_original_scores.png",
    "corrected_10_phase_marker_genes.png",
    "corrected_11_phaseportrait_hells_top2a_dynamical.png",
    "corrected_12_latent_time_and_confidence.png",
    "corrected_13_velocity_stream_latent_time.png",
]
FIGDIR = "figures_corrected"
TABLE_DIR = "tables_corrected"
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)
_captured = []
_orig_show = plt.show


def _save_current(name):
    fig = plt.gcf()
    if not fig.get_axes():
        return
    path = os.path.join(FIGDIR, name)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    _captured.append(path)
    print("[harness] wrote %s" % path)
    plt.close("all")


def _show(*args, **kwargs):
    idx = len(_captured)
    _save_current(FIGNAMES[idx] if idx < len(FIGNAMES) else "unlabeled_%02d.png" % idx)
    return None


def _flush():
    if plt.get_fignums():
        idx = len(_captured)
        _save_current(FIGNAMES[idx] if idx < len(FIGNAMES) else "unlabeled_%02d.png" % idx)


plt.show = _show
atexit.register(_flush)
# =================================== end harness =======================================

RNG = 0
random.seed(RNG)
np.random.seed(RNG)

scv.settings.verbosity = 2
scv.settings.presenter_view = True
scv.set_figure_params("scvelo")

T0 = time.time()
SUMMARY = {}


def stamp(msg):
    print("\n=== [%6.1fs] %s" % (time.time() - T0, msg))


def describe(a, tag=""):
    """VELO-06: never print a raw (sparse) matrix -- print what it is.

    `print(adata.X)` on a 3,696 x 27,998 matrix with 9.3 M stored values prints either a flood
    of coordinate triplets or (newer scipy) just a `<sparse matrix>` header; neither tells a
    student anything. Shape / dtype / nnz / density / range do.
    """

    def one(m):
        nnz = m.nnz if hasattr(m, "nnz") else int(np.count_nonzero(m))
        return "%-12s %-8s nnz=%-10d density=%5.1f%%  max=%.3f  sum=%.1f" % (
            "%s" % (getattr(m, "shape", ()),),
            str(m.dtype),
            nnz,
            100.0 * nnz / (m.shape[0] * m.shape[1]),
            float(m.max()),
            float(m.sum()),
        )

    print("[describe] %s X           %s" % (tag, one(a.X)))
    for k in a.layers:
        print("[describe] %s layer %-9s %s" % (tag, k, one(a.layers[k])))


def velocity_gene_count(a):
    """scVelo stores the flag as a Categorical ('True'/'False'/NaN) -- hence the string test."""
    return int((a.var["velocity_genes"].astype(str) == "True").sum())


# ============================== STEP 0 -- load the data ================================
stamp("STEP 0  load endocrinogenesis day 15 (scv.datasets.pancreas, served from the local copy)")
adata = scv.datasets.pancreas()
describe(adata, "raw")
print(adata)
SUMMARY["n_cells_raw"] = int(adata.n_obs)
SUMMARY["n_genes_raw"] = int(adata.n_vars)
SUMMARY["X_nnz_raw"] = int(adata.X.nnz)

# VELO-08: the object *ships* with S_score / G2M_score columns; keep an untouched copy so the
# cell-cycle figure can be drawn from the original annotation instead of a silent overwrite.
ORIG_SCORES = adata.obs[["S_score", "G2M_score"]].copy()

# VELO-09: `adata.var.head(3)` in the tutorial shows `highly_variable_genes`, a leftover
# annotation from the dataset's authors -- it is NOT the HVG set selected below. The real
# selection lands in var['highly_variable'].
print("\n[VELO-09] pre-existing var columns:", list(adata.var.columns))
print("[VELO-09] dtype of the column students actually see:",
      adata.var["highly_variable_genes"].dtype)
print("[VELO-09] its values:",
      adata.var["highly_variable_genes"].astype(str).value_counts(dropna=False).to_dict())

# ============================ STEP 1 -- filtering / normalization ======================
stamp("STEP 1  filter, normalise, log-transform, HVG, moments")
scv.pp.filter_and_normalize(adata, min_shared_counts=20)  # 27,998 -> 7,197 genes
# VELO-01: this log1p is REQUIRED, not a double transform. scv.pp.filter_and_normalize in
# scVelo 0.3.4 only filters + normalises (its docstring still claims a log transform, but the
# implementation never takes the log), so the tutorial's explicit sc.pp.log1p is what makes X
# logarithmic. The layers spliced/unspliced deliberately stay linear: moments and velocities are
# estimated on (normalised) linear counts, only X is logged.
sc.pp.log1p(adata)
describe(adata, "filtered+logged")
SUMMARY["n_genes_filtered"] = int(adata.n_vars)
print("[VELO-02] X is log-normalised while the velocity layers stay linear: "
      "X.max=%.3f  spliced.max=%.1f  unspliced.max=%.1f"
      % (float(adata.X.max()), float(adata.layers["spliced"].max()),
         float(adata.layers["unspliced"].max())))

# VELO-05: n_top_genes=2000 is a default-sized choice, not a result -- STEP 7 tests it.
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", subset=True)
print("[VELO-09] after subsetting: n_vars=%d, var['highly_variable'].sum()=%d"
      % (adata.n_vars, int((adata.var["highly_variable"].astype(str) == "True").sum())))
SUMMARY["n_genes_hvg"] = int(adata.n_vars)

# The PCA shipped with the h5ad was fitted on the full 27,998-gene matrix, so after the HVG
# subset it is stale. The tutorial silently reuses it; STEP 7 quantifies what that costs, while
# the main line stays comparable to the tutorial.
print("[note] reusing the shipped X_pca (fitted on 27,998 genes) for the neighbour graph, as the"
      " tutorial does; the recomputed-PCA variant is measured in STEP 7.")
sc.pp.neighbors(adata, n_neighbors=30, use_rep="X_pca", n_pcs=30, random_state=RNG)
# VELO-07: passing the same n_neighbors the graph already has stops scv.pp.moments from
# rebuilding the neighbours itself -- that rebuild is deprecated since scvelo 0.4.0 and is what
# produces the "Automatic neighbor calculation is deprecated" warning in the as-written run.
scv.pp.moments(adata, n_neighbors=30, n_pcs=30)

# ============================ STEP 2 -- stochastic velocity ===========================
stamp("STEP 2  stochastic (steady-state) velocity + velocity graph -- the tutorial's model")
scv.tl.velocity(adata, mode="stochastic")
# VELO-12: velocity_graph defaults to n_jobs=1 ("using 1/32 cores" in the as-written log).
scv.tl.velocity_graph(adata, n_jobs=-1)
scv.tl.velocity_confidence(adata)
adata.obs["velocity_confidence_stochastic"] = adata.obs["velocity_confidence"].values
SUMMARY["n_velocity_genes_stochastic"] = velocity_gene_count(adata)
print("[state] velocity genes (stochastic): %d / %d"
      % (SUMMARY["n_velocity_genes_stochastic"], adata.n_vars))
# The tutorial ranks velocity genes on the stochastic field; the figure panels below reuse this
# table so that they stay comparable with the reference PNGs.
scv.tl.rank_velocity_genes(adata, groupby="clusters", min_corr=0.3)
STOCH_RANK = pd.DataFrame(adata.uns["rank_velocity_genes"]["names"]).copy()
STOCH_RANK.to_csv(os.path.join(TABLE_DIR, "rank_velocity_genes_top50_stochastic.csv"), index=False)
pd.DataFrame({
    "clusters": STOCH_RANK.columns,
    "top5": [" / ".join(STOCH_RANK[c][:5]) for c in STOCH_RANK.columns],
    "n_ranked": [int(STOCH_RANK[c].notna().sum()) for c in STOCH_RANK.columns],
}).to_csv(os.path.join(TABLE_DIR, "rank_velocity_genes_top5_stochastic.csv"), index=False)

# ============================ STEP 3 -- dynamical (EM) velocity =======================
stamp("STEP 3  dynamical model (the model the theory slides teach) -- VELO-03")
# recover_dynamics fits, per gene, the full splicing kinetics
#   du/dt = alpha(t) - beta*u,   ds/dt = beta*u - gamma*s
# with a switching alpha(t) and a per-gene gamma, by EM over all cells; scv.tl.velocity then
# reports the residual of that fit -- a real velocity, not a steady-state closure.
n_jobs = min(8, os.cpu_count() or 1)
# VELO-17: recover_dynamics() defaults to var_names="velocity_genes", i.e. its gene universe is
# whatever the PREVIOUS velocity run flagged. Running stochastic and then dynamical (the natural
# order, and the one the tutorial implies) therefore fits the EM only on the steady-state model's
# 1,038 genes -- the two models are then not even fitted on the same genes. var_names="all" fits
# every gene that survived the HVG subset, which is what "run the model" should mean.
scv.tl.recover_dynamics(adata, var_names="all", n_jobs=n_jobs, show_progress_bar=False)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata, n_jobs=-1)
scv.tl.velocity_confidence(adata)
scv.tl.latent_time(adata)
dyn_cols = [c for c in ("fit_alpha", "fit_beta", "fit_gamma", "fit_t_", "fit_likelihood",
                        "fit_r2", "fit_steady_states", "fit_scaling") if c in adata.var]
print("[state] dynamical fit columns:", dyn_cols)
SUMMARY["dynamical_fit_columns"] = dyn_cols
SUMMARY["n_velocity_genes_dynamical"] = velocity_gene_count(adata)
if "fit_likelihood" in adata.var:
    SUMMARY["dynamical_median_likelihood"] = float(np.nanmedian(adata.var["fit_likelihood"]))
    SUMMARY["dynamical_frac_likelihood_gt_0"] = float(np.nanmean(adata.var["fit_likelihood"] > 0))
if "fit_gamma" in adata.var:
    g = adata.var["fit_gamma"].to_numpy(dtype=float)
    g = g[np.isfinite(g) & (g > 0)]
    SUMMARY["fit_gamma_n"] = int(g.size)
    SUMMARY["fit_gamma_median"] = float(np.median(g)) if g.size else float("nan")
    SUMMARY["fit_gamma_iqr"] = ([float(np.percentile(g, 25)), float(np.percentile(g, 75))]
                                if g.size else [])
    print("[state] per-gene gamma: n=%d median=%.3f IQR=[%.3f, %.3f]"
          % (g.size, np.median(g), np.percentile(g, 25), np.percentile(g, 75)))
    print("        (the stochastic model forces one shared rate; here it is fitted per gene)")
scv.tl.rank_velocity_genes(adata, groupby="clusters", min_corr=0.3)
DYN_RANK = pd.DataFrame(adata.uns["rank_velocity_genes"]["names"]).copy()
DYN_RANK.to_csv(os.path.join(TABLE_DIR, "rank_velocity_genes_top50_dynamical.csv"), index=False)
pd.DataFrame({
    "clusters": DYN_RANK.columns,
    "top5": [" / ".join(DYN_RANK[c][:5]) for c in DYN_RANK.columns],
    "n_ranked": [int(DYN_RANK[c].notna().sum()) for c in DYN_RANK.columns],
}).to_csv(os.path.join(TABLE_DIR, "rank_velocity_genes_top5_dynamical.csv"), index=False)
jac = []
for c in STOCH_RANK.columns:
    a_set, b_set = set(STOCH_RANK[c][:20].dropna()), set(DYN_RANK[c][:20].dropna())
    jac.append(len(a_set & b_set) / max(1, len(a_set | b_set)))
SUMMARY["top20_jaccard_stochastic_vs_dynamical"] = float(np.mean(jac))
print("[state] stochastic vs dynamical top-20 velocity genes: mean Jaccard = %.3f" % np.mean(jac))

# ============================ STEP 4 -- overview figures =============================
stamp("STEP 4  overview panels")
scv.pl.proportions(adata)
scv.pl.velocity_embedding_stream(adata, basis="umap")
scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120)
# VELO-16: `mode=` is NOT what selects the overlay. In scv.pl.velocity, mode only decides how
# many side panels are drawn (mode="stochastic" adds two variance panels); the lines inside the
# phase portrait come from `fits`, which defaults to ["velocity", "dynamics"] -- i.e. the
# tutorial's figure already shows the dynamical fit next to the steady-state line. Passing
# fits=["dynamics"] drops the steady-state estimate and leaves the fitted kinetics alone, which
# is the panel worth reading when the two models disagree.
scv.pl.velocity(adata, ["Cpe", "Gnao1", "Ins2", "Adk"], ncols=2)                        # both fits
scv.pl.velocity(adata, ["Cpe", "Gnao1", "Ins2", "Adk"], ncols=2, fits=["dynamics"])    # model only
scv.pl.scatter(adata, "Cpe", color=["clusters", "velocity"],
               add_outline="Ngn3 high EP, Pre-endocrine, Beta")

# ============================ STEP 5 -- velocity genes ===============================
stamp("STEP 5  velocity-gene panels (stochastic ranking, i.e. comparable to the reference)")
print(STOCH_RANK.head())
kwargs = dict(frameon=False, size=10, linewidth=1.5,
              add_outline="Ngn3 high EP, Pre-endocrine, Beta")
scv.pl.scatter(adata, STOCH_RANK["Ngn3 high EP"][:5], ylabel="Ngn3 high EP", **kwargs)
scv.pl.scatter(adata, STOCH_RANK["Pre-endocrine"][:5], ylabel="Pre-endocrine", **kwargs)

# ============================ STEP 6 -- cell cycle ==================================
stamp("STEP 6  cell cycle: keep the original scores, rank phase markers without scv.get_df")
# VELO-08: scv.tl.score_genes_cell_cycle overwrites S_score/G2M_score. The overwrite shifts the
# means and repaints the gradient panel, so the original annotation is restored before plotting
# and both versions are written next to each other for comparison.
tmp = adata.copy()
scv.tl.score_genes_cell_cycle(tmp)
recomputed = tmp.obs[["S_score", "G2M_score"]].copy()
cmp_ = pd.DataFrame({
    "S_score_original": ORIG_SCORES["S_score"].values,
    "S_score_recomputed": recomputed["S_score"].values,
    "G2M_score_original": ORIG_SCORES["G2M_score"].values,
    "G2M_score_recomputed": recomputed["G2M_score"].values,
    "clusters": adata.obs["clusters"].values,
}, index=adata.obs_names)
cmp_.to_csv(os.path.join(TABLE_DIR, "cellcycle_scores_original_vs_recomputed.csv"))
SUMMARY["cellcycle_pearson_S"] = float(np.corrcoef(cmp_["S_score_original"], cmp_["S_score_recomputed"])[0, 1])
SUMMARY["cellcycle_pearson_G2M"] = float(np.corrcoef(cmp_["G2M_score_original"], cmp_["G2M_score_recomputed"])[0, 1])
SUMMARY["phase_counts"] = adata.obs["phase"].value_counts().to_dict() if "phase" in adata.obs else {}
print("[VELO-08] original vs recomputed scores: Pearson S=%.4f G2M=%.4f"
      % (SUMMARY["cellcycle_pearson_S"], SUMMARY["cellcycle_pearson_G2M"]))
adata.obs["S_score"] = ORIG_SCORES["S_score"].values
adata.obs["G2M_score"] = ORIG_SCORES["G2M_score"].values
scv.pl.scatter(adata, color_gradients=["S_score", "G2M_score"], smooth=True, perc=[5, 95])

# VELO-07: scv.get_df is deprecated since 0.4.0 -- ranking a var column is one sort_values call.
s_genes, g2m_genes = scv.utils.get_phase_marker_genes(adata)
s_genes = [g for g in s_genes if g in adata.var_names]
g2m_genes = [g for g in g2m_genes if g in adata.var_names]
s_sorted = adata.var.loc[s_genes, "spearmans_score"].sort_values(ascending=False)
g2m_sorted = adata.var.loc[g2m_genes, "spearmans_score"].sort_values(ascending=False)
pd.concat([
    pd.DataFrame({"gene": s_sorted.index, "spearmans_score": s_sorted.values, "phase": "S"}),
    pd.DataFrame({"gene": g2m_sorted.index, "spearmans_score": g2m_sorted.values, "phase": "G2M"}),
]).to_csv(os.path.join(TABLE_DIR, "phase_marker_gene_ranking.csv"), index=False)
print("[state] S top5 = %s | G2M top5 = %s" % (list(s_sorted.index[:5]), list(g2m_sorted.index[:5])))
SUMMARY["cellcycle_S_top5"] = list(map(str, s_sorted.index[:5]))
SUMMARY["cellcycle_G2M_top5"] = list(map(str, g2m_sorted.index[:5]))
kwargs = dict(frameon=False, ylabel="cell cycle genes")
scv.pl.scatter(adata, list(s_sorted.index[:2]) + list(g2m_sorted.index[:3]), **kwargs)
scv.pl.velocity(adata, ["Hells", "Top2a"], ncols=2, fits=["dynamics"], add_outline=True)

# VELO-10: the tutorial's last cell re-draws the very same stream plot as cell 63 (the exported
# PNGs are byte-identical). Replaced by what the stimulus is actually used for: latent time
# (per-cell pseudotime derived from the dynamical fit) and velocity confidence.
scv.pl.scatter(adata, color=["latent_time", "velocity_confidence"], ncols=2, perc=[5, 95],
               color_map="viridis")
scv.pl.velocity_embedding_stream(adata, basis="umap", color="latent_time", color_map="viridis",
                                 title="latent time on the velocity field")

# ============================ STEP 7 -- sensitivity of the magic numbers ==============
stamp("STEP 7  what the unexplained defaults (2000 HVG, k=30, min_corr .3, shipped PCA) cost")
sens_rows = []


def run_variant(tag, n_top_genes=2000, n_neighbors=30, recompute_pca=False):
    a = scv.datasets.pancreas()
    scv.pp.filter_and_normalize(a, min_shared_counts=20)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=n_top_genes, flavor="seurat", subset=True)
    if recompute_pca:
        sc.tl.pca(a, n_comps=30, svd_solver="arpack", random_state=RNG)
    sc.pp.neighbors(a, n_neighbors=n_neighbors, use_rep="X_pca", n_pcs=30, random_state=RNG)
    scv.pp.moments(a, n_neighbors=n_neighbors, n_pcs=30)
    scv.tl.velocity(a, mode="stochastic")
    scv.tl.velocity_graph(a, n_jobs=-1)
    scv.tl.rank_velocity_genes(a, groupby="clusters", min_corr=0.3)
    names = pd.DataFrame(a.uns["rank_velocity_genes"]["names"])
    speed = np.linalg.norm(np.asarray(a.layers["velocity"], dtype=np.float64), axis=1)
    jac = []
    for c in STOCH_RANK.columns:
        a_set, b_set = set(STOCH_RANK[c][:20].dropna()), set(names[c][:20].dropna())
        jac.append(len(a_set & b_set) / max(1, len(a_set | b_set)))
    row = {
        "variant": tag,
        "n_top_genes": n_top_genes,
        "n_neighbors": n_neighbors,
        "recompute_pca": recompute_pca,
        "n_genes": int(a.n_vars),
        "n_velocity_genes": velocity_gene_count(a),
        "speed_median": float(np.median(speed)),
        "n_genes_ranked_median": float(np.median([int(names[c].notna().sum()) for c in names.columns])),
        "top20_jaccard_vs_stochastic_main": float(np.mean(jac)),
        "top3": " | ".join("%s:%s" % (c, " / ".join(names[c][:3])) for c in names.columns),
    }
    sens_rows.append(row)
    print("[sensitivity] %-34s genes=%5d velocity_genes=%4d speed_med=%6.2f top20_J=%.3f"
          % (tag, row["n_genes"], row["n_velocity_genes"], row["speed_median"],
             row["top20_jaccard_vs_stochastic_main"]))
    return a


run_variant("main (2000 HVG, k=30)")
run_variant("1000 HVG, k=30", 1000, 30)
run_variant("3000 HVG, k=30", 3000, 30)
run_variant("2000 HVG, k=15", 2000, 15)
run_variant("2000 HVG, k=50", 2000, 50)
run_variant("2000 HVG, k=30, PCA recomputed", 2000, 30, recompute_pca=True)
pd.DataFrame(sens_rows).to_csv(os.path.join(TABLE_DIR, "sensitivity_graph_and_hvg.csv"), index=False)
SUMMARY["sensitivity"] = sens_rows

# min_corr only filters the ranking (it is compared against adata.var['spearmans_score']), and the
# returned table is capped by n_genes=100 -- so on this dataset every threshold up to 0.9 still
# leaves more than 100 eligible genes and the visible table does not change. The sweep therefore
# records how many genes actually survive each threshold, next to the top-3 of two clusters.
mc_rows = []
for mc in (0.1, 0.3, 0.5, 0.7, 0.9):
    n_pass = int((adata.var["spearmans_score"] > mc).sum())
    scv.tl.rank_velocity_genes(adata, groupby="clusters", min_corr=mc, n_genes=100)
    nm = pd.DataFrame(adata.uns["rank_velocity_genes"]["names"])
    n_filled = [int(nm[c].notna().sum()) for c in nm.columns]
    row = {
        "min_corr": mc,
        "n_genes_passing_threshold": n_pass,
        "table_capped_by_n_genes": 100,
        "n_rows_median": float(np.median(n_filled)),
        "top3_Ductal": " / ".join(nm["Ductal"][:3]),
        "top3_Beta": " / ".join(nm["Beta"][:3]),
    }
    mc_rows.append(row)
    print("[sensitivity] min_corr=%.1f  genes passing=%d  rows in table=%d  Ductal top3=%s"
          % (mc, n_pass, int(np.median(n_filled)), row["top3_Ductal"]))
pd.DataFrame(mc_rows).to_csv(os.path.join(TABLE_DIR, "sensitivity_min_corr.csv"), index=False)
SUMMARY["min_corr_sweep"] = mc_rows

# the n_genes cap is the parameter that actually decides how deep the table goes
ng_rows = []
for ng in (10, 50, 100, 300):
    scv.tl.rank_velocity_genes(adata, groupby="clusters", min_corr=0.3, n_genes=ng)
    nm = pd.DataFrame(adata.uns["rank_velocity_genes"]["names"])
    ng_rows.append({"n_genes": ng, "table_rows": int(nm.shape[0]),
                    "n_rows_per_cluster_median": float(np.median([int(nm[c].notna().sum()) for c in nm.columns]))})
    print("[sensitivity] n_genes=%-4d -> table has %d rows" % (ng, nm.shape[0]))
pd.DataFrame(ng_rows).to_csv(os.path.join(TABLE_DIR, "sensitivity_n_genes_cap.csv"), index=False)
SUMMARY["n_genes_sweep"] = ng_rows
scv.tl.rank_velocity_genes(adata, groupby="clusters", min_corr=0.3)   # restore the main line

# ============================ STEP 8 -- persist =====================================
stamp("STEP 8  persist what the tutorial throws away")
adata.obs.to_csv(os.path.join(TABLE_DIR, "obs_cell_metadata.csv"))
adata.var.to_csv(os.path.join(TABLE_DIR, "var_gene_annotation.csv"))
pd.DataFrame({
    "clusters": adata.obs["clusters"].astype(str),
    "latent_time": adata.obs["latent_time"].values,
    "velocity_confidence_dynamical": adata.obs["velocity_confidence"].values,
    "velocity_confidence_stochastic": adata.obs["velocity_confidence_stochastic"].values,
}).to_csv(os.path.join(TABLE_DIR, "latent_time_and_confidence_per_cell.csv"), index=False)
SUMMARY["runtime_minutes"] = round((time.time() - T0) / 60.0, 2)
SUMMARY["versions"] = {
    "python": sys.version.split()[0],
    "scvelo": scv.__version__,
    "scanpy": sc.__version__,
    "numpy": np.__version__,
    "pandas": pd.__version__,
    "matplotlib": matplotlib.__version__,
}
with open(os.path.join(TABLE_DIR, "corrected_summary.json"), "w") as fh:
    json.dump(SUMMARY, fh, indent=2, default=str)
adata.write("corrected_adata.h5ad")
print("[harness] captured %d figures (planned %d)" % (len(_captured), len(FIGNAMES)))
print("[harness] object written: corrected_adata.h5ad (%s x %s)" % (adata.n_obs, adata.n_vars))
print("\n=== session ===")
for k, v in SUMMARY["versions"].items():
    print("%-12s %s" % (k, v))
if len(_captured) != len(FIGNAMES):
    raise SystemExit("[harness] ERROR: captured %d of %d planned figures"
                     % (len(_captured), len(FIGNAMES)))
print("\n=== done in %.1f min ===" % ((time.time() - T0) / 60.0))
