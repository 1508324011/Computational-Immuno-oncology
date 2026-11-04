# @@PURPOSE Dump the quantitative state of a finished scVelo run into plain CSVs/TSVs so the
#            numbers in the corrected edition (and the per-figure comparison) can be traced
#            back to a file instead of a console screenshot.
# @@USAGE   python scVelo_audit_run.py <adata.h5ad> <outdir> [--label aswritten|corrected]
# @@NOTE    Read-only with respect to the object; every output is a small text table.
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def dense(x):
    return x.toarray() if sp.issparse(x) else np.asarray(x)


def bool_sum(series):
    """Count True values in a boolean OR categorical ('True'/'False'/'nan') column.

    scVelo / the course h5ad store `highly_variable_genes` as a pandas Categorical, on
    which Series.sum() raises TypeError -- hence this helper.
    """
    return int((series.astype(str) == "True").sum())


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__ or "usage: script <adata.h5ad> <outdir>")
    h5ad, outdir = sys.argv[1], sys.argv[2]
    label = sys.argv[4] if len(sys.argv) > 4 and sys.argv[3] == "--label" else os.path.basename(outdir)
    print("[audit] %s -> %s (label=%s)" % (h5ad, outdir, label))

    a = ad.read_h5ad(h5ad)
    for sub in ("01_preprocess", "02_velocity", "04_gene_dynamics", "05_cellcycle"):
        os.makedirs(os.path.join(outdir, sub), exist_ok=True)

    # ---------------------------------------------------------------- 01 preprocess
    a.obs.to_csv(os.path.join(outdir, "01_preprocess", "obs_cell_metadata.csv"))
    a.var.to_csv(os.path.join(outdir, "01_preprocess", "var_gene_annotation.csv"))

    summ = []

    def kv(key, value):
        summ.append((key, value))

    kv("label", label)
    kv("n_cells", a.n_obs)
    kv("n_genes", a.n_vars)
    kv("X_nnz", a.X.nnz if sp.issparse(a.X) else int(np.count_nonzero(a.X)))
    kv("X_max", float(a.X.max()))
    kv("X_sum", float(a.X.sum()))
    for layer in ("spliced", "unspliced"):
        if layer in a.layers:
            m = a.layers[layer]
            kv("%s_nnz" % layer, m.nnz if sp.issparse(m) else int(np.count_nonzero(m)))
            kv("%s_max" % layer, float(dense(m).max()))
            kv("%s_sum" % layer, float(m.sum()))
    if "Ms" in a.layers:
        kv("Ms_sum", float(dense(a.layers["Ms"]).sum()))
        kv("Mu_sum", float(dense(a.layers["Mu"]).sum()))
    if "highly_variable" in a.var:
        kv("n_highly_variable", bool_sum(a.var["highly_variable"]))
    if "highly_variable_genes" in a.var:
        kv("n_highly_variable_genes_flag", bool_sum(a.var["highly_variable_genes"]))
    if "velocity_genes" in a.var:
        kv("n_velocity_genes", bool_sum(a.var["velocity_genes"]))
    kv("n_clusters", a.obs["clusters"].nunique())
    kv("clusters", "|".join(map(str, a.obs["clusters"].cat.categories)))
    if "velocity_params" in a.uns:
        for k, v in dict(a.uns["velocity_params"]).items():
            kv("velocity_params.%s" % k, v)
    pd.DataFrame(summ, columns=["key", "value"]).to_csv(
        os.path.join(outdir, "01_preprocess", "run_summary_%s.csv" % label), index=False
    )

    cnt = a.obs["clusters"].value_counts()
    pd.DataFrame({"clusters": cnt.index, "n_cells": cnt.values}).to_csv(
        os.path.join(outdir, "01_preprocess", "cluster_cell_counts.csv"), index=False
    )

    # What scv.pl.proportions visualises: per-cluster share of spliced/unspliced counts.
    # IMPORTANT (VELO-14): proportions() defaults to use_raw=True, i.e. it plots the *initial*
    # library sizes captured before filtering/normalisation (obs['initial_size_spliced'],
    # obs['initial_size_unspliced']), NOT the current layers. On this data set the two
    # definitions differ by 3-4 percentage points per cluster, so both are written out.
    spliced = np.asarray(a.layers["spliced"].sum(axis=1)).ravel()
    unspliced = np.asarray(a.layers["unspliced"].sum(axis=1)).ravel()
    prop = pd.DataFrame({"clusters": a.obs["clusters"].values, "spliced": spliced, "unspliced": unspliced})
    if "initial_size_spliced" in a.obs and "initial_size_unspliced" in a.obs:
        prop["spliced_initial"] = a.obs["initial_size_spliced"].values
        prop["unspliced_initial"] = a.obs["initial_size_unspliced"].values
    agg = prop.groupby("clusters", observed=True).agg(
        n_cells=("spliced", "size"),
        spliced_sum=("spliced", "sum"),
        unspliced_sum=("unspliced", "sum"),
        spliced_mean=("spliced", "mean"),
        unspliced_mean=("unspliced", "mean"),
    )
    agg["unspliced_fraction_current_layers"] = agg["unspliced_sum"] / (
        agg["spliced_sum"] + agg["unspliced_sum"])
    if "spliced_initial" in prop:
        ini = prop.groupby("clusters", observed=True)[["spliced_initial", "unspliced_initial"]].sum()
        ini["unspliced_fraction_initial_size"] = ini["unspliced_initial"] / (
            ini["spliced_initial"] + ini["unspliced_initial"])
        agg = agg.join(ini)
        tot_i = ini[["spliced_initial", "unspliced_initial"]].sum()
        kv("unspliced_fraction_all_cells_initial_size",
           float(tot_i["unspliced_initial"] / tot_i.sum()))
        kv("unspliced_fraction_all_cells_current_layers",
           float(unspliced.sum() / (unspliced.sum() + spliced.sum())))
    agg.to_csv(os.path.join(outdir, "01_preprocess", "proportions_by_cluster.csv"))

    # ---------------------------------------------------------------- 02 velocity
    if "velocity" in a.layers:
        vel = dense(a.layers["velocity"])
    elif "velocity" in a.obsm:
        vel = np.asarray(a.obsm["velocity"])
    else:
        vel = None
    if vel is not None:
        # per-cell speed = L2 norm of the velocity vector, as used by scv.tl.velocity_confidence
        speed = np.linalg.norm(vel, axis=1)
        per_cell = pd.DataFrame(
            {
                "cell": a.obs_names,
                "clusters": a.obs["clusters"].values,
                "velocity_norm": speed,
                "velocity_mean": vel.mean(axis=1),
                "velocity_absmean": np.abs(vel).mean(axis=1),
            }
        )
        per_cell.to_csv(os.path.join(outdir, "02_velocity", "velocity_per_cell.csv"), index=False)
        by = per_cell.groupby("clusters", observed=True).agg(
            n_cells=("velocity_norm", "size"),
            speed_mean=("velocity_norm", "mean"),
            speed_median=("velocity_norm", "median"),
            speed_q05=("velocity_norm", lambda s: s.quantile(0.05)),
            speed_q95=("velocity_norm", lambda s: s.quantile(0.95)),
            velocity_mean=("velocity_mean", "mean"),
        )
        # genes with the largest absolute velocity contribution per cluster
        rows = []
        for cl, idx in a.obs.groupby("clusters", observed=True).indices.items():
            sub = vel[idx]
            mean_v = sub.mean(axis=0)
            order = np.argsort(-np.abs(mean_v))[:20]
            for rank, g in enumerate(order, start=1):
                rows.append(
                    {
                        "clusters": cl,
                        "rank": rank,
                        "gene": a.var_names[g],
                        "mean_velocity": mean_v[g],
                        "n_cells_positive": int((sub[:, g] > 0).sum()),
                    }
                )
        pd.DataFrame(rows).to_csv(
            os.path.join(outdir, "02_velocity", "top_abs_velocity_genes_by_cluster.csv"), index=False
        )
        by.to_csv(os.path.join(outdir, "02_velocity", "velocity_summary_by_cluster.csv"))

    if "latent_time" in a.obs:
        lat = a.obs[["latent_time", "clusters"]].copy()
        if "velocity_confidence" in a.obs:
            lat["velocity_confidence"] = a.obs["velocity_confidence"].values
        lat.to_csv(os.path.join(outdir, "02_velocity", "latent_time_per_cell.csv"))
        by = a.obs.groupby("clusters", observed=True).agg(
            n_cells=("latent_time", "size"),
            latent_time_mean=("latent_time", "mean"),
            latent_time_median=("latent_time", "median"),
            latent_time_min=("latent_time", "min"),
            latent_time_max=("latent_time", "max"),
        )
        if "velocity_confidence" in a.obs:
            by["velocity_confidence_median"] = a.obs.groupby(
                "clusters", observed=True)["velocity_confidence"].median()
        by.to_csv(os.path.join(outdir, "02_velocity", "latent_time_summary_by_cluster.csv"))

    if "velocity_graph" in a.uns:
        g = a.uns["velocity_graph"]
        if sp.issparse(g):
            per_cell_graph = np.asarray((g != 0).sum(axis=1)).ravel()
            rows = [
                ("graph_nnz", int(g.nnz)),
                ("graph_shape", "%dx%d" % g.shape),
                ("nnz_per_cell_mean", float(per_cell_graph.mean())),
                ("nnz_per_cell_median", float(np.median(per_cell_graph))),
                ("graph_max", float(g.max())),
                ("graph_min", float(g.min())),
                ("graph_density", float(g.nnz / (g.shape[0] ** 2))),
            ]
            pd.DataFrame(rows, columns=["key", "value"]).to_csv(
                os.path.join(outdir, "02_velocity", "velocity_graph_%s.csv" % label), index=False
            )

    if "velocity_umap" in a.obsm:
        vu = np.asarray(a.obsm["velocity_umap"])
        pd.DataFrame(
            {
                "cell": a.obs_names,
                "clusters": a.obs["clusters"].values,
                "umap_x": a.obsm["X_umap"][:, 0],
                "umap_y": a.obsm["X_umap"][:, 1],
                "velocity_umap_x": vu[:, 0],
                "velocity_umap_y": vu[:, 1],
            }
        ).to_csv(os.path.join(outdir, "02_velocity", "embedding_velocity_umap.csv"), index=False)

    # ---------------------------------------------------------------- 04 gene dynamics
    if "rank_velocity_genes" in a.uns:
        names = pd.DataFrame(a.uns["rank_velocity_genes"]["names"])
        names.to_csv(os.path.join(outdir, "04_gene_dynamics", "rank_velocity_genes_top50.csv"), index=False)
        pd.DataFrame(
            {
                "clusters": names.columns,
                "top5": [" / ".join(names[c][:5]) for c in names.columns],
                "top10": [" / ".join(names[c][:10]) for c in names.columns],
                "n_genes_ranked": [int(names[c].notna().sum()) for c in names.columns],
            }
        ).to_csv(os.path.join(outdir, "04_gene_dynamics", "rank_velocity_genes_top5_summary.csv"), index=False)

    score_cols = [c for c in ("spearmans_score", "velocity_score", "velocity_r2", "velocity_gamma", "velocity_qreg_ratio") if c in a.var]
    if score_cols:
        tab = a.var[score_cols].copy()
        tab.index.name = "gene"
        if "spearmans_score" in tab:
            tab = tab.sort_values("spearmans_score", ascending=False)
        tab.to_csv(os.path.join(outdir, "04_gene_dynamics", "gene_kinetics_scores.csv"))

    focus = ["Cpe", "Gnao1", "Ins2", "Adk", "Hells", "Top2a"]
    focus = [g for g in focus if g in a.var_names]
    if focus and "velocity" in a.layers:
        rows = []
        for g in focus:
            j = a.var_names.get_loc(g)
            s = np.asarray(a.layers["spliced"][:, j].todense()).ravel() if sp.issparse(a.layers["spliced"]) else a.layers["spliced"][:, j]
            u = np.asarray(a.layers["unspliced"][:, j].todense()).ravel() if sp.issparse(a.layers["unspliced"]) else a.layers["unspliced"][:, j]
            v = vel[:, j]
            rows.append(
                {
                    "gene": g,
                    "spliced_mean": float(np.mean(s)),
                    "unspliced_mean": float(np.mean(u)),
                    "spliced_sum": float(np.sum(s)),
                    "unspliced_sum": float(np.sum(u)),
                    "velocity_mean": float(np.mean(v)),
                    "velocity_absmean": float(np.mean(np.abs(v))),
                    "frac_cells_velocity_pos": float(np.mean(v > 0)),
                    "velocity_gamma": float(a.var.loc[g, "velocity_gamma"]) if "velocity_gamma" in a.var else np.nan,
                    "velocity_r2": float(a.var.loc[g, "velocity_r2"]) if "velocity_r2" in a.var else np.nan,
                    "spearmans_score": float(a.var.loc[g, "spearmans_score"]) if "spearmans_score" in a.var else np.nan,
                }
            )
        pd.DataFrame(rows).to_csv(
            os.path.join(outdir, "04_gene_dynamics", "phase_portrait_gene_stats.csv"), index=False
        )

    # ---------------------------------------------------------------- 05 cell cycle
    if {"S_score", "G2M_score"} <= set(a.obs.columns):
        cyc = a.obs[["S_score", "G2M_score"]].copy()
        cyc["clusters"] = a.obs["clusters"].values
        if "phase" in a.obs:
            cyc["phase"] = a.obs["phase"].values
        cyc.to_csv(os.path.join(outdir, "05_cellcycle", "cellcycle_scores_per_cell.csv"))
        grp = cyc.groupby("clusters", observed=True).agg(
            n_cells=("S_score", "size"),
            S_score_mean=("S_score", "mean"),
            G2M_score_mean=("G2M_score", "mean"),
        )
        if "phase" in cyc:
            phases = pd.crosstab(cyc["clusters"], cyc["phase"])
            phases = phases.div(phases.sum(axis=1), axis=0)
            for col in phases.columns:
                grp["frac_%s" % col] = phases[col]
        grp.to_csv(os.path.join(outdir, "05_cellcycle", "cellcycle_summary_by_cluster.csv"))

    # phase marker genes: scv's own S/G2M lists ranked by spearmans_score
    try:
        import scvelo as scv

        s_genes, g2m_genes = scv.utils.get_phase_marker_genes(a)
        s_genes = [g for g in s_genes if g in a.var_names]
        g2m_genes = [g for g in g2m_genes if g in a.var_names]
        s_df = scv.get_df(a[:, s_genes], "spearmans_score", sort_values=True)
        g_df = scv.get_df(a[:, g2m_genes], "spearmans_score", sort_values=True)
        out = pd.concat(
            [
                pd.DataFrame({"gene": s_df.index, "spearmans_score": s_df["spearmans_score"].values, "phase": "S"}),
                pd.DataFrame({"gene": g_df.index, "spearmans_score": g_df["spearmans_score"].values, "phase": "G2M"}),
            ]
        )
        out.to_csv(os.path.join(outdir, "05_cellcycle", "phase_marker_gene_ranking.csv"), index=False)
        pd.DataFrame(
            {
                "phase": ["S", "G2M"],
                "n_marker_genes_in_object": [len(s_genes), len(g2m_genes)],
                "top5": [" / ".join(s_df.index[:5]), " / ".join(g_df.index[:5])],
            }
        ).to_csv(os.path.join(outdir, "05_cellcycle", "phase_marker_top5.csv"), index=False)
        print("[audit] phase markers: S top5 =", list(s_df.index[:5]), "| G2M top5 =", list(g_df.index[:5]))
    except Exception as exc:  # pragma: no cover - audit must not fail the run
        print("[audit] phase marker ranking skipped: %r" % (exc,))

    print("[audit] done ->", outdir)


if __name__ == "__main__":
    main()
