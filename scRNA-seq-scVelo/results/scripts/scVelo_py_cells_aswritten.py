# @@SOURCE class/scRNA-seq-scVelo/scRNAseq-scVelo.html (Python cells, verbatim)
# @@CELLS  the 29 code cells of the tutorial, in notebook order (Jupyter cell indices 33-85)
# @@RUN    cd /lustre1/user/bjx131_pkuhpc/velo_run_20261104 && python <this file>
# @@NOTE   1. Cell 33 is a Jupyter shell magic (`!pip install scvelo --upgrade --quiet`).
#             It is kept as a comment instead of being executed: the pinned environment
#             (velo: scvelo 0.3.4 == the current PyPI release, i.e. what --upgrade would
#             install anyway) is fixed by conda/pip and must not be mutated mid-run.
# @@NOTE   2. Bare notebook expressions (`adata`, `adata.X`, `df.head()`, ...) are wrapped
#             in print(...) so a plain script reproduces the notebook's displayed output.
#             Nothing else in those cells is changed.
# @@NOTE   3. The tutorial never saves a figure; scv.pl.* in a script prints nothing.
#             The harness below hooks matplotlib's show() to write every rendered figure
#             to ./figures/ under the names used by class/scRNA-seq-scVelo/*.png (the
#             teacher's own output, extracted from the notebook export). No plotting call
#             is modified. The order of FIGNAMES matches the order of the 11 tutorial plots.
# @@FIDEL    figures are written at dpi=200 with bbox_inches="tight", which is exactly how the
#             reference PNGs in class/scRNA-seq-scVelo/ were rendered (pHYs = 200 dpi); at dpi=150
#             every panel would still be correct but 1.33x smaller, which complicates the
#             image-level comparison in results/07_figure_comparison.
# @@CHECK  the run must end with exactly 11 captured figures, else the harness aborts.

# ============================== harness (not tutorial code) ==============================
import atexit
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGNAMES = [
    "01_proportions.png",                      # cell 50  scv.pl.proportions
    "02_velocity_embedding_stream.png",         # cell 63  scv.pl.velocity_embedding_stream
    "03_velocity_embedding.png",                # cell 66  scv.pl.velocity_embedding
    "04_phaseportrait_markers.png",             # cell 69  scv.pl.velocity(Cpe, Gnao1, Ins2, Adk)
    "05_cpe_groups_velocity.png",               # cell 71  scv.pl.scatter(Cpe, clusters+velocity)
    "06_top5_ngn3_high_ep.png",                 # cell 75  scv.pl.scatter(top5 Ngn3 high EP)
    "07_top5_pre_endocrine.png",                # cell 75  scv.pl.scatter(top5 Pre-endocrine)
    "08_cellcycle_gradients.png",               # cell 79  scv.pl.scatter(color_gradients=S/G2M)
    "09_phase_marker_genes.png",                # cell 81  scv.pl.scatter(phase marker genes)
    "10_phaseportrait_hells_top2a.png",         # cell 83  scv.pl.velocity(Hells, Top2a)
    "11_velocity_embedding_stream_repeat.png",  # cell 84  scv.pl.velocity_embedding_stream (again)
]

FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)
_captured = []
_orig_show = plt.show


def _show(*args, **kwargs):
    fig = plt.gcf()
    if not fig.get_axes():          # empty placeholder figure -> nothing to save
        return _orig_show(*args, **kwargs)
    idx = len(_captured)
    name = FIGNAMES[idx] if idx < len(FIGNAMES) else "unlabeled_%02d.png" % idx
    path = os.path.join(FIGDIR, name)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    _captured.append(path)
    print("[harness] wrote %s" % path)
    plt.close("all")
    return None


def _flush_open_figures():
    for num in list(plt.get_fignums()):
        fig = plt.figure(num)
        if fig.get_axes():
            idx = len(_captured)
            name = FIGNAMES[idx] if idx < len(FIGNAMES) else "unlabeled_%02d.png" % idx
            path = os.path.join(FIGDIR, name)
            fig.savefig(path, dpi=200, bbox_inches="tight")
            _captured.append(path)
            print("[harness] wrote %s (from an open figure at exit)" % path)


plt.show = _show
atexit.register(_flush_open_figures)
# ============================ end harness (tutorial code follows) ========================


#@@CELL 33
# update to the latest version, if not done yet.
# !pip install scvelo --upgrade --quiet        # shell magic: not executed (see @@NOTE 1)

#@@CELL 34
import scvelo as scv
import scanpy as sc

#@@CELL 35
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

#@@CELL 36
### Load the Data

#@@CELL 38
adata = scv.datasets.pancreas()
print(adata)

#@@CELL 41
# X for sorting counts matrix, typically in Compressed Sparse Row format
print(adata.X)

#@@CELL 42
# see what exactly does X have
print(adata.X[:5,:10]) # Slice to get the first 5 cells and first 10 genes.

#@@CELL 43
# convert it to a numpy array
print(adata.X.toarray())

#@@CELL 44
print(adata.X.toarray().shape)

#@@CELL 45
# different layer sorted for different gene counts,
# which is used for calculate RNA velocity
# to call this for you data, please refer to velocyto
print(adata.layers['spliced'], adata.layers['unspliced'])

#@@CELL 46
# adata.obs sorted cell meta data, include cluser and some specific score
print(adata.obs.head(3))

#@@CELL 47
# adata.var sorted gene meta data
print(adata.var.head(3))

#@@CELL 48
# adata.obsm sorted some coordinated of cell
print(adata.obsm)

#@@CELL 50
scv.pl.proportions(adata)

#@@CELL 53
scv.pp.filter_and_normalize(
    adata,
    min_shared_counts=20
)
sc.pp.log1p(adata)

#@@CELL 54
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat",
    subset=True
)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

#@@CELL 57
scv.tl.velocity(adata)

#@@CELL 59
scv.tl.velocity_graph(adata)

#@@CELL 63
scv.pl.velocity_embedding_stream(adata, basis='umap')

#@@CELL 66
scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120)

#@@CELL 69
scv.pl.velocity(adata, ['Cpe',  'Gnao1', 'Ins2', 'Adk'], ncols=2)

#@@CELL 71
scv.pl.scatter(adata, 'Cpe', color=['clusters', 'velocity'], 
               add_outline='Ngn3 high EP, Pre-endocrine, Beta')

#@@CELL 74
scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=.3)

import pandas as pd
df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
print(df.head())

#@@CELL 75
kwargs = dict(frameon=False, size=10, linewidth=1.5, 
              add_outline='Ngn3 high EP, Pre-endocrine, Beta')

scv.pl.scatter(adata, df['Ngn3 high EP'][:5], ylabel='Ngn3 high EP', **kwargs)
scv.pl.scatter(adata, df['Pre-endocrine'][:5], ylabel='Pre-endocrine', **kwargs)

#@@CELL 79
scv.tl.score_genes_cell_cycle(adata)
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])

#@@CELL 81
s_genes, g2m_genes = scv.utils.get_phase_marker_genes(adata)
s_genes = scv.get_df(adata[:, s_genes], 'spearmans_score', sort_values=True).index
g2m_genes = scv.get_df(adata[:, g2m_genes], 'spearmans_score', sort_values=True).index

kwargs = dict(frameon=False, ylabel='cell cycle genes')
scv.pl.scatter(adata, list(s_genes[:2]) + list(g2m_genes[:3]), **kwargs)

#@@CELL 83
scv.pl.velocity(adata, ['Hells', 'Top2a'], ncols=2, add_outline=True)

#@@CELL 84
scv.pl.velocity_embedding_stream(adata, basis='umap')

#@@CELL 85
import sys
import pkg_resources

print("### Session Info ###\n")

loaded_pkgs = []

for mod_name, module in sys.modules.items():
    file_path = getattr(module, "__file__", None)
    if (
        file_path
        and "site-packages" in file_path
        and not mod_name.startswith("_")
        and "." not in mod_name
    ):
        try:
            version = pkg_resources.get_distribution(mod_name).version
            loaded_pkgs.append((mod_name, version))
        except Exception:
            pass

for name, version in sorted(set(loaded_pkgs)):
    print(f"{name:<15} v{version}")


# ================================ harness post-check ====================================
# The tutorial persists nothing (no figure, no table, no object). The harness writes the
# final AnnData to the run directory so the corrected pass / audits can be checked against
# exactly the object this run produced. It is NOT part of the tutorial and never enters git.
adata.write("aswritten_adata.h5ad")
print("[harness] wrote aswritten_adata.h5ad (%s x %s)" % (adata.n_obs, adata.n_vars))

print("\n[harness] captured %d figures (expected %d)" % (len(_captured), len(FIGNAMES)))
if len(_captured) != len(FIGNAMES):
    raise SystemExit("[harness] ERROR: figure count mismatch -- see figures/ for what was produced")
# ============================== end harness post-check ==================================
