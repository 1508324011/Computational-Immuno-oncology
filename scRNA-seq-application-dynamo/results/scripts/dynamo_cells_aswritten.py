#!/usr/bin/env python
# AUTO-EXTRACTED from the tutorial HTML by extract_dynamo_cells.py -- do not hand-edit
# the cell bodies; apply deviations as a documented patch at the top instead.
# @@CELLS 34 code cells

# ===== harness: flush any figure left open by the cell that just ran =====
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as _plt
_FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figure")
os.makedirs(_FIGDIR, exist_ok=True)
_INLINE = []
def _flush(cell, tag):
    for num in _plt.get_fignums():
        fig = _plt.figure(num)
        if not fig.get_axes():
            continue
        name = "inline_cell%02d_%s_%02d.png" % (cell, tag, num)
        fig.savefig(os.path.join(_FIGDIR, name), dpi=200, bbox_inches="tight")
        _INLINE.append(name)
        print("[harness] wrote figure/%s" % name)
    _plt.close("all")


# ===================== cell 0 (patched: `!python --version` -> sys.version) =====================
# Check python version
import sys as _sys
print(_sys.version)
_flush(0, "auto")

# ============ cell 1 (patched: `!ls` teacher path -> our ./data copy) ============
# Check data exists
# Data collected from https://journals.biologists.com/dev/article/146/12/dev173849/19483/Comprehensive-single-cell-mRNA-profiling-reveals-a
print(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/endocrinogenesis_day15.h5ad"))
_flush(1, "auto")

# ===================== cell 2 (verbatim) =====================
# ! source dynamo/bin/activate # in cluster
_flush(2, "auto")

# ===================== cell 3 (verbatim) =====================
# Disable annoying warnings
import warnings
warnings.filterwarnings("ignore")

# Import core packages
import dynamo as dyn
import scanpy as sc

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))   # patched: teacher's hardcoded home dir

# Adjust visualization
# from IPython.core.display import display, HTML
# display(HTML("<style>.container { width:90% !important; }</style>"))

# Direct show figure in jupyter
# %matplotlib inline   (patched: Jupyter magic; harness already forces the Agg backend)
_flush(3, "auto")

# ===================== cell 4 (verbatim) =====================
# check the working dir
os.getcwd()
_flush(4, "auto")

# ===================== cell 5 (verbatim) =====================
# Check dynamo and all dependencies version
dyn.get_all_dependencies_version()
_flush(5, "auto")

# ===================== cell 6 (verbatim) =====================
# Set resolution/size, styling and format of figures
dyn.configuration.set_figure_params("dynamo", background="white", figsize=(3, 2))
_flush(6, "auto")

# ===================== cell 7 (verbatim) =====================
# usually using sc.read_10x_mtx or sc.read_10x_h5
data_path = './data/endocrinogenesis_day15.h5ad'
adata = sc.read_h5ad(data_path)
_flush(7, "auto")

# ===================== cell 8 (verbatim) =====================
adata
_flush(8, "auto")

# ===================== cell 9 (verbatim) =====================
# known differention related genes in pancreas
pancreas_genes = [
    "Hes1", "Nkx6-1", "Nkx2-2",
    "Neurog3", "Neurod1", "Pax4",
    "Pax6", "Arx", "Pdx1",
    "Ins1", "Ins2", "Ghrl",
    "Ptf1a", "Iapp", "Isl1",
    "Sox9", "Gcg",
]

figure_path = './figure/'
_flush(9, "auto")

# ===================== cell 10 (verbatim) =====================
# do preprocessing, analog to scanpy
# select 4000 top variavle genes and filter genes that are detected in less than 20 counts (shared)
# ensure all pancrease genes are in the list
dyn.pp.recipe_monocle(adata, n_top_genes=4000,
                      fg_kwargs={"shared_count": 20},
                      genes_to_append=pancreas_genes)
_flush(10, "auto")

# ===================== cell 11 (verbatim) =====================
# Setting the model into stochastic mode, for cells may not at steady state
dyn.tl.dynamics(adata, model="stochastic")
_flush(11, "auto")

# ===================== cell 12 (verbatim) =====================
# Dimension reduction, and then calculate UMAP basis
dyn.tl.reduceDimension(adata, n_pca_components=30)
_flush(12, "auto")

# ===================== cell 13 (verbatim) =====================
# Project high dimensional velocity vectors onto
# UMAP low dimensional embeddings
dyn.tl.cell_velocities(adata, method="pearson", other_kernels_dict={"transform": "sqrt"})
# PCA low dimensional embeddings
dyn.tl.cell_velocities(adata, basis="pca")
_flush(13, "auto")

# ===================== cell 14 (verbatim) =====================
# Plot the velocity streamline
dyn.pl.streamline_plot(adata, 
                       color=["clusters"], 
                       basis="umap", 
                       show_legend="on data", 
                       show_arrowed_spines=True,
                       save_show_or_return='save',
                       save_kwargs={
                           "prefix": "01_velocity_streamline",  # prefix
                           "ext": "png",                        # file format
                           "path": figure_path,                 # save dir (default current working dir)
                           "dpi": 300,                          # resolution
                           "transparent": True,
                           "close": True,
                           "verbose": True
                       }
                      )
_flush(14, "auto")

# ===================== cell 15 (verbatim) =====================
# Plot the velocity vector per cell
dyn.pl.cell_wise_vectors(
    adata,
    color=["clusters"],
    basis="umap",
    show_legend="on data",
    quiver_length=4,
    quiver_size=4,
    figsize=(5, 3),
    show_arrowed_spines=False,
    save_show_or_return='save',
    save_kwargs={
        "prefix": "02_cell_wise_velocity",  # prefix
        "ext": "png",                       # file format
        "path": figure_path,                # save dir (default current working dir)
        "dpi": 300,                         # resolution
        "transparent": True,
        "close": True,
        "verbose": True
    }
)
_flush(15, "auto")

# ===================== cell 16 (verbatim) =====================
# Setting the progenitor cell types
progenitor = adata.obs_names[adata.obs.clusters.isin(["Ductal"])]
len(progenitor)
_flush(16, "auto")

# ===================== cell 17 (verbatim) =====================
# Estimate the vector field
dyn.vf.VectorField(adata, basis="pca", pot_curl_div=True)
dyn.vf.VectorField(adata, basis="umap", pot_curl_div=True)

# Calculate the speed, divergence, accelaration and curl 
dyn.vf.speed(adata, basis="pca")
dyn.vf.divergence(adata, basis="pca")
dyn.vf.acceleration(adata, basis="pca")
dyn.vf.curl(adata, basis="umap")
_flush(17, "auto")

# ===================== cell 18 (verbatim) =====================
# Draw figures
dyn.configuration.set_figure_params("dynamo", background="black")

import matplotlib.pyplot as plt

fig1, fig1_axes = plt.subplots(ncols=2, nrows=2, constrained_layout=True, figsize=(10, 6))
dyn.pl.cell_wise_vectors(
    adata,
    color="pca_ddhodge_potential",
    pointsize=0.1,
    alpha=0.7,
    ax=fig1_axes[0, 0],
    quiver_length=6,
    quiver_size=6,
    save_show_or_return="return",
    background="black",
)
dyn.pl.grid_vectors(
    adata,
    color="speed_pca",
    basis="umap",
    ax=fig1_axes[0, 1],
    quiver_length=12,
    quiver_size=12,
    save_show_or_return="return",
    background="black",
)
dyn.pl.streamline_plot(
    adata, color="divergence_pca", basis="umap", ax=fig1_axes[1, 0], save_show_or_return="return", background="black"
)
dyn.pl.streamline_plot(
    adata, color="acceleration_pca", basis="umap", ax=fig1_axes[1, 1], save_show_or_return="return", background="black"
)

fig1.savefig('./figure/03_combined_plots.png')
plt.show()
_flush(18, "auto")

# ===================== cell 19 (verbatim) =====================
# Draw the potential cell states
dyn.pl.topography(
    adata, 
    basis="umap", 
    background="white", 
    color=["clusters"],
    streamline_color="black",
    show_legend="on data",
    save_show_or_return='save',
    save_kwargs={
        "prefix": "04_cell_state",  # prefix
        "ext": "png",                       # file format
        "path": figure_path,                       # save dir (default current working dir)
        "dpi": 300,                         # resolution
        "transparent": True,
        "close": True,
        "verbose": True
    }
)
_flush(19, "auto")

# ===================== cell 20 (verbatim) =====================
# Draw the state transition map
dyn.configuration.set_figure_params("dynamo", background="white")
dyn.pd.state_graph(adata, group='clusters', basis='pca', method='vf')
dyn.pl.state_graph(adata,
                   color=['clusters'],
                   group='clusters',
                   basis='umap',
                   show_legend='on data',
                   method='vf',
                   save_show_or_return='save',
                   save_kwargs={
                       "prefix": "05_cell_state_transition_map",  # prefix
                       "ext": "png",                       # file format
                       "path": figure_path,                       # save dir (default current working dir)
                       "dpi": 300,                         # resolution
                       "transparent": True,
                       "close": True,
                       "verbose": True
                   }
                  )
_flush(20, "auto")

# ===================== cell 21 (verbatim) =====================
# Calculate the jacobian on our predefined genes
dyn.vf.jacobian(adata, regulators=pancreas_genes)
_flush(21, "auto")

# ===================== cell 22 (verbatim) =====================
# Ngn3 activates Pax4 in progenitors, initiating pancreatic endocrinogenesis
dyn.pl.jacobian(
    adata,
    basis="umap",
    regulators=[
        "Neurog3",
    ],
    effectors=["Pax4"],
    alpha=1,
    save_show_or_return='save',
    save_kwargs={
        "prefix": "06_reg_Ngn3_activate_Pax4",  # prefix
        "ext": "png",                       # file format
        "path": figure_path,                       # save dir (default current working dir)
        "dpi": 300,                         # resolution
        "transparent": True,
        "close": True,
        "verbose": True
    }
)
_flush(22, "auto")

# ===================== cell 23 (verbatim) =====================
# Jacobian analyses reveal mutual inhibition of Pax4 and Arx at the bifurcation point in progenitors.
dyn.pl.jacobian(adata, basis="umap", regulators=["Arx"], effectors=["Pax4"], alpha=1)
dyn.pl.jacobian(adata, basis="umap", regulators=["Pax4"], effectors=["Arx"], alpha=1)
_flush(23, "auto")

# ===================== cell 24 (verbatim) =====================
# Pdx1 activates Ins2 in beta cells.
dyn.pl.jacobian(adata, basis="umap", regulators=["Pax4"], effectors=["Ins2"], alpha=1)
_flush(24, "auto")

# ===================== cell 25 (verbatim) =====================
adata.obs["clusters"].unique()
_flush(25, "auto")

# ===================== cell 26 (verbatim) =====================
# Rank the elements in the Jacobian
_flush(26, "auto")

# ===================== cell 27 (verbatim) =====================
dyn.vf.acceleration(adata, basis='pca')
dyn.vf.rank_acceleration_genes(adata, groups='clusters')
_flush(27, "auto")

# ===================== cell 28 (verbatim) =====================
adata.uns['rank_acceleration'][:5]
_flush(28, "auto")

# ===================== cell 29 (verbatim) =====================
# Enrichment analysis allows us to see if top ranking genes are significantly enriched in
# certain biological pathways. For example, when we pick the top 100 accelerating genes
# in Ductal cells, we found that they are highly enriched in cell cycle related pathways.
_flush(29, "auto")

# ===================== cell 30 (patched: `!pwd` -> print) =====================
print(os.getcwd())
_flush(30, "auto")

# ===================== cell 31 (verbatim) =====================
enr = dyn.ext.enrichr(
    adata.uns['rank_acceleration']['Ductal'][:100].to_list(),
    organism='mouse',
    gene_sets='GO_Biological_Process_2021',
    outdir='./figure'
)
enr.results.head(3)
_flush(31, "auto")

# ===================== cell 32 (verbatim) =====================
from gseapy.plot import barplot, dotplot
barplot(enr.res2d, column = 'Adjusted P-value',title='GO_Biological_Process_2018', cutoff=0.05)
dotplot(enr.res2d, title='KEGG_2016',cmap='viridis_r', cutoff=0.05)
plt.show()
_flush(32, "auto")

# ===================== cell 33 (verbatim) =====================
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
_flush(33, "auto")

# ===================== end of run =====================
print("[harness] inline figures captured: %d" % len(_INLINE))
for _n in _INLINE:
    print("[harness]   " + _n)
