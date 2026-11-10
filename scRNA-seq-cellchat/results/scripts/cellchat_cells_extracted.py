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


# ===================== cell 0 (verbatim) =====================
### better run in a conda enviroments
! pip install cellphonedb
_flush(0, "auto")

# ===================== cell 1 (verbatim) =====================
import pandas as pd
import glob
import os

from IPython.display import HTML, display
from cellphonedb.utils import db_releases_utils
_flush(1, "auto")

# ===================== cell 2 (verbatim) =====================
# Direct show figure in jupyter
%matplotlib inline
_flush(2, "auto")

# ===================== cell 3 (verbatim) =====================
# Check database version
display(HTML(db_releases_utils.get_remote_database_versions_html()['db_releases_html_table']))
_flush(3, "auto")

# ===================== cell 4 (verbatim) =====================
# -- Version of the databse
cpdb_version = 'v5.0.0'

# -- Path where the input files to generate the database are located
cpdb_target_dir = os.path.join('./data/cellphoneDB', cpdb_version)
_flush(4, "auto")

# ===================== cell 5 (verbatim) =====================
# donwload the database
# from cellphonedb.utils import db_utils

# db_utils.download_database(cpdb_target_dir, cpdb_version)
_flush(5, "auto")

# ===================== cell 6 (verbatim) =====================
# set paths
cpdb_file_path = 'data/cellphoneDB/v5.0.0/cellphonedb.zip'
meta_file_path = 'data/cellphoneDB/metadata.tsv'
counts_file_path = 'data/cellphoneDB/normalised_log_counts.h5ad'
microenvs_file_path = 'data/cellphoneDB/microenvironment.tsv'
out_path = 'results/cellphoneDB/score/'
_flush(6, "auto")

# ===================== cell 7 (verbatim) =====================
metadata = pd.read_csv(meta_file_path, sep = '\t')
metadata.head(3)
_flush(7, "auto")

# ===================== cell 8 (verbatim) =====================
import anndata

adata = anndata.read_h5ad(counts_file_path)
adata.shape
_flush(8, "auto")

# ===================== cell 9 (verbatim) =====================
list(adata.obs.index).sort() == list(metadata['barcode_sample']).sort()
_flush(9, "auto")

# ===================== cell 10 (verbatim) =====================
microenv = pd.read_csv(microenvs_file_path,
                       sep = '\t')
microenv.head(3)
_flush(10, "auto")

# ===================== cell 11 (verbatim) =====================
microenv.groupby('microenvironment', group_keys = False)['cell_type'] \
    .apply(lambda x : list(x.value_counts().index))
_flush(11, "auto")

# ===================== cell 12 (verbatim) =====================
from cellphonedb.src.core.methods import cpdb_analysis_method

cpdb_results = cpdb_analysis_method.call(
    cpdb_file_path = cpdb_file_path,           # mandatory: CellphoneDB database zip file.
    meta_file_path = meta_file_path,           # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,       # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    counts_data = 'hgnc_symbol',               # defines the gene annotation in counts matrix.
    microenvs_file_path = microenvs_file_path, # optional (default: None): defines cells per microenvironment.
    score_interactions = True,                 # optional: whether to score interactions or not. 
    output_path = out_path,                    # Path to save results    microenvs_file_path = None,
    separator = '|',                           # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    threads = 5,                               # number of threads to use in the analysis.
    threshold = 0.1,                           # defines the min % of cells expressing a gene for this to be employed in the analysis.
    result_precision = 3,                      # Sets the rounding for the mean values in significan_means.
    debug = False,                             # Saves all intermediate tables emplyed during the analysis in pkl format.
    output_suffix = None                       # Replaces the timestamp in the output files by a user defined string in the  (default: None)
)
_flush(12, "auto")

# ===================== cell 13 (verbatim) =====================
print(cpdb_results.keys())
_flush(13, "auto")

# ===================== cell 14 (verbatim) =====================
cpdb_results['means_result'].head(2)
_flush(14, "auto")

# ===================== cell 15 (verbatim) =====================
cpdb_results['interaction_scores'].head(2)
_flush(15, "auto")

# ===================== cell 16 (verbatim) =====================
cpdb_results['deconvoluted'].head(2)
_flush(16, "auto")

# ===================== cell 17 (verbatim) =====================
cpdb_results['deconvoluted_percents'].head(2)
_flush(17, "auto")

# ===================== cell 18 (verbatim) =====================
# set path
cpdb_file_path = 'data/cellphoneDB/v5.0.0/cellphonedb.zip'
meta_file_path = 'data/cellphoneDB/metadata.tsv'
counts_file_path = 'data/cellphoneDB/normalised_log_counts.h5ad'
microenvs_file_path = 'data/cellphoneDB/microenvironment.tsv'
active_tf_path = 'data/cellphoneDB/active_TFs.tsv'
out_path = 'results/cellphoneDB/sig/'
_flush(18, "auto")

# ===================== cell 19 (verbatim) =====================
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                 # mandatory: CellphoneDB database zip file.
    meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,             # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
    active_tfs_file_path = active_tf_path,           # optional: defines cell types and their active TFs.
    microenvs_file_path = microenvs_file_path,       # optional (default: None): defines cells per microenvironment.
    score_interactions = True,                       # optional: whether to score interactions or not. 
    iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
    threshold = 0.1,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = 5,                                     # number of threads to use in the analysis.
    debug_seed = 42,                                 # debug randome seed. To disable >=0.
    result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
    pvalue = 0.05,                                   # P-value threshold to employ for significance.
    subsampling = False,                             # To enable subsampling the data (geometri sketching).
    subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
    subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
    output_path = out_path,                          # Path to save results.
    output_suffix = None                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    )
_flush(19, "auto")

# ===================== cell 20 (verbatim) =====================
# help(cpdb_statistical_analysis_method.call)
_flush(20, "auto")

# ===================== cell 21 (verbatim) =====================
list(cpdb_results.keys())
_flush(21, "auto")

# ===================== cell 22 (verbatim) =====================
cpdb_results['pvalues'].head(2)
_flush(22, "auto")

# ===================== cell 23 (verbatim) =====================
cpdb_results['means'].head(2)
_flush(23, "auto")

# ===================== cell 24 (verbatim) =====================
cpdb_results['significant_means'].head(2)
_flush(24, "auto")

# ===================== cell 25 (verbatim) =====================
cpdb_results['interaction_scores'].head(2)
_flush(25, "auto")

# ===================== cell 26 (verbatim) =====================
cpdb_results['CellSign_active_interactions'].head(2)
_flush(26, "auto")

# ===================== cell 27 (verbatim) =====================
cpdb_results['deconvoluted'].head(2)
_flush(27, "auto")

# ===================== cell 28 (verbatim) =====================
from cellphonedb.utils import search_utils

search_results = search_utils.search_analysis_results(
    query_cell_types_1 = ['EVT_1', 'EVT_2', 'GC', 'eEVT', 'iEVT'],  # List of cells 1, will be paired to cells 2 (list or 'All').
    query_cell_types_2 = ['PV MMP11', 'PV MYH11', 'PV STEAP4'],     # List of cells 2, will be paired to cells 1 (list or 'All').
    query_genes = ['TGFBR1'],                                       # filter interactions based on the genes participating (list).
    query_interactions = ['CSF1_CSF1R'],                            # filter intereactions based on their name (list).
    significant_means = cpdb_results['significant_means'],          # significant_means file generated by CellphoneDB.
    deconvoluted = cpdb_results['deconvoluted'],                    # devonvoluted file generated by CellphoneDB.
    interaction_scores = cpdb_results['interaction_scores'],        # interaction score generated by CellphoneDB.
    query_minimum_score = 50,                                       # minimum score that an interaction must have to be filtered.
    separator = '|',                                                # separator (default: |) employed to split cells (cellA|cellB).
    long_format = True,                                             # converts the output into a wide table, removing non-significant interactions
    query_classifications = ['Signaling by Transforming growth factor']
)

search_results.head()
_flush(28, "auto")

# ===================== cell 29 (verbatim) =====================
import os
import anndata as ad
import pandas as pd
import ktplotspy as kpy
import matplotlib.pyplot as plt
%matplotlib inline
_flush(29, "auto")

# ===================== cell 30 (verbatim) =====================
kpy.plot_cpdb_heatmap(pvals = cpdb_results['pvalues'],
                      degs_analysis = False,
                      figsize = (5, 5),
                      title = "Sum of significant interactions")
_flush(30, "auto")

# ===================== cell 31 (verbatim) =====================
kpy.plot_cpdb(
    adata = adata,
    cell_type1 = "PV MYH11|PV STEAP4|PV MMPP11",
    cell_type2 = "EVT_1|EVT_2|GC|iEVT|eEVT|VCT_CCC",
    means = cpdb_results['means'],
    pvals = cpdb_results['pvalues'],
    celltype_key = "cell_labels",
    genes = ["TGFB2", "CSF1R"],
    figsize = (10, 3),
    title = "Interactions between\nPV and trophoblast",
    max_size = 3,
    highlight_size = 0.75,
    degs_analysis = False,
    standard_scale = True,
    interaction_scores = cpdb_results['interaction_scores'],
    scale_alpha_by_interaction_scores = True
)
_flush(31, "auto")

# ===================== cell 32 (verbatim) =====================
from plotnine import facet_wrap

p = kpy.plot_cpdb(
    adata = adata,
    cell_type1 = "PV MYH11",
    cell_type2 = "EVT_1|EVT_2|GC|iEVT|eEVT|VCT_CCC",
    means = cpdb_results['means'],
    pvals = cpdb_results['pvalues'],
    celltype_key = "cell_labels",
    genes = ["TGFB2", "CSF1R", "COL1A1"],
    figsize = (12, 8),
    title = "Interactions between PV and trophoblast\ns grouped by classification",
    max_size = 6,
    highlight_size = 0.75,
    degs_analysis = False,
    standard_scale = True,
)
p + facet_wrap("~ classification", ncol = 1)
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
