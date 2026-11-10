#!/usr/bin/env python3
"""cellchat_corrected.py -- the corrected/enhanced pass for the CellphoneDB module (lesson 8).

Baseline: cellchat_cells_aswritten.py, i.e. the tutorial's own 34 cells.  Nothing here replaces that
run; every fix is an *additional* computation whose delta against the baseline is reported.  Steps:

  C1  the misspelt cell type  "PV MMPP11" -> "PV MMP11" and the PV-vs-trophoblast panel re-drawn
      (the tutorial's cell 31 silently drops a third of its columns, see CPB-07)
  C2  reproducible output names (output_suffix) instead of the tutorial's timestamps (CPB-05)
  C3  persist the tables the tutorial only *displays* (cells 7-28) -- CPB-04
  C4  replace the vacuous barcode check of cell 9 (`list(...).sort() == list(...).sort()`, which is
      `None == None`) with a real set comparison (CPB-08)
  C5  threshold sensitivity: threshold 0.1 -> 0.05 and 0.20, both analyses (CPB-09)
  C6  degs_analysis=True with the shipped DEGs_inv_trophoblast.tsv, i.e. wake the dead input
      (CPB-10)
  C7  statistical analysis WITHOUT active_TFs.tsv, to size what CellSign actually adds (CPB-11)
  C8  iterations 1000 -> 100, to see whether the p-values that survive FDR-less 0.05 are stable
  C9  search_analysis_results on the IN-MEMORY frames (the only way it runs) -> search_results.tsv
      (CPB-06)

usage: python cellchat_corrected.py <run_root> <out_dir>
"""
import glob
import json
import os
import shutil
import sys
import time

import anndata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

RUN = sys.argv[1] if len(sys.argv) > 1 else '.'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'corrected'
os.makedirs(OUT, exist_ok=True)

DATA = os.path.join(RUN, 'data/cellphoneDB')
ZIP = os.path.join(DATA, 'v5.0.0/cellphonedb.zip')
META = os.path.join(DATA, 'metadata.tsv')
COUNTS = os.path.join(DATA, 'normalised_log_counts.h5ad')
MICRO = os.path.join(DATA, 'microenvironment.tsv')
TFS = os.path.join(DATA, 'active_TFs.tsv')
DEGS = os.path.join(DATA, 'DEGs_inv_trophoblast.tsv')

T0 = time.time()
log = []


def say(*a):
    line = '[%6.1fs] %s' % (time.time() - T0, ' '.join(str(x) for x in a))
    print(line, flush=True)
    log.append(line)


from cellphonedb.src.core.methods import cpdb_analysis_method            # noqa: E402
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method  # noqa: E402

adata = anndata.read_h5ad(COUNTS)
meta = pd.read_csv(META, sep='\t')
say('loaded', adata.shape, '| obs cols', list(adata.obs.columns))


# ------------------------------------------------------------------ C4
say('')
say('== C4: the barcode check of cell 9 ==')
tut = list(adata.obs.index).sort() == list(meta['barcode_sample']).sort()
say('tutorial expression  list(obs.index).sort() == list(barcode).sort()  ->  %r' % tut)
say('  (list.sort() returns None, so the cell compares None == None: it is True for ANY input)')
say('control: with the barcode column shuffled, the tutorial expression is still',
    list(adata.obs.index).sort() == list(meta['barcode_sample'].sample(frac=1, random_state=0)).sort())
ok = set(adata.obs_names) == set(meta['barcode_sample'])
say('real check  set(obs.index) == set(barcode_sample)  ->  %r' % ok)
say('  obs=%d  metadata=%d  shared=%d' % (len(adata.obs_names), len(meta),
                                          len(set(adata.obs_names) & set(meta['barcode_sample']))))

# ------------------------------------------------------------------ C2/C3 baseline rerun with stable names
say('')
say('== C2: re-run the two analyses with output_suffix so the artefacts have stable names ==')
for sub in ('simple', 'stat'):
    os.makedirs(os.path.join(OUT, sub), exist_ok=True)

if os.path.exists(os.path.join(OUT, 'table_simple_means_result.tsv')):
    say('  (reusing the persisted simple-analysis tables)')
    simple = {'means_result': pd.read_csv(os.path.join(OUT, 'table_simple_means_result.tsv'), sep='\t'),
              'deconvoluted': pd.read_csv(os.path.join(OUT, 'table_simple_deconvoluted.tsv'), sep='\t'),
              'deconvoluted_percents': pd.read_csv(os.path.join(OUT, 'table_simple_deconvoluted_percents.tsv'), sep='\t'),
              'interaction_scores': pd.read_csv(os.path.join(OUT, 'table_simple_interaction_scores.tsv'), sep='\t')}
else:
        simple = cpdb_analysis_method.call(
        cpdb_file_path=ZIP, meta_file_path=META, counts_file_path=COUNTS, counts_data='hgnc_symbol',
        microenvs_file_path=MICRO, score_interactions=True,
        output_path=os.path.join(OUT, 'simple') + '/',
        separator='|', threads=5, threshold=0.1, result_precision=3, debug=False,
        output_suffix='aswritten_t01',
    )
say('simple done (threshold 0.1):', list(simple.keys()))

# ------------------------------------------------------------------ C3 persist what cells 7-28 display
say('')
say('== C3: persist the display-only tables ==')
pairs_of_interest = ['PV MMP11|EVT_1', 'EVT_1|PV MMP11', 'PV MYH11|EVT_1', 'EVT_1|PV MYH11']
tables = {}
for name, df in simple.items():
    fn = os.path.join(OUT, 'table_simple_%s.tsv' % name)
    df.to_csv(fn, sep='\t', index=False)
    tables['simple_' + name] = list(df.shape)
    say('  %-22s %-12s -> %s' % (name, str(df.shape), os.path.basename(fn)))

if os.path.exists(os.path.join(OUT, 'table_stat_means.tsv')):
    say('  (reusing the persisted statistical tables)')
    stat = {k: pd.read_csv(os.path.join(OUT, 'table_stat_%s.tsv' % k), sep='\t')
            for k in ['deconvoluted', 'deconvoluted_percents', 'means', 'pvalues',
                      'significant_means', 'CellSign_active_interactions',
                      'CellSign_active_interactions_deconvoluted', 'interaction_scores']}
else:
    stat = cpdb_statistical_analysis_method.call(
    cpdb_file_path=ZIP, meta_file_path=META, counts_file_path=COUNTS, counts_data='hgnc_symbol',
    active_tfs_file_path=TFS, microenvs_file_path=MICRO, score_interactions=True,
    iterations=1000, threshold=0.1, threads=5, debug_seed=42, result_precision=3,
    pvalue=0.05, subsampling=False, subsampling_log=False, subsampling_num_pc=100,
    subsampling_num_cells=1000, separator='|', debug=False,
    output_path=os.path.join(OUT, 'stat') + '/', output_suffix='aswritten_t01_i1000',
    )
say('statistical done (0.1 / 1000 it):', list(stat.keys()))
for name, df in stat.items():
    fn = os.path.join(OUT, 'table_stat_%s.tsv' % name)
    df.to_csv(fn, sep='\t', index=False)
    tables['stat_' + name] = list(df.shape)
    say('  %-32s %-12s -> %s' % (name, str(df.shape), os.path.basename(fn)))


def sig_count(res):
    pcols = [c for c in res['significant_means'].columns if '|' in c]
    s = res['significant_means'][pcols].notna()
    return int(s.values.sum()), int(s.any(axis=1).sum()), len(pcols)


base_cells, base_ints, base_pairs = sig_count(stat)
say('baseline significant: %d (interaction,pair) cells | %d interactions | %d pairs'
    % (base_cells, base_ints, base_pairs))

# ------------------------------------------------------------------ C1 the typo
say('')
say('== C1: "PV MMPP11" -> "PV MMP11" and re-draw the PV-vs-trophoblast panel ==')
import ktplotspy as kpy  # noqa: E402

means, pvals = stat['means'], stat['pvalues']
CT2 = 'EVT_1|EVT_2|GC|iEVT|eEVT|VCT_CCC'
pairs_in_tables = set(p for p in stat['means'].columns if '|' in p)


def expected_pairs(ct1):
    """The ordered cell pairs ktplotspy will put on the x axis for this cell_type1 x cell_type2."""
    a, b = ct1.split('|'), CT2.split('|')
    want = ['%s|%s' % (x, y) for x in a for y in b] + ['%s|%s' % (y, x) for y in b for x in a]
    return [w for w in want if w in pairs_in_tables], [w for w in want if w not in pairs_in_tables]


for label, ct1 in [('aswritten_typo', 'PV MYH11|PV STEAP4|PV MMPP11'),
                   ('corrected', 'PV MYH11|PV STEAP4|PV MMP11')]:
    gg = kpy.plot_cpdb(
        adata=adata, cell_type1=ct1, cell_type2=CT2,
        means=means, pvals=pvals, celltype_key='cell_labels',
        genes=['TGFB2', 'CSF1R'], figsize=(10, 3),
        title='Interactions between\nPV and trophoblast', max_size=3, highlight_size=0.75,
        degs_analysis=False, standard_scale=True,
        interaction_scores=stat['interaction_scores'], scale_alpha_by_interaction_scores=True)
    fig = gg.draw()
    fig.savefig(os.path.join(OUT, '31_PV_trophoblast_%s.png' % label), dpi=200, bbox_inches='tight')
    fig.clf()
    ok_pairs, missing = expected_pairs(ct1)
    say('  %-14s %-34s pairs drawn: %d   silently dropped: %r'
        % (label, ct1, len(ok_pairs), missing))
    pd.DataFrame({'cell_pair': ok_pairs}).to_csv(
        os.path.join(OUT, '31_xaxis_pairs_%s.tsv' % label), sep='\t', index=False)

# is the misspelt label even present in the data?
labels = set(adata.obs['cell_labels'].astype(str))
say('  labels in matrix starting with PV:', sorted(l for l in labels if l.startswith('PV')))
say('  "PV MMPP11" in matrix: %r   "PV MMP11" in matrix: %r'
    % ('PV MMPP11' in labels, 'PV MMP11' in labels))

# ------------------------------------------------------------------ C5 threshold sensitivity
say('')
say('== C5: threshold sensitivity ==')
thr_rows = []
for thr in (0.05, 0.20):
    res = cpdb_statistical_analysis_method.call(
        cpdb_file_path=ZIP, meta_file_path=META, counts_file_path=COUNTS, counts_data='hgnc_symbol',
        active_tfs_file_path=TFS, microenvs_file_path=MICRO, score_interactions=True,
        iterations=1000, threshold=thr, threads=5, debug_seed=42, result_precision=3,
        pvalue=0.05, subsampling=False, subsampling_log=False, subsampling_num_pc=100,
        subsampling_num_cells=1000, separator='|', debug=False,
        output_path=os.path.join(OUT, 'stat') + '/',
        output_suffix='corrected_t%03d_i1000' % int(thr * 100))
    c, i_, p_ = sig_count(res)
    thr_rows.append({'threshold': thr, 'sig_cells': c, 'sig_interactions': i_, 'n_pairs': p_})
    say('  threshold %.2f -> %d (interaction,pair) cells, %d interactions' % (thr, c, i_))
    if thr == 0.05:
        res['significant_means'].to_csv(os.path.join(OUT, 'table_stat_significant_means_t005.tsv'),
                                        sep='\t', index=False)
thr_rows.insert(0, {'threshold': 0.10, 'sig_cells': base_cells,
                    'sig_interactions': base_ints, 'n_pairs': base_pairs})
pd.DataFrame(thr_rows).to_csv(os.path.join(OUT, 'threshold_sensitivity.tsv'), sep='\t', index=False)

# ------------------------------------------------------------------ C6 degs_analysis
say('')
say('== C6: degs_analysis=True with the shipped DEGs table (the dead input) ==')
degs = pd.read_csv(DEGS, sep='\t')
say('  degs file: %s columns=%s' % (str(degs.shape), list(degs.columns)))
degs_attempt = {}
try:
    gg = kpy.plot_cpdb(
        adata=adata, cell_type1='PV MYH11', cell_type2='EVT_1|EVT_2|GC|iEVT|eEVT|VCT_CCC',
        means=degs, pvals=None, celltype_key='cell_labels',
        genes=['TGFB2', 'CSF1R', 'COL1A1'], figsize=(12, 8),
        title='degs_analysis=True', max_size=6, highlight_size=0.75,
        degs_analysis=True, standard_scale=True)
    fig = gg.draw()
    fig.savefig(os.path.join(OUT, '32_degs_analysis.png'), dpi=200, bbox_inches='tight')
    fig.clf()
    degs_attempt = {'ok': True, 'x_pairs': int(gg.data['x'].nunique()),
                    'shapes': sorted(map(str, gg.data['shape'].unique()))}
    say('  -> worked; x pairs=%d shapes=%s' % (degs_attempt['x_pairs'], degs_attempt['shapes']))
except Exception as e:                                             # noqa: BLE001
    degs_attempt = {'ok': False, 'error': '%s: %s' % (type(e).__name__, e)}
    say('  -> FAILED: %s: %s' % (type(e).__name__, e))
    say('  (expected: plot_cpdb wants `meandf` in CellphoneDB DEG layout -- one row per'
        ' cell_type/gene with mean+pval columns; the shipped table is cluster/gene/logFC/P.Value)')

# ------------------------------------------------------------------ C7 no active TFs
say('')
say('== C7: statistical analysis without active_TFs.tsv ==')
stat_notf = cpdb_statistical_analysis_method.call(
    cpdb_file_path=ZIP, meta_file_path=META, counts_file_path=COUNTS, counts_data='hgnc_symbol',
    active_tfs_file_path=None, microenvs_file_path=MICRO, score_interactions=True,
    iterations=1000, threshold=0.1, threads=5, debug_seed=42, result_precision=3,
    pvalue=0.05, subsampling=False, subsampling_log=False, subsampling_num_pc=100,
    subsampling_num_cells=1000, separator='|', debug=False,
    output_path=os.path.join(OUT, 'stat') + '/', output_suffix='corrected_notf_t01_i1000')
c_notf = sig_count(stat_notf)
say('  keys: %s' % list(stat_notf.keys()))
say('  significant: %d cells / %d interactions  (with active TFs: %d / %d)'
    % (c_notf[0], c_notf[1], base_cells, base_ints))

# ------------------------------------------------------------------ C8 iterations
say('')
say('== C8: iterations 1000 -> 100 ==')
stat_100 = cpdb_statistical_analysis_method.call(
    cpdb_file_path=ZIP, meta_file_path=META, counts_file_path=COUNTS, counts_data='hgnc_symbol',
    active_tfs_file_path=TFS, microenvs_file_path=MICRO, score_interactions=True,
    iterations=100, threshold=0.1, threads=5, debug_seed=42, result_precision=3,
    pvalue=0.05, subsampling=False, subsampling_log=False, subsampling_num_pc=100,
    subsampling_num_cells=1000, separator='|', debug=False,
    output_path=os.path.join(OUT, 'stat') + '/', output_suffix='corrected_t01_i0100')
c100 = sig_count(stat_100)
pcols = [c for c in stat['significant_means'].columns if '|' in c]
s1000 = stat['significant_means'][pcols].notna()
s100 = stat_100['significant_means'][pcols].notna()
inter = int((s1000.values & s100.values).sum())
union = int((s1000.values | s100.values).sum())
say('  100 it: %d cells / %d interactions (1000 it: %d / %d)  Jaccard %.3f'
    % (c100[0], c100[1], base_cells, base_ints, inter / union if union else float('nan')))
p1000 = stat['pvalues'][pcols].values.ravel()
p100 = stat_100['pvalues'][pcols].values.ravel()
say('  p-value agreement: both <0.05 in %d cells; 1000<0.05 & 100>=0.05 in %d; vice versa in %d'
    % (int(((p1000 < .05) & (p100 < .05)).sum()),
       int(((p1000 < .05) & (p100 >= .05)).sum()),
       int(((p1000 >= .05) & (p100 < .05)).sum())))

# ------------------------------------------------------------------ C9 search from memory
say('')
say('== C9: search_analysis_results on the in-memory frames ==')
from cellphonedb.utils import search_utils  # noqa: E402
res = search_utils.search_analysis_results(
    query_cell_types_1=['EVT_1', 'EVT_2', 'GC', 'eEVT', 'iEVT'],
    query_cell_types_2=['PV MMP11', 'PV MYH11', 'PV STEAP4'],   # C1: correct spelling
    query_genes=['TGFBR1'], query_interactions=['CSF1_CSF1R'],
    significant_means=stat['significant_means'], deconvoluted=stat['deconvoluted'],
    interaction_scores=stat['interaction_scores'], query_minimum_score=50, separator='|',
    long_format=True, query_classifications=['Signaling by Transforming growth factor'])
say('  shape %s  columns %s' % (str(res.shape), list(res.columns)[:8]))
res.to_csv(os.path.join(OUT, 'search_results.tsv'), sep='\t', index=False)
if res.shape[0]:
    say(res.head(10).to_string())

# ------------------------------------------------------------------ summary
summary = {
    'baseline': {'sig_cells': base_cells, 'sig_interactions': base_ints, 'n_pairs': base_pairs},
    'no_active_tfs': {'sig_cells': c_notf[0], 'sig_interactions': c_notf[1]},
    'iterations_100': {'sig_cells': c100[0], 'sig_interactions': c100[1],
                       'jaccard_vs_1000': round(inter / union, 4) if union else None},
    'degs_analysis': degs_attempt,
    'barcode_check': {'tutorial_expr': tut, 'real_check': bool(ok)},
    'tables': tables,
}
json.dump(summary, open(os.path.join(OUT, 'corrected_summary.json'), 'w'), indent=2)
open(os.path.join(OUT, 'corrected_log.txt'), 'w').write('\n'.join(log) + '\n')
say('')
say('total %.1f s' % (time.time() - T0))
print('->', OUT)
