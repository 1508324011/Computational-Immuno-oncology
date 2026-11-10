#!/usr/bin/env python3
"""cellchat_degs.py -- run the CellphoneDB DEG analysis the tutorial skips, and draw the panel that
`degs_analysis=True` is supposed to draw.

The tutorial ships DEGs_inv_trophoblast.tsv (cluster, gene, logFC, P.Value, adj.P.Val) and then
hard-codes `degs_analysis=False` in both plot_cpdb calls, so the file is never read.  It is however
exactly the input `cpdb_degs_analysis_method.call(degs_file_path=...)` wants: file_utils._load_degs
takes the FIRST TWO COLUMNS as (cluster, gene) and prints "DEGs expects 2 columns and got 5" for the
rest.  So the file is usable as shipped -- only the call is missing.

Two variants are run, because the shipped table is not a DEG *list*, it is a full limma result table
(13,809 cluster-gene rows, no filtering):
  A  as shipped        -- all 13,809 (cluster, gene) pairs, i.e. every gene is "differentially
                          expressed", which cannot be what the authors meant
  B  filtered          -- adj.P.Val < 0.05, which is what "DEGs" normally means
For each: run the DEG analysis, count significant interactions, draw the degs_analysis=True panel,
and compare the surviving gene sets with the limma ranking.

usage: python cellchat_degs.py <run_root> <out_dir>
"""
import os
import sys
import time

import anndata
import matplotlib
matplotlib.use('Agg')
import pandas as pd

RUN = sys.argv[1] if len(sys.argv) > 1 else '.'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'degs'
os.makedirs(OUT, exist_ok=True)

DATA = os.path.join(RUN, 'data/cellphoneDB')
T0 = time.time()


def say(*a):
    print('[%5.1fs] %s' % (time.time() - T0, ' '.join(str(x) for x in a)), flush=True)


from cellphonedb.src.core.methods import cpdb_degs_analysis_method  # noqa: E402

shipped = pd.read_csv(os.path.join(DATA, 'DEGs_inv_trophoblast.tsv'), sep='\t')
say('shipped DEGs table %s columns=%s' % (str(shipped.shape), list(shipped.columns)))
say('clusters in file: %d  (metadata labels: %d)'
    % (shipped['cluster'].nunique(), pd.read_csv(os.path.join(DATA, 'metadata.tsv'), sep='\t')['cell_type'].nunique()))
sig = shipped[shipped['adj.P.Val'] < 0.05]
say('rows with adj.P.Val < 0.05: %d (%.1f%%) across %d clusters'
    % (len(sig), 100.0 * len(sig) / len(shipped), sig['cluster'].nunique()))
say('per-cluster DEG counts (top 8): %s'
    % sig['cluster'].value_counts().head(8).to_dict())

# NOTE: the shipped table is already a filtered DEG *list* -- all 13,809 rows have
# adj.P.Val < 0.05, so "filter at 0.05" is a no-op and a second variant would be byte-identical.
variants = {'A_asis': shipped}
summary = {}
for tag, df in variants.items():
    ddir = os.path.join(OUT, tag)
    os.makedirs(ddir, exist_ok=True)
    fp = os.path.join(ddir, 'degs_used.tsv')
    df[['cluster', 'gene']].to_csv(fp, sep='\t', index=False)
    say('')
    say('== variant %s: %d (cluster, gene) rows -> %s ==' % (tag, len(df), fp))
    res = cpdb_degs_analysis_method.call(
        cpdb_file_path=os.path.join(DATA, 'v5.0.0/cellphonedb.zip'),
        meta_file_path=os.path.join(DATA, 'metadata.tsv'),
        counts_file_path=os.path.join(DATA, 'normalised_log_counts.h5ad'),
        counts_data='hgnc_symbol', degs_file_path=fp,
        microenvs_file_path=os.path.join(DATA, 'microenvironment.tsv'),
        active_tfs_file_path=os.path.join(DATA, 'active_TFs.tsv'), score_interactions=True,
        threshold=0.1, result_precision=3, threads=5, debug=False, separator='|',
        output_path=ddir + '/', output_suffix=tag)
    say('  keys: %s' % list(res.keys()))
    row = {}
    for k, v in res.items():
        if hasattr(v, 'shape'):
            v.to_csv(os.path.join(OUT, 'degs_%s_%s.tsv' % (tag, k)), sep='\t', index=False)
            row[k] = list(v.shape)
            say('    %-32s %s' % (k, str(v.shape)))
    if 'significant_means' in res:
        P = [c for c in res['significant_means'].columns if '|' in c]
        sm = res['significant_means'][P]
        row['sig_cells'] = int(sm.notna().values.sum())
        row['sig_interactions'] = int(sm.notna().any(axis=1).sum())
        row['n_pairs'] = len(P)
        say('  significant: %d (interaction,pair) cells / %d interactions / %d pairs'
            % (row['sig_cells'], row['sig_interactions'], row['n_pairs']))
        # the panel the tutorial never draws
        import ktplotspy as kpy  # noqa: E402
        # degs_analysis=True takes `pvals` = the relevant_interactions table (not None: ktplotspy
        # calls .copy() on it), and fills the shape gap from `means` itself.
        for ct1 in ('EVT_1', 'PV MYH11'):
            try:
                gg = kpy.plot_cpdb(
                    adata=anndata.read_h5ad(os.path.join(DATA, 'normalised_log_counts.h5ad')),
                    cell_type1=ct1, cell_type2='EVT_1|EVT_2|GC|iEVT|eEVT|VCT_CCC|PV MMP11|PV MYH11|PV STEAP4',
                    means=res['means'], pvals=res['relevant_interactions'],
                    celltype_key='cell_labels',
                    genes=['TGFB2', 'CSF1R', 'COL1A1'], figsize=(12, 6),
                    title='degs_analysis=True, %s (%s)' % (ct1, tag),
                    max_size=6, highlight_size=0.75, degs_analysis=True, standard_scale=True)
                fig = gg.draw()
                fn = os.path.join(OUT, '32_degs_analysis_%s_%s.png' % (tag, ct1.replace(' ', '_')))
                fig.savefig(fn, dpi=200, bbox_inches='tight')
                fig.clf()
                say('  wrote %s' % os.path.basename(fn))
            except Exception as e:  # noqa: BLE001
                say('  degs panel for cell_type1=%r FAILED: %s: %s' % (ct1, type(e).__name__, e))
    summary[tag] = row

# which genes survived into the DEG analysis, versus the limma ranking
if 'deconvoluted' in summary.get('A_asis', {}):
    d = pd.read_csv(os.path.join(OUT, 'degs_A_asis_deconvoluted.tsv'), sep='\t')
    say('')
    say('== DEG analysis vs the limma table ==')
    for ct in ('EVT_1', 'PV MYH11'):
        if ct not in d.columns:
            continue
        used = set(d[d[ct].notna()]['gene_name'])
        lim_sig = set(shipped[shipped['cluster'] == ct]['gene'])
        lim_top = set(shipped[shipped['cluster'] == ct].sort_values('adj.P.Val').head(200)['gene'])
        if not lim_sig:
            say('  %-10s no DEGs in the shipped file for this cluster' % ct)
            continue
        say('  %-10s cellphonedb kept %d genes; limma adj.P<0.05 %d genes; overlap %d (jaccard %.3f);'
            ' limma top-200 overlap %d'
            % (ct, len(used), len(lim_sig), len(used & lim_sig),
               len(used & lim_sig) / len(used | lim_sig) if used | lim_sig else float('nan'),
               len(used & lim_top)))

pd.DataFrame(summary).T.to_csv(os.path.join(OUT, 'degs_summary.tsv'), sep='\t')
say('')
say('total %.1f s' % (time.time() - T0))
