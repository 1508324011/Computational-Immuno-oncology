#!/usr/bin/env python3
"""cellchat_evidence.py -- turn the CellphoneDB run's output tables into the numbers this module's
docs quote.

The tutorial's cells 7-28 each end with a bare expression (`metadata.head(3)`,
`cpdb_results['means'].head(2)`, `search_results.head()`, ...), which a notebook displays through its
implicit-display hook.  A script drops them silently, so the as-written console contains none of
those tables (CPB-04).  This script reads the files the analysis actually wrote and records what
those cells would have shown, plus the aggregates the tutorial never computes at all.

usage: python cellchat_evidence.py <run_root> <out_dir>
"""
import glob
import json
import os
import sys

import pandas as pd

RUN = sys.argv[1] if len(sys.argv) > 1 else '.'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'results_tables'
os.makedirs(OUT, exist_ok=True)

DATA = os.path.join(RUN, 'data/cellphoneDB')
SCORE = sorted(glob.glob(os.path.join(RUN, 'results/cellphoneDB/score/*.txt')))
SIG = sorted(glob.glob(os.path.join(RUN, 'results/cellphoneDB/sig/*.txt')))


def load(paths, kind, exclude=None):
    """Pick the file of one kind; its name carries a timestamp (CPB-05), so glob + assert unique."""
    hits = [p for p in paths if kind in os.path.basename(p)]
    if exclude:
        hits = [p for p in hits if exclude not in os.path.basename(p)]
    if len(hits) != 1:
        raise SystemExit('expected exactly one %r file, got %r' % (kind, hits))
    return pd.read_csv(hits[0], sep='\t'), os.path.basename(hits[0])


log = []


def say(*a):
    line = ' '.join(str(x) for x in a)
    print(line)
    log.append(line)


# ---------------------------------------------------------------- input side
meta = pd.read_csv(os.path.join(DATA, 'metadata.tsv'), sep='\t')
import anndata  # noqa: E402
adata = anndata.read_h5ad(os.path.join(DATA, 'normalised_log_counts.h5ad'), backed='r')
obs = adata.obs['cell_labels'].astype(str)
micro = pd.read_csv(os.path.join(DATA, 'microenvironment.tsv'), sep='\t')
tfs = pd.read_csv(os.path.join(DATA, 'active_TFs.tsv'), sep='\t')
degs = pd.read_csv(os.path.join(DATA, 'DEGs_inv_trophoblast.tsv'), sep='\t')

say('== input files ==')
say('metadata.tsv            %s  columns=%s' % (meta.shape, list(meta.columns)))
say('microenvironment.tsv    %s  columns=%s' % (micro.shape, list(micro.columns)))
say('active_TFs.tsv          %s  columns=%s' % (tfs.shape, list(tfs.columns)))
say('DEGs_inv_trophoblast.tsv %s columns=%s' % (degs.shape, list(degs.columns)))

say('')
say('== who is in which microenvironment ==')
rows = []
for env, g in micro.groupby('microenvironment', sort=False):
    types = list(g['cell_type'])
    n = int(obs.isin(types).sum())
    rows.append({'microenvironment': env, 'n_types': len(types), 'cell_types': '|'.join(types),
                 'n_cells_from_metadata': n})
    say('%-18s %d types  %-55s %d cells' % (env, len(types), '|'.join(types), n))
micro_u = micro.drop_duplicates('cell_type')
pd.DataFrame(rows).to_csv(os.path.join(OUT, 'microenvironments.tsv'), sep='\t', index=False)

say('')
say('matrix obs cell_labels: %d distinct; microenvironment cell_types: %d distinct'
    % (obs.nunique(), micro_u.shape[0]))
say('cell types present in the matrix but NOT in any microenvironment: %d'
    % len(set(obs.unique()) - set(micro_u['cell_type'])))
say('  -> those cells are still in the count matrix; CellphoneDB scores their means but, with'
    ' microenvs given, never pairs them.')
say('cells inside the microenvironment types: %d of %d (%.1f%%)'
    % (int(obs.isin(set(micro_u['cell_type'])).sum()), len(obs),
       100.0 * obs.isin(set(micro_u['cell_type'])).sum() / len(obs)))
counts = obs.value_counts()
pd.DataFrame({'cell_type': counts.index, 'n_cells': counts.values}).to_csv(
    os.path.join(OUT, 'cell_type_counts.tsv'), sep='\t', index=False)

# ---------------------------------------------------------------- simple scoring
say('')
say('== simple (non-statistical) scoring ==')
t_means, f_means = load(SCORE, 'means_result')
t_deconv, f_deconv = load(SCORE, 'deconvoluted', exclude='percents')
t_pct, f_pct = load(SCORE, 'deconvoluted_percents')
t_scores, f_scores = load(SCORE, 'interaction_scores')
for name, df, fn in [('means_result', t_means, f_means), ('deconvoluted', t_deconv, f_deconv),
                     ('deconvoluted_percents', t_pct, f_pct), ('interaction_scores', t_scores, f_scores)]:
    say('%-24s %-8s from %s' % (name, str(df.shape), fn))

cell_cols = [c for c in t_means.columns if '|' in c]
say('cell-pair columns: %d  (per microenvironment: %d types -> %d ordered pairs)'
    % (len(cell_cols), micro_u.shape[0], micro_u.shape[0] ** 2))

# ---------------------------------------------------------------- statistical
say('')
say('== statistical analysis (1000 permutations, p<=0.05) ==')
t_p, f_p = load(SIG, 'pvalues')
t_m, f_m = load(SIG, 'statistical_analysis_means', exclude='result')
t_sm, f_sm = load(SIG, 'significant_means')
t_s, f_s = load(SIG, 'statistical_analysis_interaction_scores')
t_cs, f_cs = load(SIG, 'CellSign_active_interactions', exclude='deconvoluted')
say('pvalues            %s' % str(t_p.shape))
say('means              %s' % str(t_m.shape))
say('significant_means  %s' % str(t_sm.shape))
say('interaction_scores %s' % str(t_s.shape))
say('CellSign_active_interactions %s' % str(t_cs.shape))

pcols = [c for c in t_p.columns if '|' in c]
sig = t_sm[pcols].notna()
say('')
say('interactions scored          %d' % t_p.shape[0])
say('cell pairs                   %d' % len(pcols))
say('significant (interaction,pair) cells: %d of %d (%.1f%%)'
    % (int(sig.values.sum()), sig.size, 100.0 * sig.values.sum() / sig.size))
say('interactions with >=1 significant pair: %d / %d'
    % (int(sig.any(axis=1).sum()), sig.shape[0]))
say('cell pairs with >=1 significant interaction: %d / %d'
    % (int(sig.any(axis=0).sum()), sig.shape[1]))

per_pair = sig.sum(axis=0).sort_values(ascending=False)
per_int = sig.sum(axis=1).sort_values(ascending=False)
pd.DataFrame({'cell_pair': per_pair.index, 'n_significant': per_pair.values}).to_csv(
    os.path.join(OUT, 'sig_per_cell_pair.tsv'), sep='\t', index=False)
# NOTE (CPB-19): significant_means is NOT in the same row order as means/pvalues/
# interaction_scores -- only 0.0% of rows hold the same id_cp_interaction at the same
# position. Joining positionally (or by row label) silently attaches the counts of one
# interaction to the name of another. Always merge on id_cp_interaction.
top_int = pd.DataFrame({
    'id_cp_interaction': sig.index,
    'n_sig_pairs': per_int.values,
}).merge(t_p[['id_cp_interaction', 'interacting_pair', 'classification']],
        on='id_cp_interaction', how='left').head(200)
top_int = top_int[['interacting_pair', 'classification', 'n_sig_pairs']]
top_int.to_csv(os.path.join(OUT, 'sig_top_interactions.tsv'), sep='\t', index=False)

# row-order audit across the four statistical tables (evidence for CPB-19)
_ref = t_p['id_cp_interaction'].values
order_audit = []
for nm, tb in (('pvalues', t_p), ('means', t_m), ('significant_means', t_sm),
               ('interaction_scores', t_s)):
    order_audit.append({'table': nm,
                        'same_row_same_id': float((tb['id_cp_interaction'].values == _ref).mean())})
pd.DataFrame(order_audit).to_csv(os.path.join(OUT, 'table_row_order.tsv'), sep='\t', index=False)

# 9x9 matrix of significant-interaction counts, in microenv order
order = list(micro_u['cell_type'])
mat = pd.DataFrame(0, index=order, columns=order, dtype=int)
for col in pcols:
    a, b = col.split('|')
    if a in mat.index and b in mat.columns:
        mat.loc[a, b] = int(sig[col].sum())
mat.to_csv(os.path.join(OUT, 'sig_matrix_9x9.tsv'), sep='\t')

say('')
say('top 10 cell pairs by number of significant interactions:')
for k, v in per_pair.head(10).items():
    say('  %-30s %4d' % (k, v))
say('top 10 interactions by number of significant cell pairs:')
for i in per_int.index[:10]:
    say('  %-40s %-42s %3d' % (t_p.loc[i, 'interacting_pair'],
                               str(t_p.loc[i, 'classification'])[:40], per_int[i]))

# ---------------------------------------------------------------- CellSign (active TFs)
say('')
say('== CellSign_active_interactions (the active_TFs.tsv input) ==')
say(str(t_cs.head(10).to_string()))
t_cs.to_csv(os.path.join(OUT, 'cellsign_active_interactions.tsv'), sep='\t', index=False)

# ---------------------------------------------------------------- search_utils
say('')
say('== search_analysis_results (tutorial cell 28) ==')
say('Cannot be replayed from the saved files: search_utils.search_analysis_results() needs the')
say('LONG-format `deconvoluted` frame (it indexes deconvoluted["gene_name"]), while the file')
say('cellphonedb writes (statistical_analysis_deconvoluted_*.txt, %d columns) has no gene_name --'
    % t_deconv.shape[1])
say('it carries gene_a/gene_b instead.  In the notebook the call works because it is handed the')
say('IN-MEMORY frame; the tutorial never writes that frame to disk, so cell 28\'s table is')
say('unrecoverable from the artefacts (CPB-06).  Reproduced by re-running the analysis in memory:')
say('see cellchat_corrected.py step C9 -> search_results.tsv.')

# ---------------------------------------------------------------- DEGs (dead input)
say('')
say('== DEGs_inv_trophoblast.tsv (shipped, never used: degs_analysis=False in every plot) ==')
say('shape %s' % str(degs.shape))
say('first 3 columns: %s' % list(degs.columns[:3]))
say('deprecated pct.1/pct.2 style columns present: %s'
    % [c for c in degs.columns if c.startswith('pct')])

with open(os.path.join(OUT, 'evidence_log.txt'), 'w') as fh:
    fh.write('\n'.join(log) + '\n')
json.dump({'n_cell_pairs': len(pcols), 'n_interactions': int(t_p.shape[0]),
           'n_sig_cells': int(sig.values.sum()), 'n_sig_interactions': int(sig.any(axis=1).sum()),
           'n_types_microenv': int(micro_u.shape[0]), 'n_types_metadata': int(obs.nunique()),
           'n_cells': int(len(obs)),
           'n_cells_in_microenv': int(obs.isin(set(micro_u['cell_type'])).sum())},
          open(os.path.join(OUT, 'evidence_summary.json'), 'w'), indent=2)
print('\n-> %s' % OUT)
