#!/usr/bin/env python
"""Recompute the two evidence outputs that the positional-join bug (CPB-19) had poisoned.

Does not re-run CellphoneDB: reads the tables the as-written pass already wrote.
"""
import os
import pandas as pd

SIG = 'results/cellphoneDB/sig/statistical_analysis_'
OUT = 'evidence'
names = ['pvalues', 'means', 'significant_means', 'interaction_scores']
T = {n: pd.read_csv(SIG + n + '_11_09_2026_211716.txt', sep='\t') for n in names}
P = [c for c in T['pvalues'].columns if '|' in c]
sig = T['significant_means'][P + ['id_cp_interaction']].copy()
sig['n'] = sig[P].notna().sum(axis=1)
top = (sig[['id_cp_interaction', 'n']]
       .merge(T['pvalues'][['id_cp_interaction', 'interacting_pair', 'classification']],
              on='id_cp_interaction', how='left')
       .sort_values('n', ascending=False)[['interacting_pair', 'classification', 'n']]
       .rename(columns={'n': 'n_sig_pairs'}).head(200))
top.to_csv(os.path.join(OUT, 'sig_top_interactions.tsv'), sep='\t', index=False)

_ref = T['pvalues']['id_cp_interaction'].values
audit = [{'table': n, 'same_row_same_id': float((T[n]['id_cp_interaction'].values == _ref).mean())}
         for n in names]
pd.DataFrame(audit).to_csv(os.path.join(OUT, 'table_row_order.tsv'), sep='\t', index=False)
print(top.head(12).to_string(index=False))
print()
print(pd.DataFrame(audit).to_string(index=False))
