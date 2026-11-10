#!/usr/bin/env python3
"""cellchat_determinism.py -- is the 1000-permutation p-value reproducible?

Three 1000-iteration runs with the *same* debug_seed=42 and threads=5 already returned 4248, 4251 and
4248 significant (interaction, pair) cells, so something is not reproducible.  This script isolates the
variable: threads.  It runs the identical statistical analysis twice more, both with threads=1 (and
once with threads=5 for a same-configuration repeat), then hashes significant_means.

usage: python cellchat_determinism.py <run_root> <out_dir>
"""
import hashlib
import json
import os
import sys

import pandas as pd

RUN = sys.argv[1] if len(sys.argv) > 1 else '.'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'determinism'
os.makedirs(OUT, exist_ok=True)

DATA = os.path.join(RUN, 'data/cellphoneDB')
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method  # noqa: E402


def run(threads, tag):
    res = cpdb_statistical_analysis_method.call(
        cpdb_file_path=os.path.join(DATA, 'v5.0.0/cellphonedb.zip'),
        meta_file_path=os.path.join(DATA, 'metadata.tsv'),
        counts_file_path=os.path.join(DATA, 'normalised_log_counts.h5ad'),
        counts_data='hgnc_symbol',
        active_tfs_file_path=os.path.join(DATA, 'active_TFs.tsv'),
        microenvs_file_path=os.path.join(DATA, 'microenvironment.tsv'),
        score_interactions=True, iterations=1000, threshold=0.1, threads=threads,
        debug_seed=42, result_precision=3, pvalue=0.05, subsampling=False, subsampling_log=False,
        subsampling_num_pc=100, subsampling_num_cells=1000, separator='|', debug=False,
        output_path=os.path.join(OUT, 'runs') + '/',
        output_suffix='threads%d_%s' % (threads, tag))
    sm = res['significant_means'].sort_values('id_cp_interaction').reset_index(drop=True)
    pv = res['pvalues'].sort_values('id_cp_interaction').reset_index(drop=True)
    h = hashlib.md5(sm.to_csv(index=False).encode()).hexdigest()
    hp = hashlib.md5(pv.to_csv(index=False).encode()).hexdigest()
    pcols = [c for c in sm.columns if '|' in c]
    n = int(sm[pcols].notna().values.sum())
    print('%s threads=%d  sig_cells=%d  md5(sig)=%s  md5(pval)=%s' % (tag, threads, n, h[:12], hp[:12]),
          flush=True)
    return {'tag': tag, 'threads': threads, 'sig_cells': n, 'md5_sig': h, 'md5_pval': hp}


rows = [run(1, 'a'), run(1, 'b'), run(5, 'c')]
df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, 'determinism.tsv'), sep='\t', index=False)
same_1 = df.loc[df.threads == 1, 'md5_sig'].nunique() == 1
print('\nthreads=1 reproducible across 2 runs: %r' % same_1)
print('threads=5 differs from threads=1:     %r'
      % (df.loc[df.threads == 5, 'md5_sig'].iloc[0] != df.loc[df.threads == 1, 'md5_sig'].iloc[0]))
json.dump({'threads1_reproducible': bool(same_1), 'rows': rows},
          open(os.path.join(OUT, 'determinism.json'), 'w'), indent=2)
