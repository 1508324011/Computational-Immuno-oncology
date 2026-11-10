#!/usr/bin/env python
"""Figures 05-08 and 10 of the corrected edition: the quantities the tutorial only prints.

Reads the tables produced by the as-written pass (results/03_statistical, results/01_input),
the corrected pass (corrected/) and the DEG pass (degs/) -- it does not re-run CellphoneDB.

    python cellchat_plots.py <cellchat_run_20261109> <out_dir>
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

RUN = sys.argv[1] if len(sys.argv) > 1 else '.'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'plots'
RES = os.path.join(RUN, 'results', 'cellphoneDB')
SIG = os.path.join(RUN, 'results', 'cellphoneDB', 'sig')
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({'figure.dpi': 110, 'savefig.dpi': 200, 'font.size': 9,
                     'axes.titlesize': 10, 'axes.labelsize': 9.5,
                     'axes.spines.top': False, 'axes.spines.right': False})


def out(name):
    return os.path.join(OUT, name)


def say(msg):
    print('[plot] %s' % msg, flush=True)


# ---------------------------------------------------------------- fig 05: 9x9 matrix
mat = pd.read_csv(os.path.join(RUN, 'evidence', 'sig_matrix_9x9.tsv'), sep='\t', index_col=0)
order = list(mat.index)
fig, ax = plt.subplots(figsize=(6.6, 5.6))
im = ax.imshow(mat.values, cmap='magma', aspect='auto')
ax.set_xticks(range(len(order)))
ax.set_yticks(range(len(order)))
ax.set_xticklabels(order, rotation=45, ha='right')
ax.set_yticklabels(order)
for i in range(len(order)):
    for j in range(len(order)):
        v = mat.values[i, j]
        ax.text(j, i, str(v), ha='center', va='center', fontsize=7.5,
                color='white' if v < mat.values.max() * 0.62 else 'black')
ax.set_xlabel('receiver cell type (cell_type_b)')
ax.set_ylabel('sender cell type (cell_type_a)')
ax.set_title('Significant interactions, 9x9 ordered pairs\n'
             '{:,} of {:,} tested (interaction, pair) cells = {:.1f}%'
             .format(mat.values.sum(), 2440 * 81, 100.0 * mat.values.sum() / (2440 * 81)))
ax.set_xticks(np.arange(-.5, len(order), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(order), 1), minor=True)
ax.grid(which='minor', color='white', linewidth=0.6)
ax.tick_params(which='minor', length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
cb.set_label('significant interactions')
fig.text(0.02, 0.005,
         'row sums (sent): %s   |   column sums (received): %s'
         % (', '.join('%s %d' % (o, s) for o, s in zip(order, mat.values.sum(1))),
            ', '.join('%s %d' % (o, s) for o, s in zip(order, mat.values.sum(0)))),
         fontsize=5.6, color='0.35')
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(out('fig_05_sig_matrix_9x9.png'), bbox_inches='tight')
fig.clf()
say('fig_05_sig_matrix_9x9.png')

# ---------------------------------------------------------------- fig 06: p-value distribution
pv = pd.read_csv(os.path.join(SIG, 'statistical_analysis_pvalues_11_09_2026_211716.txt'), sep='\t')
P = [c for c in pv.columns if '|' in c]
v = pv[P].values.astype(float).ravel()
v = v[np.isfinite(v)]
counts = pd.Series({
    'p = 0': int((v == 0).sum()),
    '0<p\u22640.01': int(((v > 0) & (v <= 0.01)).sum()),
    '0.01<p\u22640.05': int(((v > 0.01) & (v <= 0.05)).sum()),
    '0.05<p<1': int(((v > 0.05) & (v < 1)).sum()),
    'p = 1': int((v == 1).sum()),
})
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6), gridspec_kw={'width_ratios': [1.35, 1]})
axes[0].hist(v, bins=np.linspace(0, 1, 101), color='#31688e')
axes[0].set_yscale('log')
axes[0].set_xlabel('permutation p-value (1000 label shuffles)')
axes[0].set_ylabel('(interaction, pair) cells (log)')
axes[0].set_title('{:,} tested cells; {:.2f}% sit at p = 1 because at least\n'
                  'one partner is below the 10% expression cut in that pair'
                  .format(len(v), 100.0 * (v == 1).mean()))
axes[0].axvline(0.05, color='#b5361c', lw=1.2, ls='--')
axes[0].text(0.055, axes[0].get_ylim()[1] * 0.25, 'p = 0.05', color='#b5361c', fontsize=8)
bars = axes[1].bar(range(len(counts)), counts.values,
                   color=['#440154', '#3b528b', '#21918c', '#5ec962', '#b8b8b8'])
axes[1].set_yscale('log')
axes[1].set_xticks(range(len(counts)))
axes[1].set_xticklabels(counts.index, fontsize=7, rotation=28, ha='right')
for b, c in zip(bars, counts.values):
    axes[1].text(b.get_x() + b.get_width() / 2, c * 1.15, '%d' % c, ha='center', fontsize=7)
axes[1].set_ylabel('cells (log)')
axes[1].set_ylim(1, counts.max() * 6)
axes[1].set_title('p-value buckets')
fig.tight_layout()
fig.savefig(out('fig_06_pvalue_distribution.png'), bbox_inches='tight')
fig.clf()
say('fig_06_pvalue_distribution.png')

# ---------------------------------------------------------------- fig 07: threshold sweep
th = pd.read_csv(os.path.join(RUN, 'corrected', 'threshold_sensitivity.tsv'), sep='\t')
th = th.sort_values('threshold')
fig, ax = plt.subplots(figsize=(5.6, 3.6))
x = np.arange(len(th))
b1 = ax.bar(x - 0.2, th.sig_cells.values, width=0.4, color='#31688e',
            label='significant (interaction, pair) cells')
b2 = ax.bar(x + 0.2, th.sig_interactions.values, width=0.4, color='#e8842c',
            label='interactions with >= 1 significant pair')
ax.set_xticks(x)
ax.set_xticklabels(['threshold %.2f' % t for t in th.threshold])
ax.set_ylabel('count')
ax.set_yscale('log')
ax.set_ylim(50, th.sig_cells.max() * 3)
for bars in (b1, b2):
    for b in bars:
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() * 1.08,
                '%d' % b.get_height(), ha='center', fontsize=7.5)
ax.legend(frameon=False, fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=1)
ax.set_title('Expression threshold drives almost every number that follows\n'
             '(everything else held fixed: 1000 permutations, p <= 0.05)')
fig.tight_layout()
fig.savefig(out('fig_07_threshold_sensitivity.png'), bbox_inches='tight')
fig.clf()
say('fig_07_threshold_sensitivity.png')

# ---------------------------------------------------------------- fig 08: DEG vs permutation
perm = pd.read_csv(os.path.join(SIG, 'statistical_analysis_significant_means_11_09_2026_211716.txt'),
                   sep='\t').set_index('id_cp_interaction')
degs = pd.read_csv(os.path.join(RUN, 'degs', 'degs_A_asis_significant_means.tsv'),
                   sep='\t').set_index('id_cp_interaction')
degs = degs.reindex(perm.index)
sp = perm[P].notna().values
sd = degs[P].notna().values
both_c, perm_c, degs_c = int((sp & sd).sum()), int((sp & ~sd).sum()), int((~sp & sd).sum())
both_i, perm_i, degs_i = int((sp.any(1) & sd.any(1)).sum()), int((sp.any(1) & ~sd.any(1)).sum()), \
    int((~sp.any(1) & sd.any(1)).sum())

fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.5))
sets = [(both_c, perm_c, degs_c, int(sp.sum()), int(sd.sum()), 'significant (interaction, pair) cells'),
        (both_i, perm_i, degs_i, int(sp.any(1).sum()), int(sd.any(1).sum()), 'significant interactions')]
for ax, (a, b, c, tot_a, tot_b, unit) in zip(axes, sets):
    vals = [a, b, c]
    bars = ax.bar([0, 1, 2], vals, color=['#440154', '#21918c', '#e8842c'], width=0.62)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, v * 1.02, '{:,}'.format(v),
                ha='center', va='bottom', fontsize=8.5)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['both', 'permutation\nonly', 'DEG list\nonly'], fontsize=8.5)
    ax.set_ylim(0, max(vals) * 1.22)
    ax.set_ylabel(unit.replace('signficant', 'significant'), fontsize=8)
    ax.set_title('%s\npermutation n = %s, DEG-based n = %s, Jaccard %.3f'
                 % (unit.split(' significant')[0].capitalize() if unit.startswith('significant (')
                    else 'Interactions',
                    '{:,}'.format(tot_a), '{:,}'.format(tot_b), a / float(a + b + c)), fontsize=9)
fig.suptitle('Two independent admissibility criteria agree only partially '
             '(both read the same 2,440 interactions x 81 pairs)', fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig(out('fig_08_degs_vs_permutation.png'), bbox_inches='tight')
fig.clf()
say('fig_08_degs_vs_permutation.png')

# ---------------------------------------------------------------- fig 10: top interactions
top = pd.read_csv(os.path.join(RUN, 'evidence', 'sig_top_interactions.tsv'), sep='\t').head(20)
top = top.iloc[::-1]
fam = top.classification.fillna('(unclassified)')
fams = list(dict.fromkeys(fam))
cmap = plt.get_cmap('tab10')
colmap = {f: cmap(i % 10) for i, f in enumerate(fams)}
fig, ax = plt.subplots(figsize=(7.2, 5.4))
ax.barh(np.arange(len(top)), top.n_sig_pairs.values, color=[colmap[f] for f in fam])
ax.set_yticks(np.arange(len(top)))
ax.set_yticklabels(top.interacting_pair.values, fontsize=7.5)
ax.set_xlabel('number of the 81 ordered cell-type pairs with p <= 0.05')
ax.set_xlim(0, 81)
ax.axvline(40.5, color='0.4', ls=':', lw=1)
ax.text(40.5, len(top) - 1.4, 'half of all pairs', color='0.4', fontsize=7.5, ha='center')
handles = [plt.Rectangle((0, 0), 1, 1, color=colmap[f]) for f in fams]
ax.legend(handles, fams, frameon=False, fontsize=6.6, loc='upper center',
          bbox_to_anchor=(0.5, -0.16), ncol=3, columnspacing=1.2, handlelength=1.1)
sm_all = pd.read_csv(os.path.join(SIG, 'statistical_analysis_significant_means_11_09_2026_211716.txt'),
                     sep='\t')
per_int_all = sm_all[P].notna().sum(1)
ax.set_title('Most widely significant interactions\n'
             '%d of %d interactions reach significance in at least one pair '
             '(median 0, mean %.1f)' % (int((per_int_all > 0).sum()), len(per_int_all),
                                        per_int_all.mean()), fontsize=9.5)
fig.tight_layout()
fig.savefig(out('fig_10_top_interactions.png'), bbox_inches='tight')
fig.clf()
say('fig_10_top_interactions.png')

# ---------------------------------------------------------------- fig 09: DEG coverage
sh = pd.read_csv(os.path.join(RUN, 'degs', 'A_asis', 'degs_used.tsv'), sep='\t')
cnt = sh.groupby('cluster').size()
types = list(mat.index)
fig, ax = plt.subplots(figsize=(6.0, 3.4))
vals = [int(cnt.get(t, 0)) for t in types]
cols = ['#31688e' if v else '#c9c9c9' for v in vals]
bars = ax.bar(range(len(types)), vals, color=cols)
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width() / 2, v + 60, '0' if v == 0 else '%d' % v,
            ha='center', fontsize=7.5)
ax.set_xticks(range(len(types)))
ax.set_xticklabels(types, rotation=45, ha='right')
ax.set_ylabel('genes in the shipped DEG list')
ax.set_ylim(0, max(vals) * 1.18)
ax.set_title('degs_analysis=True needs a DEG list per cell type.\n'
             'The shipped file covers %d of the 9 microenvironment types'
             % sum(1 for v in vals if v))
fig.tight_layout()
fig.savefig(out('fig_09_degs_coverage.png'), bbox_inches='tight')
fig.clf()
say('fig_09_degs_coverage.png')
say('done')
