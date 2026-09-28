"""Scratch (s14): total error vs real-data SNR for the candidate recipes, both subjects, with the cohort's SNR
distribution underneath. -> outputs/tune_gas/snr_crossover.png"""
import csv, os
import numpy as np
import matplotlib
import matplotlib.ticker
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O = os.path.expanduser('~/Hooman/Work/Codes/2026_ASAP_Recon/workspace/outputs/tune_gas')
INK, INK2, GRID, BG = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
SER = [('prod', 1, 'MS 240 + deapod', '#2a78d6'), ('MS200', 1, 'MS 200 + deapod', '#eb6834'),
       ('MS160', 1, 'MS 160 + deapod', '#1baf7a'), ('best', 1, 'best LB + deapod (per SNR)', '#eda100')]
GL = ['gplb100', 'gplb150', 'gplb200', 'prod', 'gplb450', 'gplb700']
def load(d): return list(csv.DictReader(open(os.path.join(O, d, 'metrics.csv'))))
def agg(R, n, d, c):
    return np.mean([float(r[c]) for r in R if r['setting'] == n and r['deapod'] == str(d)])
coh = np.array([float(r['snr_med']) for r in csv.DictReader(open(os.path.join(O, 'cohort_gas_snr.csv')))])
fig = plt.figure(figsize=(12, 6), facecolor=BG)
gs = fig.add_gridspec(2, 2, height_ratios=[4, 1], hspace=0.08)
for j, (k, lab) in enumerate((('2024-11-13_025JC', '025JC (healthy)'), ('2024-10-02_011CN', '011CN (EBV)'))):
    ax = fig.add_subplot(gs[0, j]); axh = fig.add_subplot(gs[1, j], sharex=ax)
    xs, ys = [], {s[2]: [] for s in SER}
    for nm, suf in ((1, '_snr'), (1.5, '_snr_nm1.5'), (2, '_snr_nm2'), (3, '_snr_nm3'), (4, '_snr_nm4')):
        R = load(k + suf); p = agg(R, 'prod', 0, 'total')
        xs.append(agg(R, 'prod', 0, 'real_snr') / nm)
        for n, d, l, c in SER:
            v = min(agg(R, g, 1, 'total') for g in GL) if n == 'best' else agg(R, n, d, 'total')
            ys[l].append(v / p)
    for n, d, l, c in SER:
        ax.plot(xs, ys[l], marker='o', ms=6, lw=2, color=c, label=l, mec=BG, mew=1.2)
    ax.axhline(1, color=INK2, lw=1, ls='--')
    ax.annotate('production (MS 240, LB 300, no deapod)', (xs[0], 1), textcoords='offset points', xytext=(-210, 5), fontsize=8, color=INK2)
    ax.set_xscale('log'); ax.set_xlim(4.3, 26)
    ax.set_title(f'{lab}: noise scaled ×1 … ×4 on the real trajectory', fontsize=10, color=INK)
    ax.set_ylabel('total in-lung error / production' if j == 0 else '', color=INK)
    ax.set_ylim(0.75, 1.42)
    plt.setp(ax.get_xticklabels(), visible=False)
    if j == 0:
        ax.legend(frameon=False, fontsize=8, loc='upper right')
    axh.hist(coh, bins=np.logspace(np.log10(4.3), np.log10(26), 22), color='#9ec5f4', edgecolor=BG)
    axh.set_xlabel('real-data gas SNR (lung mean / far-corner σ, median of bins)', color=INK)
    axh.set_ylabel('sessions' if j == 0 else '', color=INK)
    axh.annotate(f'cohort n = {len(coh)}, median {np.median(coh):.1f}', (4.5, axh.get_ylim()[1] * 0.75), fontsize=8, color=INK2)
    for a in (ax, axh):
        a.set_facecolor(BG); a.grid(color=GRID, lw=0.8, which='both')
        for s in ('top', 'right'): a.spines[s].set_visible(False)
        for s in ('left', 'bottom'): a.spines[s].set_color(INK2)
        a.tick_params(colors=INK2, labelsize=8)
    axh.set_xticks([5, 7, 10, 13, 20]); axh.set_xticklabels(['5', '7', '10', '13', '20'])
    axh.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
fig.savefig(os.path.join(O, 'snr_crossover.png'), dpi=130, facecolor=BG, bbox_inches='tight')
print(os.path.join(O, 'snr_crossover.png'))
