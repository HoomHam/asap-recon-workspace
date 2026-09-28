"""Scratch: side-by-side summary of tune_gas.py runs (one row per setting, columns per subject).
usage: python _tune_gas_combine.py KEY1 KEY2 ...  -> outputs/tune_gas/combined_{table.md, total.png}"""
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = os.path.expanduser('~/Hooman/Work/Codes/2026_ASAP_Recon/workspace/outputs/tune_gas')
keys = sys.argv[1:]
data = {}
for k in keys:
    with open(os.path.join(OUT, k, 'metrics.csv')) as f:
        data[k] = list(csv.DictReader(f))
names = list(dict.fromkeys(r['setting'] for r in data[keys[0]]))


def agg(k, n, dea, col):
    v = [float(r[col]) for r in data[k] if r['setting'] == n and r['deapod'] == str(dea) and r.get(col, '') != '']
    return float(np.mean(v)) if v else np.nan


cols = [('fwhm_bin0', 'FWHM', '{:.2f}'), ('noise', 'noise', '{:.3f}'), ('alias', 'alias', '{:.3f}'),
        ('blur', 'blur', '{:.3f}'), ('total', 'total', '{:.3f}'), ('holes', 'holes', '{:.2f}')]
with open(os.path.join(OUT, 'combined_table.md'), 'w') as f:
    f.write('| setting | deapod | ' + ' | '.join(f'{lab} {k[-5:]}' for c, lab, _ in cols for k in keys) + ' |\n')
    f.write('|---|---|' + '---|' * (len(cols) * len(keys)) + '\n')
    for n in names:
        for dea in (0, 1):
            f.write(f'| {n} | {dea} | ' + ' | '.join(fmt.format(agg(k, n, dea, c)) for c, lab, fmt in cols for k in keys) + ' |\n')

fig, ax = plt.subplots(1, 2, figsize=(14, 5))
x = np.arange(len(names))
w = 0.8 / (2 * len(keys))
series = ['#2563EB', '#EA580C', '#16A34A', '#6B7280']
for i, k in enumerate(keys):
    for dea in (0, 1):
        j = 2 * i + dea
        v = [agg(k, n, dea, 'total') for n in names]
        ax[0].bar(x - 0.4 + (j + 0.5) * w, v, w, color=series[i], alpha=1.0 if dea == 0 else 0.45,
                  label=f'{k[-5:]}' + (' deapod' if dea else ''))
    p = agg(k, 'prod', 0, 'total')
    ax[0].axhline(p, color=series[i], lw=0.8, ls='--')
ax[0].set_xticks(x)
ax[0].set_xticklabels(names, rotation=45, ha='right', fontsize=8)
ax[0].set_ylabel('total in-lung rms error / lung mean (mean of bins)')
ax[0].set_title('T-image criterion (dashed = production)')
ax[0].legend(frameon=False, fontsize=8)
for i, k in enumerate(keys):
    xs = [agg(k, n, 0, 'fwhm_bin0') for n in names]
    ys = [agg(k, n, 0, 'noise') for n in names]
    ax[1].plot(xs, ys, 'o', color=series[i], label=k[-5:])
    for n, xx, yy in zip(names, xs, ys):
        if i == 0:
            ax[1].annotate(n, (xx, yy), textcoords='offset points', xytext=(4, 3), fontsize=7)
ax[1].set_xlabel('PSF FWHM bin 0 (voxels of 3.5 mm)')
ax[1].set_ylabel('noise sd / lung mean')
ax[1].set_title('resolution vs noise')
ax[1].legend(frameon=False, fontsize=8)
for a in ax:
    a.grid(alpha=0.25)
    for s in ('top', 'right'):
        a.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'combined_total.png'), dpi=120)
print('->', OUT)
