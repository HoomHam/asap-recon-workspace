"""Scratch: figure for 2steve/09 — sharpness bought vs noise paid (knobs), and CG iterations (same trade).
-> outputs/tune_gas/sharpness_vs_noise.png"""
import csv, os
import numpy as np
import matplotlib
import matplotlib.ticker
matplotlib.use('Agg')
import matplotlib.pyplot as plt

O = os.path.expanduser('~/Hooman/Work/Codes/2026_ASAP_Recon/workspace/outputs/tune_gas')
SUBJ = [('2024-11-13_025JC', '025JC (healthy)', '#2a78d6'), ('2024-10-02_011CN', '011CN (EBV)', '#eb6834')]
FAM = {'gplb': (['gplb150', 'gplb200', 'prod', 'gplb450', 'gplb700', 'gplb0'], 'o'),
       'kernel': (['k0.5', 'k1.0', 'k2.0'], 's'), 'MS': (['MS200', 'MS160'], '^')}
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'


def agg(rows, n, col):
    return np.mean([float(r[col]) for r in rows if r['setting'] == n and r['deapod'] == '0' and r.get(col, '') != ''])


fig, ax = plt.subplots(1, 2, figsize=(13, 5.2), facecolor='#fcfcfb')
fits = []
for key, lab, c in SUBJ:
    rows = list(csv.DictReader(open(os.path.join(O, key, 'metrics.csv'))))
    f0, n0 = agg(rows, 'prod', 'fwhm_bin0'), agg(rows, 'prod', 'noise')
    for fam, (names, mk) in FAM.items():
        x = np.array([agg(rows, n, 'fwhm_bin0') / f0 for n in names])
        y = np.array([agg(rows, n, 'noise') / n0 for n in names])
        ax[0].plot(x, y, marker=mk, ms=8, color=c, lw=1.5 if fam == 'gplb' else 0, mec='#fcfcfb', mew=1.5,
                   label=lab if fam == 'gplb' else None)
        if fam == 'gplb':
            fits.append(np.polyfit(np.log(x), np.log(y), 1)[0])
        if key == SUBJ[0][0]:
            for n, xx, yy in zip(names, x, y):
                t = {'gplb0': 'gplb off', 'prod': 'prod (gplb 300)'}.get(n, n)
                off = {'gplb450': (-52, -12), 'k0.5': (-30, 8), 'k1.0': (-10, 9), 'k2.0': (-34, -3), 'MS160': (8, 4),
                       'MS200': (8, -3)}.get(n, (7, 4))
                ax[0].annotate(t, (xx, yy), textcoords='offset points', xytext=off, fontsize=8, color=INK2)
xs = np.linspace(0.8, 1.4, 50)
ax[0].plot(xs, xs ** -1.5, ls='--', color=INK2, lw=1)
ax[0].annotate('uniform-density expectation\nnoise ∝ FWHM$^{-1.5}$', (0.86, 0.86 ** -1.5), textcoords='offset points',
               xytext=(-10, -34), fontsize=8, color=INK2)
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_xticks([0.8, 0.9, 1.0, 1.2, 1.4])
ax[0].set_xticklabels(['0.8', '0.9', '1.0', '1.2', '1.4'])
ax[0].set_yticks([0.5, 0.7, 1, 1.5, 2])
ax[0].set_yticklabels(['0.5', '0.7', '1', '1.5', '2'])
ax[0].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[0].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[0].set_xlabel('PSF FWHM relative to production (bin 0)', color=INK)
ax[0].set_ylabel('noise sd in lung relative to production', color=INK)
ax[0].set_title(f'Knobs: readout window (line), kernel (■), MS (▲)\nlog-log slope along gplb {np.mean(fits):.2f} (uniform density: −1.5)', fontsize=10, color=INK)
ax[0].legend(frameon=False, fontsize=8, loc='upper right')

for key, lab, c in SUBJ:
    rows = list(csv.DictReader(open(os.path.join(O, key, 'cg_control.csv'))))
    for case, ls in (('bin 0, noise-free', '--'), ('bin 0, matched noise', '-')):
        r = [x for x in rows if x['case'] == case and int(x['iters']) >= 5]
        ax[1].plot([int(x['iters']) for x in r], [100 * float(x['rms_err']) for x in r], ls=ls, marker='o', ms=5, color=c, lw=1.8,
                   label=f'{lab}, {"noise-free" if "free" in case else "with noise"}')
    m = list(csv.DictReader(open(os.path.join(O, key, 'metrics.csv'))))
    steve = 100 * float(next(x for x in m if x['setting'] == 'prod' and x['deapod'] == '0' and x['bin'] == '0')['total'])
    ax[1].axhline(steve, color=c, lw=1, ls=':')
    ax[1].annotate(f'Steve prod, bin 0: {steve:.1f} %', (5.2, steve), textcoords='offset points', xytext=(0, 3), fontsize=8, color=INK2)
ax[1].set_xscale('log')
ax[1].set_xticks([5, 10, 20, 40, 80])
ax[1].set_xticklabels(['5', '10', '20', '40', '80'])
ax[1].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[1].set_ylim(0, 45)
ax[1].set_xlabel('CG iterations (weighted least squares, no regulariser)', color=INK)
ax[1].set_ylabel('in-lung rms error vs object (% of lung mean)', color=INK)
ax[1].set_title('Same samples, exact operator: sharpening by iterating\nfinds the object without noise, amplifies noise with it', fontsize=10, color=INK)
ax[1].legend(frameon=False, fontsize=8, loc='upper right')
for a in ax:
    a.set_facecolor('#fcfcfb')
    a.grid(color=GRID, lw=0.8, which='both')
    for s in ('top', 'right'):
        a.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        a.spines[s].set_color(INK2)
    a.tick_params(colors=INK2, labelsize=8)
fig.tight_layout()
p = os.path.join(O, 'sharpness_vs_noise.png')
fig.savefig(p, dpi=130, facecolor='#fcfcfb')
print(p, 'gplb slopes', np.round(fits, 2))
