#!/usr/bin/env python3
"""
Scratch (2026-09-13): evidence figure for 2steve note 04 (RBC/TP split).
  (a) logged R at the "solved" phase vs the spectral target RBCTPratio, all split sessions
  (b) 045VS EI bin: background (corner) noise of aRBC vs aTP -> correlation
  (c) 045VS EI bin: coronal slice of aRBC and aTP (TP negative inside lung)
Output: workspace/outputs/snr_2026-09-13/04_rbctp_split_evidence.png
"""
import csv
import importlib.util
import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy import ndimage

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / 'outputs' / 'snr_2026-09-13'
DYN = Path('/Volumes/HoomHamExt/AIkill_Dynamic')
spec = importlib.util.spec_from_file_location('snr', HERE / 'snr_calc.py')
S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)

rows = {r['key']: r for r in csv.DictReader(open(OUT / 'snr_table.csv'))}
pat = re.compile(r'RBC/TP phase solved at ph=([0-9.]+) rad \(R=([-0-9.]+) crossed target RBCTPratio=([-0-9.]+)\)')

# (a) logged R vs target
tgt, Rs, cls = [], [], []
for k, r in rows.items():
    if not r['dp_method'].startswith('split'):
        continue
    sol = [(float(R), float(t)) for _, R, t in pat.findall((DYN / k / 'd' / 'tyger.log').read_text(errors='replace'))]
    near = sum(abs(R / t - 1) < 0.5 for R, t in sol)
    c = 'valid' if near == len(sol) else ('invalid' if near == 0 else 'partial')
    for R, t in sol:
        tgt.append(t); Rs.append(R); cls.append(c)
tgt, Rs, cls = np.array(tgt), np.array(Rs), np.array(cls)

# (b,c) 045VS
key = '2024-09-10_045VS'
m = sio.loadmat(DYN / key / 'd' / 'recon.mat',
                variable_names=['gas_phase', 'dissolved_phase_real', 'dissolved_phase_imag'])
g = m['gas_phase'].astype(float); R = m['dissolved_phase_real'].astype(float); T = m['dissolved_phase_imag'].astype(float)
corners = S.corner_box(g.shape[1:]); gm = g.mean(0)
bg = corners & ~ndimage.binary_dilation(S.lung_mask(gm, gm[corners].std()), iterations=3)
b = int(rows[key]['insp_bin'])
mk = S.lung_mask(g[b], g[b][bg].std())
corr = np.corrcoef(R[b][bg], T[b][bg])[0, 1]

fig = plt.figure(figsize=(16, 4.8))
ax = fig.add_subplot(1, 4, 1)
colors = {'valid': 'tab:green', 'partial': 'tab:orange', 'invalid': 'tab:red'}
for c in ('invalid', 'partial', 'valid'):
    sel = cls == c
    ax.scatter(tgt[sel], np.clip(Rs[sel], 1e-3, None), s=6, alpha=0.5, color=colors[c],
               label=f'{c} ({len(set(zip(tgt[sel])))} target values)')
xx = np.array([0.05, 5]); ax.plot(xx, xx, 'k--', lw=1, label='R = target')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('spectral target RBCTPratio'); ax.set_ylabel('logged R at "solved" phase')
ax.set_title(f'(a) {len(Rs)} bins, 86 split sessions'); ax.legend(fontsize=7)

ax = fig.add_subplot(1, 4, 2)
ax.scatter(R[b][bg], T[b][bg], s=1, alpha=0.3)
ax.set_xlabel('aRBC (background voxels)'); ax.set_ylabel('aTP (background voxels)')
ax.set_title(f'(b) {key} bin {b}: noise corr {corr:+.3f}')
ax.set_aspect('equal', 'datalim')

y = int(np.argmax(mk.sum(axis=(0, 2))))
for i, (img, name) in enumerate(((R[b], 'aRBC (real)'), (T[b], 'aTP (imag)'))):
    ax = fig.add_subplot(1, 4, 3 + i)
    lim = np.percentile(np.abs(img[:, y, :]), 99.5)
    im = ax.imshow(img[:, y, :], cmap='RdBu_r', vmin=-lim, vmax=lim, origin='upper')
    ax.contour(mk[:, y, :], levels=[0.5], colors='k', linewidths=0.6)
    ax.set_title(f'({"c" if i == 0 else "d"}) {name}, lung mean {img[mk].mean():.0f}')
    ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(im, ax=ax, fraction=0.046)
fig.suptitle('Steve results.py RBC/TP split: phase stop criterion and conditioning (DIAPHRAGM recon outputs)')
fig.tight_layout()
OUT.mkdir(parents=True, exist_ok=True)
p = OUT / '04_rbctp_split_evidence.png'
fig.savefig(p, dpi=110)
print(f'bins plotted {len(Rs)}; {key} bin {b}: corr {corr:+.3f}, lung mean aRBC {R[b][mk].mean():.0f}, aTP {T[b][mk].mean():.0f}')
print(f'figure -> {p}')
