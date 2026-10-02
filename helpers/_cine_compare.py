"""Scratch (2026-10-01): side-by-side cine of two binned recons (rows) x 10 coronal slices (columns).

recon.mat 'gas_phase' (bins, SI, AP, LR); |.| displayed; 10 AP slices spread over the lung span of the
first recon; each row gets its own window (99.5th pct over the shown slices). Loops the bins.

usage: _cine_compare.py <out.mp4> <label1> <recon1.mat> <label2> <recon2.mat> [--loops 3 --fps 4]
"""
import sys
import numpy as np
import scipy.io as sio
from scipy import ndimage as ndi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

out = sys.argv[1]
pairs = [(sys.argv[2], sys.argv[3]), (sys.argv[4], sys.argv[5])]
loops = int(sys.argv[sys.argv.index('--loops') + 1]) if '--loops' in sys.argv else 3
fps = int(sys.argv[sys.argv.index('--fps') + 1]) if '--fps' in sys.argv else 4
vols = [np.abs(sio.loadmat(f)['gas_phase']).astype(np.float32) for _, f in pairs]
nb, NSI, NAP, NLR = vols[0].shape

# lung span along AP from the first recon: dark voxels inside the per-coronal filled body hull
m = ndi.gaussian_filter(vols[0].mean(0), 1.5)
p99 = np.percentile(m, 99)
lungcount = []
for ap in range(NAP):
    hull = ndi.binary_fill_holes(m[:, ap, :] > 0.33 * p99)
    lungcount.append((hull & (m[:, ap, :] < 0.28 * p99)).sum())
lungcount = np.array(lungcount)
apc = np.where(lungcount > 0.25 * lungcount.max())[0]
aps = np.linspace(apc.min(), apc.max(), 10).round().astype(int)
# crop to the body bounding box (SI, LR) for bigger panels
body = (m > 0.2 * p99).any(axis=1)
si_idx, lr_idx = np.where(body)
si0, si1 = max(si_idx.min() - 2, 0), min(si_idx.max() + 3, NSI)
lr0, lr1 = max(lr_idx.min() - 2, 0), min(lr_idx.max() + 3, NLR)
vmax = [np.percentile(v[:, si0:si1, aps, lr0:lr1], 99.9) for v in vols]

fig, ax = plt.subplots(2, 10, figsize=(20, 5.6), gridspec_kw=dict(wspace=0.02, hspace=0.08))
ims = []
for r, v in enumerate(vols):
    row = []
    for c, ap in enumerate(aps):
        im = ax[r, c].imshow(v[0, si0:si1, ap, lr0:lr1], cmap='gray', vmin=0, vmax=vmax[r], origin='upper')
        ax[r, c].set_xticks([]); ax[r, c].set_yticks([])
        if r == 0:
            ax[r, c].set_title(f'AP {ap}', fontsize=8)
        row.append(im)
    ax[r, 0].set_ylabel(pairs[r][0], fontsize=10)
    ims.append(row)
ttl = fig.suptitle('', fontsize=11)
seq = list(range(nb)) * loops


def upd(i):
    b = seq[i]
    for r, v in enumerate(vols):
        for c, ap in enumerate(aps):
            ims[r][c].set_data(v[b, si0:si1, ap, lr0:lr1])
    ph = 'end-exp' if b in (0, 1, 2, 13, 14, 15) else ('end-insp' if b in (7, 8) else ('insp' if b < 7 else 'exp'))
    ttl.set_text(f'006KL 1H free-breathing cine, coronal slices — bin {b:2d}/15 ({ph})')
    return [im for row in ims for im in row] + [ttl]


FuncAnimation(fig, upd, frames=len(seq)).save(out, writer=FFMpegWriter(fps=fps, bitrate=6000))
print('wrote', out, 'AP slices', aps.tolist())
