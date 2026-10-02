"""Scratch (h1 challenge): side-by-side cine of 16-phase recons in the Tyger frame (SI, AP, LR).
One row per recon; columns: coronal at three AP positions, sagittal through each lung, axial. Each row has its own
window (0 .. win x 99.5th percentile). usage: _h1c_cine.py <out.mp4> <label>=<file> [...] [--win 0.5] [--fps 5] [--loops 3]"""
import sys, numpy as np, scipy.io as sio
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
args = sys.argv[1:]
def opt(name, default):
    if name in args:
        i = args.index(name); v = float(args[i + 1]); del args[i:i + 2]; return v
    return default
win, fps, loops = opt('--win', 0.5), int(opt('--fps', 5)), int(opt('--loops', 3))
out, items = args[0], [a.split('=', 1) for a in args[1:]]
def load(f):
    if f.endswith('.mat'): return np.abs(sio.loadmat(f)['gas_phase']).astype(np.float32)
    if f.endswith('.npz'): return np.abs(np.load(f)['truth_ksphere']).astype(np.float32)
    return np.abs(np.load(f)).astype(np.float32)
vols = [load(f) for _, f in items]
views = [(lambda v, b: v[b][18:96, 45, 4:96], 'coronal AP 45'), (lambda v, b: v[b][18:96, 53, 4:96], 'coronal AP 53'), (lambda v, b: v[b][18:96, 61, 4:96], 'coronal AP 61'),
         (lambda v, b: v[b][18:96, 22:92, 31], 'sagittal LR 31'), (lambda v, b: v[b][18:96, 22:92, 68], 'sagittal LR 68'), (lambda v, b: v[b][48, 18:92, 4:96], 'axial SI 48')]
fig, ax = plt.subplots(len(vols), len(views), figsize=(3.3 * len(views), 3.0 * len(vols)), squeeze=False, gridspec_kw=dict(wspace=0.02, hspace=0.10))
ims = []
for r, v in enumerate(vols):
    vm = win * np.percentile(v, 99.5); row = []
    for c, (f, t) in enumerate(views):
        row.append(ax[r, c].imshow(f(v, 0), cmap='gray', vmin=0, vmax=vm, interpolation='nearest')); ax[r, c].set_xticks([]); ax[r, c].set_yticks([])
        if r == 0: ax[r, c].set_title(t, fontsize=9)
    ax[r, 0].set_ylabel(items[r][0], fontsize=9); ims.append(row)
ttl = fig.suptitle('', fontsize=10)
seq = list(range(16)) * loops
def upd(i):
    b = seq[i]
    for r, v in enumerate(vols):
        for c, (f, _) in enumerate(views): ims[r][c].set_data(f(v, b))
    ttl.set_text(f'phase {b} of 16 (0, 15 = end expiration; 7, 8 = end inspiration)')
    return [m for row in ims for m in row] + [ttl]
FuncAnimation(fig, upd, frames=len(seq), blit=False).save(out, writer=FFMpegWriter(fps=fps, bitrate=6000), dpi=80)
print('wrote', out)
