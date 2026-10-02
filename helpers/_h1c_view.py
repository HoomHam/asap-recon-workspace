"""Scratch: montage of 16-phase recons in the Tyger frame (SI, AP, LR). One row per recon:
coronal / sagittal / axial at phase 0 (end expiration) and phase 8 (end inspiration), and phase 8 - phase 0.
usage: _h1c_view.py <out.png> <label>=<file.npy|recon.mat> [...]   (--ap 55 --lr 32 --si 50 --win 0.6)"""
import sys, numpy as np, scipy.io as sio
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
args = sys.argv[1:]
opt = dict(ap=55, lr=32, si=50, win=0.6)
for k_ in list(opt):
    if f'--{k_}' in args:
        i = args.index(f'--{k_}'); opt[k_] = float(args[i + 1]); del args[i:i + 2]
out, items = args[0], [a.split('=', 1) for a in args[1:]]
ap, lr, si = int(opt['ap']), int(opt['lr']), int(opt['si'])
def load(f):
    if f.endswith('.mat'): return np.abs(sio.loadmat(f)['gas_phase']).astype(np.float32)
    if f.endswith('.npz'): return np.abs(np.load(f)['truth_ksphere']).astype(np.float32)
    return np.abs(np.load(f)).astype(np.float32)
fig, ax = plt.subplots(len(items), 7, figsize=(26, 3.9 * len(items)), squeeze=False)
for r, (lab, f) in enumerate(items):
    v = load(f); ref = np.percentile(v[0], 99.5); vm = opt['win'] * ref
    views = [(v[0][:, ap, :], 'coronal, phase 0'), (v[8][:, ap, :], 'coronal, phase 8'), (v[0][:, :, lr], 'sagittal, phase 0'), (v[8][:, :, lr], 'sagittal, phase 8'),
             (v[0][si, :, :], 'axial, phase 0'), (v[8][si, :, :], 'axial, phase 8')]
    for c, (im, t) in enumerate(views):
        ax[r, c].imshow(im, cmap='gray', vmin=0, vmax=vm); ax[r, c].set_title(f'{lab}: {t}', fontsize=8); ax[r, c].axis('off')
    ax[r, 6].imshow((v[8] - v[0])[:, ap, :], cmap='RdBu', vmin=-0.5 * vm, vmax=0.5 * vm); ax[r, 6].set_title('coronal, phase 8 - phase 0', fontsize=8); ax[r, 6].axis('off')
fig.tight_layout(); fig.savefig(out, dpi=55)
