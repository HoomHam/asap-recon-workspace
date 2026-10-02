"""Scratch: body support mask + sensitivity ratio on the working full-field grid (LR 176 x AP 120 x SI 160 at 3.5 mm,
local axis order x = LR, y = AP, z = SI). Low-resolution estimates."""
import os, sys, json, numpy as np
from scipy import ndimage as ndi
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H
G = (176, 120, 160)
ybar = np.load(os.path.join(H.EXT, 'ybar_train.npy')).astype(np.complex128)
D = H.load(); k = D['k']; op = H.Nufft(k, np.arange(H.NILV), n=G)
a1 = [np.cos(np.pi * (np.arange(n) - n / 2) / n) ** 2 for n in G]
h = a1[0][:, None, None] * a1[1][None, :, None] * a1[2][None, None, :]
w = np.ones(op.shape)
for i in range(25):
    w = w / np.abs(op.fwd(h * op.adj(w)))
dcf = w / w.sum(); np.save(os.path.join(H.EXT, 'dcf_g.npy'), dcf)
kr = np.linalg.norm(k[:, H.KILL:], axis=-1)
gf = lambda z, s: ndi.gaussian_filter(z.real, s) + 1j * ndi.gaussian_filter(z.imag, s)
g = np.stack([op.adj(ybar[c] * dcf * np.exp(-(kr / 0.04) ** 2)) for c in range(2)])
a = np.sqrt((np.abs(g) ** 2).sum(0)); p = np.percentile(a, 99)
print('low-res SOS: p99', p, 'background median / p99', np.median(a) / p)
body = a > 0.22 * p
lab, n = ndi.label(body); sz = ndi.sum(body, lab, range(1, n + 1))
body = np.isin(lab, np.where(sz > 1500)[0] + 1)
print('components kept', (sz > 1500).sum(), 'sizes', np.sort(sz)[::-1][:6])
# lungs / airways: closing with a large ball in a padded array, then slice-wise fill in the axial plane
pad = 12; bp = np.pad(body, pad)
bp = ndi.binary_closing(bp, structure=ndi.generate_binary_structure(3, 1), iterations=10)
body = bp[pad:-pad, pad:-pad, pad:-pad]
for i in range(G[2]): body[:, :, i] = ndi.binary_fill_holes(body[:, :, i])
body = ndi.binary_dilation(body, iterations=3)
ext = {}
for ax, nm in enumerate(('LR', 'AP', 'SI')):
    idx = np.where(body.any(axis=tuple(b for b in range(3) if b != ax)))[0]
    ext[nm] = (float((idx.min() - G[ax] / 2) * H.DX), float((idx.max() - G[ax] / 2) * H.DX))
print('support voxels', int(body.sum()), 'of', body.size, 'extent mm', ext)
num = g[1] * np.conj(g[0]) * body; den = np.abs(g[0]) ** 2 * body
R = gf(num, 3.0) / (ndi.gaussian_filter(den, 3.0) + 1e-9 * den.max())
print('|R| in support: median %.3f, 5-95 pct %.3f-%.3f' % tuple(np.percentile(np.abs(R)[body], [50, 5, 95])))
nrm = np.sqrt(1 + np.abs(R) ** 2)
S = np.stack([1 / nrm, R / nrm]).astype(np.complex64)
np.save(os.path.join(H.EXT, 'sens_g.npy'), S); np.save(os.path.join(H.EXT, 'support_g.npy'), body)
json.dump(dict(grid=G, extent_mm=ext, support_vox=int(body.sum())), open(os.path.join(H.OUT, 'support_g.json'), 'w'), indent=1)
fig, ax = plt.subplots(3, 3, figsize=(15, 13))
for j, (nm, im, kw) in enumerate((('low-res SOS + support', a, dict(cmap='gray', vmax=0.7 * p)), ('|R|', np.abs(R) * body, dict(cmap='viridis', vmax=1.0)), ('arg R', np.angle(R) * body, dict(cmap='twilight', vmin=-np.pi, vmax=np.pi)))):
    for r_, sl in enumerate((np.s_[G[0] // 2, :, :], np.s_[:, G[1] // 2, :], np.s_[:, :, G[2] // 2 - 15])):
        ax[r_, j].imshow(im[sl], **kw); ax[r_, j].contour(body[sl], [0.5], colors='y', linewidths=.5); ax[r_, j].set_title(nm, fontsize=9); ax[r_, j].axis('off')
fig.tight_layout(); fig.savefig(os.path.join(H.OUT, 'support_g.png'), dpi=55)
