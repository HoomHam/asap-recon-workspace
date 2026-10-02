"""Scratch (h1 challenge): full-field (700 mm, 200^3 at 3.5 mm) static image, two whitened channels.

All passes replay the same 832 interleaves, so the static problem lives on the unique set: data = per-interleave
mean over the training lines (steady state minus the shared held-out lines), weight = repeat count x DCF.
Steps: NUFFT Pipe DCF for the 700 mm field -> gridded image per channel -> smooth relative sensitivities ->
DCF-weighted least squares (CG) for one complex image. Prints the object's extent (to size the dynamic grid).
"""
import os, sys, json, time
import numpy as np
from scipy import ndimage as ndi
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H

NB = 200
KIT = np.load('/Volumes/HoomHamExt/Work/Codes/2026_XeCS_Recon/h1challenge/shared/judging_kit_v1.npz')
D = H.load(); k, vol = D['k'], D['vol']
v, info = H.prewhiten(D['raw'])
np.savez(os.path.join(H.EXT, 'whiten.npz'), noise_cov=info['noise_cov'], W=info['W'], U=info['U'])
L = v.shape[1]
train = np.setdiff1d(np.arange(H.NILV, L), KIT['heldout_lines'])
il = train % H.NILV
cnt = np.bincount(il, minlength=H.NILV).astype(float)
ybar = np.zeros((2, H.NILV, H.NPTS - H.KILL), np.complex128)
for c in range(2):
    np.add.at(ybar[c], il, v[c, train, H.KILL:])
ybar /= cnt[None, :, None]
print('train lines', train.size, 'repeats per interleave', cnt.min(), cnt.max())

op = H.Nufft(k, np.arange(H.NILV), n=NB)
# Pipe DCF by NUFFT: w <- w / |A h A^H w|, h = separable cos^2 image window (positive, compact kernel in k)
f = os.path.join(H.EXT, 'dcf_n200.npy')
if os.path.exists(f):
    dcf = np.load(f)
else:
    a = np.cos(np.pi * (np.arange(NB) - NB / 2) / NB) ** 2
    h = a[:, None, None] * a[None, :, None] * a[None, None, :]
    w = np.ones(op.shape)
    t0 = time.time()
    for i in range(25):
        den = np.abs(op.fwd(h * op.adj(w)))
        w = w / den
        if i % 5 == 4:
            print(f'  dcf it {i+1}: flatness sd/mean of A h A^H w = {np.std(den)/np.mean(den):.4f}  ({time.time()-t0:.0f} s)', flush=True)
    dcf = w / w.sum()
    np.save(f, dcf)
Wt = dcf * cnt[:, None]; Wt = Wt / Wt.sum()

t0 = time.time()
g = np.stack([op.adj(ybar[c] * Wt) for c in range(2)])                     # gridded, 700 mm
print('gridded 2 ch %.1f s' % (time.time() - t0))
# relative sensitivities: channel images referenced to the phase of channel 0, Gaussian low-pass, RSS-normalised
ref = np.exp(-1j * np.angle(g[0]))
lp = np.stack([ndi.gaussian_filter((g[c] * ref).real, 4.0) + 1j * ndi.gaussian_filter((g[c] * ref).imag, 4.0) for c in range(2)])
rss = np.sqrt((np.abs(lp) ** 2).sum(0))
S = lp / (rss + 1e-3 * rss.max())
body = ndi.gaussian_filter(np.sqrt((np.abs(g) ** 2).sum(0)), 2.0); body = body > 0.15 * np.percentile(body, 99.5)
lab, n = ndi.label(body); sz = ndi.sum(body, lab, range(1, n + 1)); body = np.isin(lab, np.where(sz > 2000)[0] + 1)
print('|S1|/|S0| inside body: median %.3f, 5-95 pct %.3f-%.3f' % tuple(np.percentile((np.abs(S[1]) / np.abs(S[0]))[body], [50, 5, 95])))
ext = {}
for ax, nm in enumerate('xyz'):
    idx = np.where(body.any(axis=tuple(a for a in range(3) if a != ax)))[0]
    ext[nm] = [float((idx.min() - NB / 2) * H.DX), float((idx.max() - NB / 2) * H.DX)]
print('body extent mm (local x, y, z):', ext, ' body voxels', int(body.sum()))

lam = 1e-3
def normal(x):
    out = lam * x
    for c in range(2):
        out = out + np.conj(S[c]) * op.adj(Wt * op.fwd(S[c] * x))
    return out
rhs = sum(np.conj(S[c]) * op.adj(Wt * ybar[c]) for c in range(2))
sc = np.vdot(rhs, normal(rhs)).real / np.vdot(rhs, rhs).real
lam = 1e-3 * sc                                                             # Tikhonov relative to the operator scale
res = []
x = H.cg(normal, rhs, it=15, cb=lambda i, x, r: (res.append(r), print(f'  cg {i+1}: rel resid {r:.4f}', flush=True)))
np.save(os.path.join(H.EXT, 'static700_ls.npy'), x.astype(np.complex64))
np.save(os.path.join(H.EXT, 'static700_grid.npy'), g.astype(np.complex64))
np.save(os.path.join(H.EXT, 'sens700.npy'), S.astype(np.complex64))
np.save(os.path.join(H.EXT, 'body700.npy'), body)
# data fit in noise units per readout band (training data, averaged repeats: noise var = 1/cnt per sample)
for c in range(2):
    r = op.fwd(S[c] * x) - ybar[c]
    for a, b in ((0, 48), (48, 198), (198, 510)):
        print(f'  ch{c} samples {a+2}-{b+2}: mean |resid|^2 x cnt = {(np.abs(r[:, a:b])**2 * cnt[:, None]).mean():.2f} (1 = noise)')
json.dump(dict(extent_mm=ext, cg_resid=res, train=int(train.size)), open(os.path.join(H.OUT, 'static700.json'), 'w'), indent=1)

# figure: centre slices, gridded-350 (channel SOS) vs gridded-700 vs LS-700, all cropped to the 350 box
o350 = H.Nufft(k, np.arange(H.NILV), n=100)
d100 = np.load(os.path.join(H.EXT, 'dcf_n100.npy')) * cnt[:, None]
g350 = np.sqrt(sum(np.abs(o350.adj(ybar[c] * d100)) ** 2 for c in range(2)))
c0, c1 = NB // 2 - 50, NB // 2 + 50
crop = lambda a: a[c0:c1, c0:c1, c0:c1]
ims = [('gridded, 350 mm grid (SOS)', g350), ('gridded, 700 mm grid, cropped (SOS)', crop(np.sqrt((np.abs(g) ** 2).sum(0)))), ('least squares 700 mm, 2-ch, cropped', crop(np.abs(x)))]
fig, ax = plt.subplots(3, 3, figsize=(13, 13))
for r, (nm, im) in enumerate(ims):
    t = H.to_tyger(im); vm = np.percentile(t, 99.5)
    for cidx, (sl, lab_) in enumerate(((np.s_[:, 55, :], 'coronal AP 55'), (np.s_[50, :, :], 'axial SI 50'), (np.s_[:, :, 35], 'sagittal LR 35'))):
        ax[r, cidx].imshow(t[sl], cmap='gray', vmin=0, vmax=0.6 * vm); ax[r, cidx].set_title(f'{nm}\n{lab_}', fontsize=9); ax[r, cidx].axis('off')
fig.tight_layout(); fig.savefig(os.path.join(H.OUT, 'static_350_vs_700.png'), dpi=80)
fig, ax = plt.subplots(1, 3, figsize=(15, 5.4))
for a_, axn in zip(ax, range(3)):
    a_.imshow(np.sqrt((np.abs(x) ** 2).mean(axis=axn)), cmap='gray'); a_.set_title(f'LS 700 mm, rms projection along local axis {axn}'); a_.axis('off')
fig.tight_layout(); fig.savefig(os.path.join(H.OUT, 'static700_proj.png'), dpi=70)
