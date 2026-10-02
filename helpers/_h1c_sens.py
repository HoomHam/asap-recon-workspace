"""Scratch: relative sensitivity of the two whitened virtual channels from LOW-RESOLUTION full-field images
(k-space Gaussian apodisation: high SNR, densely sampled, no fold-in), body mask and object extent.
Check: two-channel static least squares, channel-1 fit of the high-k data (the static part) in noise units."""
import os, sys, json, numpy as np
from scipy import ndimage as ndi
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H
NB = 200
ybar = np.load(os.path.join(H.EXT, 'ybar_train.npy')).astype(np.complex128)
D = H.load(); k = D['k']; op = H.Nufft(k, np.arange(H.NILV), n=NB)
dcf = np.load(os.path.join(H.EXT, 'dcf_n200.npy'))
kr = np.linalg.norm(k[:, H.KILL:], axis=-1)
bands = ((8, 48), (48, 198), (198, 510))
gf = lambda z, s: ndi.gaussian_filter(z.real, s) + 1j * ndi.gaussian_filter(z.imag, s)

def lowres(kc):
    ap = np.exp(-(kr / kc) ** 2)
    return np.stack([op.adj(ybar[c] * dcf * ap) for c in range(2)])

g = lowres(0.03)
a = np.sqrt((np.abs(g) ** 2).sum(0)); p = np.percentile(a, 99.5)
body = a > 0.10 * p
lab, n = ndi.label(body); sz = ndi.sum(body, lab, range(1, n + 1)); body = np.isin(lab, np.where(sz > 3000)[0] + 1)
body = ndi.binary_closing(body, iterations=2)
for ax in range(3):
    for i in range(NB):
        sl = [slice(None)] * 3; sl[ax] = i
        body[tuple(sl)] |= ndi.binary_fill_holes(body[tuple(sl)])
ext = {}
for ax, nm in enumerate('xyz'):
    idx = np.where(body.any(axis=tuple(b for b in range(3) if b != ax)))[0]
    ext[nm] = (float((idx.min() - NB / 2) * H.DX), float((idx.max() - NB / 2) * H.DX))
print('body voxels', int(body.sum()), 'extent mm local x/y/z', ext, flush=True)
np.save(os.path.join(H.EXT, 'body700.npy'), body)

def sens(kc, s):
    g = lowres(kc)
    num = g[1] * np.conj(g[0]) * body; den = np.abs(g[0]) ** 2 * body
    R = gf(num, s) / (ndi.gaussian_filter(den, s) + 1e-9 * den.max())      # normalised convolution: smooth fill outside the body
    nrm = np.sqrt(1 + np.abs(R) ** 2)
    return np.stack([1 / nrm, R / nrm]), R

def fit(S, it=30):
    Wt = dcf / dcf.sum()
    rhs = sum(np.conj(S[c]) * op.adj(Wt * ybar[c]) for c in range(2))
    sc = np.vdot(rhs, sum(np.conj(S[c]) * op.adj(Wt * op.fwd(S[c] * rhs)) for c in range(2))).real / np.vdot(rhs, rhs).real
    nrm = lambda x: sum(np.conj(S[c]) * op.adj(Wt * op.fwd(S[c] * x)) for c in range(2)) + 1e-3 * sc * x
    x = H.cg(nrm, rhs, it=it)
    e = [[(np.abs((op.fwd(S[c] * x) - ybar[c])[:, a_:b_]) ** 2).mean() * 13 for a_, b_ in bands] for c in range(2)]
    return x, e

best = None
for kc, s in ((0.03, 3.0), (0.03, 6.0), (0.015, 6.0), (0.05, 3.0), (0.05, 1.5)):
    S, R = sens(kc, s)
    x, e = fit(S)
    print(f'kc {kc} sigma {s}: |R| body median {np.median(np.abs(R)[body]):.3f} (5-95 pct {np.percentile(np.abs(R)[body], 5):.3f}-{np.percentile(np.abs(R)[body], 95):.3f}) '
          f'| resid ch0 {np.round(e[0], 2)} ch1 {np.round(e[1], 2)}', flush=True)
    score = e[1][1] + e[1][2] + e[0][1] + e[0][2]
    if best is None or score < best[0]: best = (score, kc, s, S, R, x)
_, kc, s, S, R, x = best
print('chosen kc', kc, 'sigma', s)
np.save(os.path.join(H.EXT, 'sens700.npy'), S.astype(np.complex64)); np.save(os.path.join(H.EXT, 'static700_ls.npy'), x.astype(np.complex64))
json.dump(dict(extent_mm=ext, kc=kc, sigma=s), open(os.path.join(H.OUT, 'sens700.json'), 'w'), indent=1)
fig, ax = plt.subplots(2, 4, figsize=(18, 9)); c = NB // 2; vm = 0.6 * np.percentile(np.abs(x), 99.8)
for j, (nm, im, kw) in enumerate((('low-res SOS + body mask', a, dict(cmap='gray', vmax=0.6 * p)), ('2-ch least squares |x|', np.abs(x), dict(cmap='gray', vmax=vm)),
                                 ('|R| = |S1/S0|', np.abs(R), dict(cmap='viridis', vmax=1.0)), ('arg R', np.angle(R), dict(cmap='twilight', vmin=-np.pi, vmax=np.pi)))):
    for r_, sl in enumerate((np.s_[c + 5, :, :], np.s_[:, c, :])):
        ax[r_, j].imshow(im[sl], **kw); ax[r_, j].contour(body[sl], [0.5], colors='y', linewidths=.4); ax[r_, j].set_title(nm, fontsize=9); ax[r_, j].axis('off')
fig.tight_layout(); fig.savefig(os.path.join(H.OUT, 'sens700.png'), dpi=60)
