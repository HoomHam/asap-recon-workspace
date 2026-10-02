"""Scratch (h1 challenge): shared setup on the corrected model (thermal-noise whitening + receiver filter).
Pass-averaged training data, sensitivity ratio from low-resolution images, static two-channel check fit."""
import os, sys, json, numpy as np
from scipy import ndimage as ndi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H
G = (176, 120, 160)
KIT = np.load('/Volumes/HoomHamExt/Work/Codes/2026_XeCS_Recon/h1challenge/shared/judging_kit_v2.npz')
D = H.load(); k = D['k']; v, info = H.prewhiten(D['raw'])
np.savez(os.path.join(H.EXT, 'whiten.npz'), noise_cov=info['noise_cov'], W=info['W'], U=info['U'])
np.savez(os.path.join(H.EXT, 'shared', 'noise_whitening.npz'), noise_cov=info['noise_cov'], W=info['W'], U=info['U'])
train = KIT['scored_train_lines']; il = train % H.NILV
keep_il = np.setdiff1d(np.arange(H.NILV), KIT['heldout_interleaves'])
cnt = np.bincount(il, minlength=H.NILV).astype(float)
ybar = np.zeros((2, H.NILV, H.NPTS - H.KILL), np.complex128)
for c in range(2): np.add.at(ybar[c], il, v[c, train, H.KILL:])
ybar[:, keep_il] /= cnt[keep_il][None, :, None]
np.save(os.path.join(H.EXT, 'ybar_train.npy'), ybar.astype(np.complex64)); np.save(os.path.join(H.EXT, 'cnt_train.npy'), cnt)
op = H.RxNufft(k, n=G, nt=1, single=False, eps=1e-5)
dcf = np.load(os.path.join(H.EXT, 'dcf_g.npy')); M = np.load(os.path.join(H.EXT, 'support_g.npy'))
wsel = (cnt > 0)[:, None] * dcf; Wt = wsel / wsel.sum()
kr = np.linalg.norm(k[:, H.KILL:], axis=-1)
gf = lambda z, s: ndi.gaussian_filter(z.real, s) + 1j * ndi.gaussian_filter(z.imag, s)
g = np.stack([op.adj(ybar[c] * Wt * np.exp(-(kr / 0.04) ** 2)) for c in range(2)])
num = g[1] * np.conj(g[0]) * M; den = np.abs(g[0]) ** 2 * M
R = gf(num, 3.0) / (ndi.gaussian_filter(den, 3.0) + 1e-9 * den.max())
print('|R| in support: median %.3f, 5-95 pct %.3f-%.3f' % tuple(np.percentile(np.abs(R)[M], [50, 5, 95])))
nrm = np.sqrt(1 + np.abs(R) ** 2); S = np.stack([1 / nrm, R / nrm])
np.save(os.path.join(H.EXT, 'sens_g.npy'), S.astype(np.complex64))
bands = ((8, 48), (48, 198), (198, 510))
def fit(chs, tag, it=40, use_S=True):
    Sx = [S[c] if use_S else np.ones(G) for c in chs]
    rhs = M * sum(np.conj(s_) * op.adj(Wt * ybar[c]) for s_, c in zip(Sx, chs)); n0 = lambda z: M * sum(np.conj(s_) * op.adj(Wt * op.fwd(s_ * z)) for s_ in Sx)
    sc = np.vdot(rhs, n0(rhs)).real / np.vdot(rhs, rhs).real
    x = H.cg(lambda z: n0(z) + 1e-4 * sc * z, rhs, it=it)
    for s_, c in zip(Sx, chs):
        r = (op.fwd(s_ * x) - ybar[c])[keep_il]; y = ybar[c][keep_il]
        print(f'{tag:30s} ch{c}: misfit/signal power ' + ' '.join(f'{(np.abs(r[:, a:b])**2).sum()/(np.abs(y[:, a:b])**2).sum():.4f}' for a, b in bands) +
              '   misfit in thermal units x repeats ' + ' '.join(f'{(np.abs(r[:, a:b])**2 * cnt[keep_il][:, None]).mean():9.1f}' for a, b in bands), flush=True)
    return x
fit([0], 'ch0 alone', use_S=False); fit([1], 'ch1 alone', use_S=False)
x = fit([0, 1], 'two channels with R')
np.save(os.path.join(H.EXT, 'static_g_ls.npy'), x.astype(np.complex64))
