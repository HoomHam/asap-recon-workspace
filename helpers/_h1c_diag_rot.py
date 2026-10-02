"""Scratch: per-interleave rotation / scale / affine error of the trajectory as the fixed misfit? (errors that grow with |k|)"""
import os, sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY
G = DY.G
ybar = np.load(os.path.join(H.EXT, 'ybar_train.npy')).astype(np.complex128)[0]
k = H.load()['k']; op = H.Nufft(k, np.arange(H.NILV), n=G)
x = np.load(os.path.join(H.EXT, 'static_g_ls_ch0.npy')).astype(np.complex128)
yh = op.fwd(x); r = ybar - yh
ax_ = [np.arange(n) - n // 2 for n in G]
Dk = np.stack([op.fwd(-1j * ax_[a].reshape([-1 if i == a else 1 for i in range(3)]) * x) for a in range(3)])
kap = (2 * np.pi * H.DX * k[:, H.KILL:]).transpose(2, 0, 1)        # (3, 832, 510)
wt = 1 / np.sqrt((np.abs(ybar) ** 2).mean(0)); sl = slice(8, 440); bands = ((8, 48), (48, 198), (198, 440))
E = np.eye(3)
rot = [sum(Dk[a] * np.cross(E[j], kap, axisa=0, axisb=0, axisc=0)[a] for a in range(3)) for j in range(3)]
scl = [sum(Dk[a] * kap[a] for a in range(3))]
aff = [Dk[a] * kap[b] for a in range(3) for b in range(3)]
def per_il(cols, name):
    res = np.zeros_like(r); par = []
    for il in range(H.NILV):
        A_ = np.stack([np.concatenate([(c[il, sl] * wt[sl]).real, (c[il, sl] * wt[sl]).imag]) for c in cols], axis=1)
        b_ = np.concatenate([(r[il, sl] * wt[sl]).real, (r[il, sl] * wt[sl]).imag])
        p, *_ = np.linalg.lstsq(A_, b_, rcond=None); par.append(p)
        res[il] = r[il] - sum(pi * c[il] for pi, c in zip(p, cols))
    fr = [1 - (np.abs(res[:, a:b]) ** 2).sum() / (np.abs(r[:, a:b]) ** 2).sum() for a, b in bands]
    print(f'{name:44s} explains {fr[0]*100:5.1f} / {fr[1]*100:5.1f} / {fr[2]*100:5.1f} % (samples 10-50 / 50-200 / 200-440)', flush=True)
    return np.array(par)
# chance level: the same number of random real regressors with the same envelope
rng = np.random.default_rng(0)
rnd = [yh * np.exp(2j * np.pi * rng.random(yh.shape)) for _ in range(9)]
per_il(rnd[:3], 'chance level, 3 random regressors')
per_il(rnd, 'chance level, 9 random regressors')
pr = per_il(rot, 'rotation per interleave (3)')
ps = per_il(scl, 'scale per interleave (1)')
pa = per_il(aff, 'full 3x3 affine per interleave (9)')
print('fitted rotation (deg): rms per axis', np.round(np.degrees(np.sqrt((pr ** 2).mean(0))), 3), ' mean', np.round(np.degrees(pr.mean(0)), 3))
print('fitted scale: mean %.4f sd %.4f' % (ps.mean(), ps.std()))
np.savez(os.path.join(H.EXT, 'diag_rot_fits.npz'), pr=pr, ps=ps, pa=pa)
fig, ax = plt.subplots(1, 3, figsize=(18, 4.5))
for j in range(3): ax[0].plot(np.degrees(pr[:, j]), lw=.6, label=f'axis {j}')
ax[0].legend(); ax[0].set_title('fitted rotation per interleave (deg)'); ax[1].plot(ps, lw=.6); ax[1].set_title('fitted scale error per interleave')
ax[2].imshow(np.degrees(np.linalg.norm(pr, axis=1)).reshape(26, 32), aspect='auto'); ax[2].set_title('|rotation| (deg), 26 x 32')
fig.tight_layout(); fig.savefig(os.path.join(H.OUT, 'diag_rot.png'), dpi=60)
