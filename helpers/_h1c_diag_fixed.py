"""Scratch: structure of the fixed-pattern misfit (channel 0, support-constrained static model, pass-averaged data).
Per-interleave fits of candidate calibration errors: complex gain, k-offset (3), readout time shift, phase ramp in time."""
import os, sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY
G = DY.G
ybar = np.load(os.path.join(H.EXT, 'ybar_train.npy')).astype(np.complex128)[0]
k = H.load()['k']; op = H.Nufft(k, np.arange(H.NILV), n=G)
x = np.load(os.path.join(H.EXT, 'static_g_ls_ch0.npy')).astype(np.complex128)
yh = op.fwd(x); r = ybar - yh
kk = k[:, H.KILL:]; kr = np.linalg.norm(kk, axis=-1)
print('|k| (cycles/mm) at samples 10/50/100/200/300/400/500: mean', np.round(kr[:, [8, 48, 98, 198, 298, 398, 498]].mean(0), 4), ' sd over interleaves', np.round(kr[:, [8, 48, 98, 198, 298, 398, 498]].std(0), 5))
print('k at sample 2 (first kept): mean vector', np.round(kk[:, 0].mean(0), 5), 'sd', np.round(kk[:, 0].std(0), 5), ' 1/FOV350 =', round(1 / 350, 5))
ax_ = [np.arange(n) - n // 2 for n in G]
Dk = []
for a in range(3):
    ra = ax_[a].reshape([-1 if i == a else 1 for i in range(3)])
    Dk.append(op.fwd(-1j * ra * x))                           # d yhat / d kappa_a (radians per voxel)
Dk = np.stack(Dk)
Dt = np.gradient(yh, axis=1)                                   # d yhat / d sample
s = np.arange(510, dtype=float)
bands = ((8, 48), (48, 198), (198, 510))
def explain(basis_fn, name, real=True):
    res = np.zeros_like(r); par = []
    for il in range(H.NILV):
        B = basis_fn(il)                                       # (nparam, 510) complex
        sl = slice(8, 510)
        if real:
            A_ = np.concatenate([B[:, sl].real, B[:, sl].imag], axis=1).T; b_ = np.concatenate([r[il, sl].real, r[il, sl].imag])
            p, *_ = np.linalg.lstsq(A_, b_, rcond=None)
        else:
            p, *_ = np.linalg.lstsq(B[:, sl].T, r[il, sl], rcond=None)
        res[il] = r[il] - p @ B; par.append(p)
    par = np.array(par)
    fr = [1 - (np.abs(res[:, a:b]) ** 2).sum() / (np.abs(r[:, a:b]) ** 2).sum() for a, b in bands]
    print(f'{name:38s} explains {fr[0]*100:5.1f} % / {fr[1]*100:5.1f} % / {fr[2]*100:5.1f} % of the misfit power (samples 10-50 / 50-200 / 200-512)')
    return par, res
print('misfit / signal power per band:', [round((np.abs(r[:, a:b]) ** 2).sum() / (np.abs(ybar[:, a:b]) ** 2).sum(), 4) for a, b in bands])
pg, _ = explain(lambda il: yh[il][None], 'complex gain per interleave', real=False)
pk, rk = explain(lambda il: Dk[:, il], 'k-offset per interleave (3)')
pt, _ = explain(lambda il: Dt[il][None], 'time shift per interleave')
pp, _ = explain(lambda il: np.stack([1j * yh[il], 1j * s * yh[il]]), 'phase offset + ramp per interleave')
pa, _ = explain(lambda il: np.concatenate([Dk[:, il], Dt[il][None], np.stack([1j * yh[il], 1j * s * yh[il], yh[il]])]), 'all of the above (7 real)')
dkap = 2 * np.pi * H.DX                                        # radians/voxel per (cycle/mm)
print('fitted k-offset (units of 1/FOV350): rms per axis', np.round(np.sqrt((pk ** 2).mean(0)) / dkap * 350, 3), ' mean', np.round(pk.mean(0) / dkap * 350, 3))
print('fitted time shift (samples): mean %.3f sd %.3f' % (pt.mean(), pt.std()))
print('fitted |gain-1|: |g| rms %.4f ; phase ramp over readout (rad): mean %.3f sd %.3f' % (np.sqrt((np.abs(pg) ** 2).mean()), (pp[:, 1] * 510).mean(), (pp[:, 1] * 510).std()))
np.savez(os.path.join(H.EXT, 'diag_fixed_fits.npz'), pk=pk, pt=pt, pg=pg, pp=pp, pa=pa)
e = (np.abs(r[:, 48:198]) ** 2).sum(1) / (np.abs(ybar[:, 48:198]) ** 2).sum(1)
fig, ax = plt.subplots(2, 3, figsize=(17, 8))
ax[0, 0].imshow(e.reshape(26, 32), aspect='auto'); ax[0, 0].set_title('relative misfit (samples 50-200) per interleave, 26 x 32')
ax[0, 1].imshow(e.reshape(32, 26), aspect='auto'); ax[0, 1].set_title('same, 32 x 26')
ax[0, 2].plot(e); ax[0, 2].set_title('relative misfit against interleave index')
for a in range(3): ax[1, 0].plot(pk[:, a] / dkap * 350, lw=.6, label='xyz'[a])
ax[1, 0].legend(); ax[1, 0].set_title('fitted k-offset per interleave (1/FOV350)')
ax[1, 1].plot(pt, lw=.6); ax[1, 1].set_title('fitted time shift per interleave (samples)')
ax[1, 2].semilogy(np.sqrt((np.abs(r) ** 2).mean(0)), label='misfit rms'); ax[1, 2].semilogy(np.sqrt((np.abs(ybar) ** 2).mean(0)), label='signal rms'); ax[1, 2].legend(); ax[1, 2].set_title('against readout sample')
fig.tight_layout(); fig.savefig(os.path.join(H.OUT, 'diag_fixed.png'), dpi=60)
