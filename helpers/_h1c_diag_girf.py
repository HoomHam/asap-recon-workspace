"""Scratch: can a per-axis linear filter on the measured trajectory (delay, gain, eddy-current-like terms) or a
B0 eddy phase explain the fixed misfit? Global least squares on the residual of the support-constrained static fit
(channel 0), samples weighted by 1 / signal rms(sample) so that every part of the readout counts."""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY
G = DY.G
ybar = np.load(os.path.join(H.EXT, 'ybar_train.npy')).astype(np.complex128)[0]
k = H.load()['k']; op = H.Nufft(k, np.arange(H.NILV), n=G)
x = np.load(os.path.join(H.EXT, 'static_g_ls_ch0.npy')).astype(np.complex128)
yh = op.fwd(x); r = ybar - yh
ax_ = [np.arange(n) - n // 2 for n in G]
Dk = np.stack([op.fwd(-1j * ax_[a].reshape([-1 if i == a else 1 for i in range(3)]) * x) for a in range(3)])
kap = 2 * np.pi * H.DX * k                                      # radians/voxel, (832, 512, 3)
wt = 1 / np.sqrt((np.abs(ybar) ** 2).mean(0))                   # (510,)
sl = slice(8, 440); bands = ((8, 48), (48, 198), (198, 440))
lags = range(-6, 7)
def shifted(a, j):                                              # kappa_a(il, s - j) on the kept samples (sample index + KILL)
    idx = np.clip(np.arange(H.KILL, H.NPTS) - j, 0, H.NPTS - 1)
    return kap[:, idx, a]
def fit(cols, name):
    A_ = np.stack([np.concatenate([(c[:, sl] * wt[sl]).real.ravel(), (c[:, sl] * wt[sl]).imag.ravel()]) for c in cols], axis=1)
    b_ = np.concatenate([(r[:, sl] * wt[sl]).real.ravel(), (r[:, sl] * wt[sl]).imag.ravel()])
    p, *_ = np.linalg.lstsq(A_, b_, rcond=None)
    res = r - sum(pi * c for pi, c in zip(p, cols))
    fr = [1 - (np.abs(res[:, a:b]) ** 2).sum() / (np.abs(r[:, a:b]) ** 2).sum() for a, b in bands]
    print(f'{name:52s} explains {fr[0]*100:5.1f} / {fr[1]*100:5.1f} / {fr[2]*100:5.1f} % (samples 10-50 / 50-200 / 200-440)', flush=True)
    return p
print('misfit / signal power:', [round((np.abs(r[:, a:b]) ** 2).sum() / (np.abs(ybar[:, a:b]) ** 2).sum(), 4) for a, b in bands])
fit([Dk[a] * shifted(a, 0) for a in range(3)], 'gain per axis (3)')
fit([Dk[a] * (shifted(a, 1) - shifted(a, -1)) / 2 for a in range(3)], 'delay per axis (3)')
p = fit([Dk[a] * shifted(a, j) for a in range(3) for j in lags], 'FIR per axis, lags -6..6 (39)')
fit([Dk[a] * shifted(b, j) for a in range(3) for b in range(3) for j in (-1, 0, 1)], 'cross-axis 3x3, lags -1..1 (27)')
fit([1j * yh * shifted(a, j) for a in range(3) for j in lags], 'B0 eddy phase from each axis, lags -6..6 (39)')
fit([Dk[a] * shifted(a, j) for a in range(3) for j in lags] + [1j * yh * shifted(a, j) for a in range(3) for j in lags] + [Dk[a] for a in range(3)] + [yh, 1j * yh], 'all + constant k-offset + global gain (83)')
s = np.arange(510.0)
fit([1j * yh * s, 1j * yh * s ** 2, yh * s, yh * s ** 2], 'global phase / amplitude drift along the readout (4)')
fit([Dk[a] * s for a in range(3)] + [Dk[a] * s ** 2 for a in range(3)] + [Dk[a] for a in range(3)], 'common k-drift along the readout (9)')
