"""Scratch: receiver-bandwidth hypothesis for the fixed misfit. The trajectory moves ~1.2/350 mm^-1 per 5 us sample, so
the readout band (+-0.5 cycles/sample) covers only +-150 mm along the instantaneous gradient; the body reaches +-300 mm.
1) speed along the readout; 2) noise and data spectra along the readout (receiver filter shape);
3) model prediction on a 4x finer time grid: how much of its energy lies outside the receiver band."""
import os, sys, json, numpy as np, finufft
from scipy.interpolate import CubicSpline
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY
G = DY.G; D = H.load(); k = D['k']; raw = D['raw']
sp = np.linalg.norm(np.diff(k, axis=1), axis=-1)                     # (832, 511) cycles/mm per sample
print('|dk/ds| x 350 mm (samples per 1/FOV350) at samples 2, 5, 10, 20, 50, 100, 200, 300, 400, 440, 480, 505:', np.round(sp.mean(0)[[2, 5, 10, 20, 50, 100, 200, 300, 400, 440, 480, 505]] * 350, 3))
print('-> half-band distance 0.5 / |dk/ds| (mm):', np.round(0.5 / sp.mean(0)[[10, 20, 50, 100, 200, 300, 400]], 0))
# 2) spectra along the readout, samples 60..440 (constant speed part), channel 0 raw
L = raw.shape[1]; npass = L // H.NILV
a = raw[0, :npass * H.NILV].reshape(npass, H.NILV, H.NPTS)[1:, :, 60:444]
nz = (a[1:] - a[:-1]).reshape(-1, a.shape[-1]); sg = a.reshape(-1, a.shape[-1])
win = np.hanning(a.shape[-1])[None, :]
Pn = (np.abs(np.fft.fftshift(np.fft.fft(nz[::7] * win, axis=1), axes=1)) ** 2).mean(0); Ps = (np.abs(np.fft.fftshift(np.fft.fft(sg[::7] * win, axis=1), axes=1)) ** 2).mean(0)
f = np.fft.fftshift(np.fft.fftfreq(a.shape[-1]))
edges = np.linspace(-0.5, 0.5, 11)
print('band (cycles/sample)      pass-difference PSD   data PSD   (each normalised to its mean)')
for lo, hi in zip(edges[:-1], edges[1:]):
    s_ = (f >= lo) & (f < hi); print(f'  {lo:+.2f} .. {hi:+.2f}:   {Pn[s_].mean()/Pn.mean():8.3f}          {Ps[s_].mean()/Ps.mean():8.3f}')
# 3) model on a 4x finer time grid
x = np.load(os.path.join(H.EXT, 'static_g_ls_ch0.npy')).astype(np.complex128)
OS = 4; s0 = np.arange(H.NPTS); sf = np.arange(0, H.NPTS - 1 + 1e-9, 1 / OS)
kf = CubicSpline(s0, k, axis=1)(sf)                                    # (832, nf, 3)
kap = 2 * np.pi * H.DX * kf.reshape(-1, 3)
yf = finufft.nufft3d2(kap[:, 0].copy(), kap[:, 1].copy(), kap[:, 2].copy(), x, isign=-1, eps=1e-6).reshape(H.NILV, -1)
seg = yf[:, 60 * OS:444 * OS]; ff = np.fft.fftshift(np.fft.fftfreq(seg.shape[1], 1 / OS))
Pm = (np.abs(np.fft.fftshift(np.fft.fft(seg * np.hanning(seg.shape[1])[None, :], axis=1), axes=1)) ** 2).mean(0)
out = Pm[np.abs(ff) > 0.5].sum() / Pm.sum()
print(f'model (support-limited image, continuous time): fraction of readout-spectrum energy OUTSIDE +-0.5 cycles/sample = {out:.4f}; in 0.4-0.5: {Pm[(np.abs(ff) > 0.4) & (np.abs(ff) <= 0.5)].sum()/Pm.sum():.4f}')
np.save(os.path.join(H.EXT, 'model_fine_ch0.npy'), yf.astype(np.complex64))
fig, ax = plt.subplots(1, 3, figsize=(18, 4.6))
ax[0].plot(sp.mean(0) * 350); ax[0].set_title('trajectory speed |dk/ds| (1/FOV350 per sample)'); ax[0].set_xlabel('sample')
ax[1].semilogy(f, Pn / Pn.mean(), label='pass difference (noise + motion)'); ax[1].semilogy(f, Ps / Ps.mean(), label='data'); ax[1].legend(); ax[1].set_title('spectrum along the readout (samples 60-444)'); ax[1].set_xlabel('cycles/sample')
ax[2].semilogy(ff, Pm / Pm.max()); ax[2].axvline(-0.5, c='r'); ax[2].axvline(0.5, c='r'); ax[2].set_title('model on 4x finer time grid (red = receiver band)'); ax[2].set_xlabel('cycles per 5 us sample')
fig.tight_layout(); fig.savefig(os.path.join(H.OUT, 'diag_bw.png'), dpi=60)
