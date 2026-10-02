"""Scratch: estimate the receiver filter from the data (FIR on the 4x fine time grid), alternating with the image.
Also the thermal noise from within-line differences in the stopped tail (the pass-to-pass differences there are motion)."""
import os, sys, json, numpy as np, finufft
from scipy.interpolate import CubicSpline
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY
G = DY.G; D = H.load(); k = D['k']; raw = D['raw']
# ---- thermal noise: differences along the readout inside the stopped tail (same k, same line) ----
tl = raw[:, H.NILV:, 462:512]
for lag in (1, 2, 4, 8):
    d = (tl[:, :, lag:] - tl[:, :, :-lag]).reshape(2, -1) / np.sqrt(2)
    C = d @ d.conj().T / d.shape[1]
    print(f'tail within-line difference, lag {lag}: var ch1 {C[0,0].real:.3e} ch2 {C[1,1].real:.3e} |corr| {abs(C[0,1])/np.sqrt(C[0,0].real*C[1,1].real):.3f}')
print('(pass-to-pass estimate used so far: 1.40e-10, 1.70e-11, 0.30)')
# ---- receiver filter ----
OS, PAD, J = 4, 128, 24
s0 = np.arange(H.NPTS); sf = np.arange(0, H.NPTS - 1 + 1e-9, 1 / OS); nf = sf.size
kf = CubicSpline(s0, k, axis=1)(sf); kap = (2 * np.pi * H.DX * kf.reshape(-1, 3))
pf = finufft.Plan(2, G, eps=1e-5, isign=-1); pa = finufft.Plan(1, G, eps=1e-5, isign=+1)
for p in (pf, pa): p.setpts(kap[:, 0].copy(), kap[:, 1].copy(), kap[:, 2].copy())
ff = np.fft.fftfreq(nf + 2 * PAD, 1 / OS); keep = np.arange(H.KILL, H.NPTS) * OS
def make(Hf):
    def filt(z):
        zp = np.pad(z, ((0, 0), (PAD, PAD)), mode='edge')
        return np.fft.ifft(np.fft.fft(zp, axis=1) * Hf[None, :], axis=1)[:, PAD:PAD + nf]
    def filt_adj(z):
        zp = np.zeros((z.shape[0], nf + 2 * PAD), complex); zp[:, PAD:PAD + nf] = z
        zp = np.fft.ifft(np.fft.fft(zp, axis=1) * np.conj(Hf)[None, :], axis=1)
        out = zp[:, PAD:PAD + nf].copy(); out[:, 0] += zp[:, :PAD].sum(1); out[:, -1] += zp[:, PAD + nf:].sum(1)
        return out
    fwd = lambda x: filt(pf.execute(x).reshape(H.NILV, nf))[:, keep]
    def adj(y):
        z = np.zeros((H.NILV, nf), complex); z[:, keep] = y
        return pa.execute(filt_adj(z).ravel())
    return fwd, adj
ybar = np.load(os.path.join(H.EXT, 'ybar_train.npy')).astype(np.complex128)[0]
dcf = np.load(os.path.join(H.EXT, 'dcf_g.npy')); Wt = dcf / dcf.sum(); M = np.load(os.path.join(H.EXT, 'support_g.npy'))
bands = ((8, 48), (48, 198), (198, 510)); env = np.sqrt((np.abs(ybar) ** 2).mean(0))
def fit(fwd, adj, it=40, x0=None):
    rhs = M * adj(Wt * ybar); n0 = lambda z: M * adj(Wt * fwd(z)); sc = np.vdot(rhs, n0(rhs)).real / np.vdot(rhs, rhs).real
    x = H.cg(lambda z: n0(z) + 1e-4 * sc * z, rhs, it=it, x0=x0)
    return x, [(np.abs((fwd(x) - ybar)[:, a_:b_]) ** 2).sum() / (np.abs(ybar[:, a_:b_]) ** 2).sum() for a_, b_ in bands]
Hf = np.clip((0.5 - np.abs(ff)) / 0.1, 0, 1); Hf = 0.5 - 0.5 * np.cos(np.pi * Hf)
x = None; lags = np.arange(-J, J + 1)
for rnd in range(4):
    fwd, adj = make(Hf)
    x, e = fit(fwd, adj, x0=x)
    print(f'round {rnd}: misfit/signal power {np.round(e, 4)}', flush=True)
    yfine = np.pad(pf.execute(x).reshape(H.NILV, nf), ((0, 0), (PAD, PAD)), mode='edge')
    sl = slice(8, 440)                                            # fit taps where the trajectory moves
    cols = np.stack([yfine[:, PAD + keep[sl] - j] for j in lags], axis=-1) / env[None, sl, None]      # (832, n, taps)
    A_ = cols.reshape(-1, lags.size); b_ = (ybar[:, sl] / env[None, sl]).ravel()
    h, *_ = np.linalg.lstsq(A_, b_, rcond=None)
    pred = A_ @ h; print(f'   FIR fit: relative residual power {np.vdot(b_ - pred, b_ - pred).real / np.vdot(b_, b_).real:.4f}; sum of taps {h.sum():.3f}')
    hz = np.zeros(nf + 2 * PAD, complex); hz[(lags) % (nf + 2 * PAD)] = h          # y[n] = sum_j h_j yfine[n - j]
    Hf = np.fft.fft(hz)
fr = np.fft.fftshift(ff); Hs = np.fft.fftshift(Hf)
sel = [np.argmin(np.abs(fr - v)) for v in (0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.6, 0.8, 1.0, 1.5)]
print('estimated |H| at f = 0, .1, .2, .3, .4, .45, .5, .6, .8, 1.0, 1.5 cycles/sample:', np.round(np.abs(Hs[sel]), 3))
print('phase (deg):', np.round(np.degrees(np.angle(Hs[sel])), 1), ' and at -f:', np.round(np.degrees(np.angle(Hs[[np.argmin(np.abs(fr + v)) for v in (0.1, 0.2, 0.3, 0.4)]])), 1))
np.savez(os.path.join(H.EXT, 'receiver_fir.npz'), h=h, lags=lags, OS=OS, Hf=Hf, ff=ff)
np.save(os.path.join(H.EXT, 'static_g_ls_ch0_rx.npy'), x.astype(np.complex64))
fig, ax = plt.subplots(1, 3, figsize=(17, 4.5))
ax[0].plot(fr, np.abs(Hs)); ax[0].set_xlim(-2, 2); ax[0].axvline(0.5, c='r', lw=.5); ax[0].axvline(-0.5, c='r', lw=.5); ax[0].set_title('estimated receiver response |H(f)|'); ax[0].set_xlabel('cycles per 5 us sample')
ax[1].plot(fr, np.degrees(np.angle(Hs))); ax[1].set_xlim(-0.6, 0.6); ax[1].set_title('phase (deg)')
ax[2].stem(lags / OS, np.abs(h)); ax[2].set_title('|FIR taps| against lag (samples)')
fig.tight_layout(); fig.savefig(os.path.join(H.OUT, 'receiver_filter.png'), dpi=60)
