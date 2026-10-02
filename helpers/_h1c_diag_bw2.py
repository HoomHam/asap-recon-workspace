"""Scratch: does a receiver band-limit in the forward model remove the fixed misfit?
Model: y = decimate( h * (A_fine x) ), A_fine on a 4x finer time grid (trajectory spline), h = low-pass along time.
Static, channel 0, body support, pass-averaged data, 40 CG iterations; misfit / signal power per band."""
import os, sys, numpy as np, finufft
from scipy.interpolate import CubicSpline
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY
G = DY.G; D = H.load(); k = D['k']; raw = D['raw']
# receiver filter shape from noise: pass-to-pass differences in the stopped tail (samples 448..511)
L = raw.shape[1]; npass = L // H.NILV
a = raw[0, :npass * H.NILV].reshape(npass, H.NILV, H.NPTS)[1:, :, 448:512]
nz = (a[1:] - a[:-1]).reshape(-1, 64)
Pn = (np.abs(np.fft.fftshift(np.fft.fft(nz, axis=1), axes=1)) ** 2).mean(0); fn = np.fft.fftshift(np.fft.fftfreq(64))
print('tail noise PSD (normalised to max) at |f| = 0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.48:', np.round([Pn[np.argmin(np.abs(fn - v))] / Pn.max() for v in (0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.484)], 3))
ac = [np.real((nz[:, :-l] * np.conj(nz[:, l:])).mean()) / np.real((np.abs(nz) ** 2).mean()) for l in (1, 2, 3)]
print('tail noise autocorrelation at lags 1, 2, 3:', np.round(ac, 3))
OS, PAD = 4, 128
s0 = np.arange(H.NPTS); sf = np.arange(0, H.NPTS - 1 + 1e-9, 1 / OS); nf = sf.size
kf = CubicSpline(s0, k, axis=1)(sf); kap = (2 * np.pi * H.DX * kf.reshape(-1, 3))
pf = finufft.Plan(2, G, eps=1e-5, isign=-1); pa = finufft.Plan(1, G, eps=1e-5, isign=+1)
for p in (pf, pa): p.setpts(kap[:, 0].copy(), kap[:, 1].copy(), kap[:, 2].copy())
ff = np.fft.fftfreq(nf + 2 * PAD, 1 / OS)                                  # cycles per 5 us sample
keep = np.arange(H.KILL, H.NPTS) * OS
def make(Hf):
    def filt(z):
        zp = np.pad(z, ((0, 0), (PAD, PAD)), mode='edge')
        return np.fft.ifft(np.fft.fft(zp, axis=1) * Hf[None, :], axis=1)[:, PAD:PAD + nf]
    def filt_adj(z):                                                        # adjoint of pad(edge) + filter + crop
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
bands = ((8, 48), (48, 198), (198, 510))
def fit(fwd, adj, tag, its=(15, 40)):
    rhs = M * adj(Wt * ybar); n0 = lambda z: M * adj(Wt * fwd(z)); sc = np.vdot(rhs, n0(rhs)).real / np.vdot(rhs, rhs).real
    def cb(i, z, r):
        if (i + 1) in its:
            e = [(np.abs((fwd(z) - ybar)[:, a_:b_]) ** 2).sum() / (np.abs(ybar[:, a_:b_]) ** 2).sum() for a_, b_ in bands]
            print(f'{tag:44s} it {i+1:3d}: misfit/signal power {np.round(e, 4)}', flush=True)
    return H.cg(lambda z: n0(z) + 1e-4 * sc * z, rhs, it=max(its), cb=cb)
# adjoint check
rng = np.random.default_rng(0); fw, ad = make((np.abs(ff) <= 0.5).astype(float))
u = rng.normal(size=G) + 1j * rng.normal(size=G); w = rng.normal(size=(H.NILV, 510)) + 1j * rng.normal(size=(H.NILV, 510))
print('adjoint check (should be ~0):', abs(np.vdot(fw(u), w) - np.vdot(u, ad(w))) / abs(np.vdot(fw(u), w)))
fit(*make(np.ones_like(ff)), 'no filter (reference, 4x time grid)', its=(40,))
x = fit(fw, ad, 'brick wall at +-0.5 cycles/sample')
np.save(os.path.join(H.EXT, 'static_g_ls_ch0_bw.npy'), x.astype(np.complex64))
for fc, wd in ((0.45, 0.05), (0.40, 0.10), (0.30, 0.20)):
    Hc = np.clip((fc + wd - np.abs(ff)) / (2 * wd), 0, 1); Hc = 0.5 - 0.5 * np.cos(np.pi * Hc)
    fit(*make(Hc), f'raised-cosine edge, half-amplitude at {fc}, +-{wd}', its=(40,))
