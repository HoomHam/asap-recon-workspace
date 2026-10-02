"""Scratch (2026-10-02, h1 challenge, ASAP side = no sparsity prior): toolkit for the 006KL 1H free-breathing cine.

Shared inputs (rule 1): XeCS cache on Ext (raw (2, 13000, 512), measured water trajectory (832, 512, 3) cycles/mm,
k0 surrogate vol, high = inspiration). Grid: 100^3, 3.5 mm, FOV 350, delivered as (bins, SI, AP, LR) like p/recon.mat.

load()            cache -> dict
prewhiten()       2-channel noise covariance from pass-to-pass differences in the readout tail -> whitened,
                  SVD-rotated virtual channels (channel 0 carries the signal)
cyc_coord()       Steve's bin(): rank-ordered cyclic breathing coordinate in [0, 1) (expiration->inspiration->back)
soft_w()          Gaussian cyclic bin weights (lines x bins), as cudarecon (exp(-d^2/s2))
pipe_dcf()        Pipe-Menon DCF on the 832 unique interleaves
Nufft             finufft plan pair (type 2 forward, type 1 adjoint) on a fixed set of lines, grid n^3 of voxel dx
to_tyger()        local (x, y, z) array -> (SI, AP, LR) orientation of p/recon.mat
"""
import os
import numpy as np
import finufft

HERE = os.path.dirname(os.path.abspath(__file__))
WS = os.path.dirname(HERE)
CACHE = '/Volumes/HoomHamExt/Work/Codes/2026_XeCS_Recon/h1challenge/cache_006KL_1H.npz'
EXT = '/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/h1challenge'
OUT = os.path.join(WS, 'outputs', 'h1challenge_2026-10-02')
TYGER_P = '/Volumes/HoomHamExt/AIkill_Dynamic_1H/2024-03-12_006KL/p/recon.mat'
NILV, NPTS, KILL, DWELL = 832, 512, 2, 5e-6
N, DX = 100, 3.5


def load():
    z = np.load(CACHE)
    raw = z['raw']
    tr = float(z['tr'])
    return dict(raw=raw, k=z['k'], tr=tr, t=np.arange(raw.shape[1]) * tr, vol=z['vol'].astype(np.float64))


def prewhiten(raw, tail=slice(462, 512), lag=2):
    """Thermal noise covariance from differences ALONG the readout inside the stopped tail of each line (same k, same
    line; the pass-to-pass differences there are breathing, not noise), whiten, then rotate so that virtual channel 0
    holds the common signal."""
    tl = raw[:, NILV:, tail]
    d = (tl[:, :, lag:] - tl[:, :, :-lag]).reshape(2, -1) / np.sqrt(2)
    C = d @ d.conj().T / d.shape[1]
    Lc = np.linalg.cholesky(C)
    W = np.linalg.inv(Lc)
    rw = np.einsum('ab,bls->als', W, raw)
    s = rw[:, NILV:, KILL:64].reshape(2, -1)                      # signal-dominated samples, steady state
    S = s @ s.conj().T / s.shape[1]
    ev, U = np.linalg.eigh(S)
    U = U[:, ::-1]
    v = np.einsum('ba,bls->als', U.conj(), rw).astype(np.complex64)
    return v, dict(noise_cov=C, signal_eig=ev[::-1], U=U, W=W)


def cyc_coord(vol, lines=None):
    """Steve's raw.bin(): falling samples are mirrored above the rising ones, then rank / count."""
    v = np.asarray(vol, float).copy()
    v = (v - v.min()) / (v.max() - v.min())
    fall = np.zeros(v.size, bool)
    fall[1:] = v[1:] < v[:-1]
    fall[0] = v[1] < v[0]
    maxf = v[fall].max() * 1.00001
    maxr = v[~fall].max()
    u = v.copy()
    u[fall] = maxf + maxr - v[fall]
    c = np.empty(v.size)
    if lines is None:
        lines = np.arange(v.size)
    order = np.sort(u[lines])
    c = np.searchsorted(order, u) / lines.size
    return c.clip(0, 1 - 1e-9)


def soft_w(c, nbins=16, s2=0.44):
    """(lines, bins) weights exp(-d^2 / s2), d = cyclic distance in bins between line coordinate and bin index."""
    x = c[:, None] * nbins
    b = np.arange(nbins)[None, :]
    d = np.abs(x - b)
    d = np.minimum(d, nbins - d)
    return np.exp(-d * d / s2)


def pipe_dcf(k, n=N, dx=DX, it=30):
    import sigpy.mri as smr
    w = smr.pipe_menon_dcf(k.reshape(-1, 3) * n * dx, img_shape=(n, n, n), max_iter=it, show_pbar=False)
    return np.abs(w).reshape(k.shape[:2])


class Nufft:
    """A x: image (n1, n2, n3) -> samples (L, S) on lines `lines` (samples KILL..), and its adjoint.
    nt > 1: stacks of nt images / sample sets in one call (single precision by default for the 4D engine)."""

    def __init__(self, k, lines, n=N, dx=DX, eps=1e-5, nt=1, single=False):
        kk = k[np.asarray(lines) % NILV][:, KILL:]
        self.shape = kk.shape[:2]
        self.nt = nt
        ft = np.float32 if single else np.float64
        self.ct = 'complex64' if single else 'complex128'
        kx = np.ascontiguousarray(2 * np.pi * dx * kk.reshape(-1, 3)).astype(ft)
        self.n = (n, n, n) if np.isscalar(n) else tuple(n)
        self.f = finufft.Plan(2, self.n, n_trans=nt, eps=eps, isign=-1, dtype=self.ct)
        self.a = finufft.Plan(1, self.n, n_trans=nt, eps=eps, isign=+1, dtype=self.ct)
        xs = [np.ascontiguousarray(kx[:, i]) for i in range(3)]
        self.f.setpts(*xs)
        self.a.setpts(*xs)

    def fwd(self, x):
        y = self.f.execute(np.ascontiguousarray(x, dtype=self.ct))
        return y.reshape(self.shape) if self.nt == 1 else y.reshape((self.nt,) + self.shape)

    def adj(self, y):
        y = np.ascontiguousarray(y, dtype=self.ct)
        return self.a.execute(y.ravel() if self.nt == 1 else y.reshape(self.nt, -1))


class RxNufft:
    """Forward model with the receiver band-limit: y = decimate( h * (A on an os-times finer time grid) ).
    The trajectory moves ~1.2/350 mm^-1 per sample, so the receiver band covers only +-150 mm along the gradient;
    h (estimated from the data, shared/receiver_filter_v1.npz) removes what the scanner's decimation filter removed.
    Same interface as Nufft: fwd (nt, n1, n2, n3) -> (nt, 832, 510), adj the reverse."""
    PAD = 64

    def __init__(self, k, n=N, dx=DX, eps=1e-4, nt=1, single=True, os_=2, filt=True):
        from scipy.interpolate import CubicSpline
        self.nt, self.os = nt, os_
        sf = np.arange(0, NPTS - 1 + 1e-9, 1 / os_)
        self.nf = sf.size
        kf = CubicSpline(np.arange(NPTS), k, axis=1)(sf)
        ft = np.float32 if single else np.float64
        self.ct = 'complex64' if single else 'complex128'
        kx = np.ascontiguousarray(2 * np.pi * dx * kf.reshape(-1, 3)).astype(ft)
        self.n = (n, n, n) if np.isscalar(n) else tuple(n)
        self.f = finufft.Plan(2, self.n, n_trans=nt, eps=eps, isign=-1, dtype=self.ct)
        self.a = finufft.Plan(1, self.n, n_trans=nt, eps=eps, isign=+1, dtype=self.ct)
        xs = [np.ascontiguousarray(kx[:, i]) for i in range(3)]
        self.f.setpts(*xs)
        self.a.setpts(*xs)
        self.keep = np.arange(KILL, NPTS) * os_
        self.shape = (NILV, NPTS - KILL)
        ff = np.fft.fftfreq(self.nf + 2 * self.PAD, 1 / os_)
        if filt:
            z = np.load(os.path.join(EXT, 'shared', 'receiver_filter_v1.npz'))
            self.Hf = (np.interp(ff, z['f'], z['H'].real) + 1j * np.interp(ff, z['f'], z['H'].imag)).astype(self.ct)
        else:
            self.Hf = np.ones(ff.size, self.ct)

    def fwd(self, x):
        import scipy.fft as sfft
        P = self.PAD
        y = self.f.execute(np.ascontiguousarray(x, dtype=self.ct)).reshape(-1, NILV, self.nf)
        yp = np.pad(y, ((0, 0), (0, 0), (P, P)), mode='edge')
        yp = sfft.ifft(sfft.fft(yp, axis=-1, workers=-1) * self.Hf, axis=-1, workers=-1)
        out = yp[:, :, P + self.keep]
        return out[0] if self.nt == 1 else out

    def adj(self, y):
        import scipy.fft as sfft
        P = self.PAD
        y = np.asarray(y, dtype=self.ct).reshape(-1, NILV, NPTS - KILL)
        zp = np.zeros((y.shape[0], NILV, self.nf + 2 * P), self.ct)
        zp[:, :, P + self.keep] = y
        zp = sfft.ifft(sfft.fft(zp, axis=-1, workers=-1) * np.conj(self.Hf), axis=-1, workers=-1)
        z = np.ascontiguousarray(zp[:, :, P:P + self.nf])
        z[:, :, 0] += zp[:, :, :P].sum(-1)
        z[:, :, -1] += zp[:, :, P + self.nf:].sum(-1)
        return self.a.execute(z.reshape(self.nt, -1) if self.nt > 1 else z.ravel())


_T = None


def to_tyger(x):
    """local (x, y, z) -> (SI, AP, LR) of p/recon.mat. Transform found once by _h1c_orient.py (perm, flips)."""
    global _T
    if _T is None:
        _T = np.load(os.path.join(EXT, 'orient.npz'))
    y = np.transpose(x, tuple(_T['perm']) if x.ndim == 3 else (0,) + tuple(1 + _T['perm']))
    for a, f in enumerate(_T['flip']):
        if f:                                              # flipping an even grid moves the centre by one voxel
            y = np.roll(np.flip(y, axis=a + (x.ndim - 3)), 1, axis=a + (x.ndim - 3))
    return np.ascontiguousarray(y)


def cg(Aop, b, x0=None, it=20, tol=1e-6, cb=None):
    """Conjugate gradients on the Hermitian positive system Aop(x) = b."""
    x = np.zeros_like(b) if x0 is None else x0.copy()
    r = b - Aop(x)
    p = r.copy()
    rs = np.vdot(r, r).real
    b2 = np.vdot(b, b).real
    for i in range(it):
        Ap = Aop(p)
        al = rs / np.vdot(p, Ap).real
        x += al * p
        r -= al * Ap
        rn = np.vdot(r, r).real
        if cb is not None:
            cb(i, x, np.sqrt(rn / b2))
        if rn < tol * tol * b2:
            break
        p = r + (rn / rs) * p
        rs = rn
    return x
