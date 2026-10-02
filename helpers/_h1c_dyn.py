"""Scratch (h1 challenge, ASAP = no sparsity prior): linear 4D least-squares engine on the full field.

Model: image at breathing coordinate c (cyclic, 0..1) = linear interpolation between nb node images x_b (b/nb), so
every line uses its own breathing state and nothing is blurred inside a bin. Forward operator = NUFFT on a 2x finer
time grid + receiver band-limit + decimation (H.RxNufft). Channels: the whitened virtual channels listed in `chs`
(default: only the principal one; channel 1 adds ~8 % SNR but its sensitivity ratio is not known well enough).
Body support, Tikhonov in space (lam_s) and on cyclic node-to-node differences (lam_t, quadratic = linear smoothing).
All passes replay the same 832 interleaves, so the normal operator lives on the unique set:
  N = sum_c S_c^H A^H [w * (Omega (A S_c x))],  Omega[b, b', interleave] = sum over lines of h_b h_b'.
DCF-weighted (w), solved by CG. Grid LR 176 x AP 120 x SI 160 at 3.5 mm; crop() gives the 100^3 Tyger frame.
"""
import os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H

G = (176, 120, 160)
KIT = '/Volumes/HoomHamExt/Work/Codes/2026_XeCS_Recon/h1challenge/shared/judging_kit_v2.npz'
BANDS = ((0, 48), (48, 198), (198, 510))


class Dyn:
    def __init__(self, nb=16, support=True, coord=None, chs=(0,), filt=True, data=None, k=None, vol=None):
        D = H.load()
        self.k = D['k'] if k is None else k
        self.v = H.prewhiten(D['raw'])[0] if data is None else data          # (2, lines, 512) whitened virtual channels
        self.c = H.cyc_coord(D['vol'] if vol is None else vol) if coord is None else coord
        self.L = self.v.shape[1]
        self.nb, self.chs = nb, tuple(chs)
        self.wpow = 1.0
        S = np.load(os.path.join(H.EXT, 'sens_g.npy')).astype(np.complex64)
        self.S = [S[c] for c in self.chs] if len(self.chs) > 1 else [None]
        self.M = np.load(os.path.join(H.EXT, 'support_g.npy')) if support else np.ones(G, bool)
        self.dcf = np.load(os.path.join(H.EXT, 'dcf_g.npy'))
        self.op = H.RxNufft(self.k, n=G, nt=nb, filt=filt)
        kit = np.load(KIT)
        self.kit = kit
        self.ss = np.arange(H.NILV, self.L)
        self.ho_ilv_lines = self.ss[np.isin(self.ss % H.NILV, kit['heldout_interleaves'])]

    def hat(self, lines):
        u = self.c[lines] * self.nb
        if getattr(self, 'hard', False):                   # hard assignment to the nearest node (no interpolation)
            return np.rint(u).astype(int) % self.nb, np.zeros(u.size)
        b0 = np.floor(u).astype(int) % self.nb
        return b0, u - np.floor(u)

    def assemble(self, lines, win=None):
        """Omega0[b, il] = sum h_b^2, Omega1[b, il] = sum h_b h_{b+1}; Y[c, b, il, s] = sum h_b y. win: readout weights (510,)."""
        nb = self.nb
        il = lines % H.NILV
        b0, f = self.hat(lines)
        b1 = (b0 + 1) % nb
        O0 = np.zeros((nb, H.NILV)); O1 = np.zeros((nb, H.NILV))
        np.add.at(O0, (b0, il), (1 - f) ** 2); np.add.at(O0, (b1, il), f ** 2)
        np.add.at(O1, (b0, il), (1 - f) * f)
        Y = np.zeros((len(self.chs), nb, H.NILV, H.NPTS - H.KILL), np.complex64)
        for i, c in enumerate(self.chs):
            y = self.v[c, lines, H.KILL:]
            np.add.at(Y[i], (b0, il), (1 - f)[:, None] * y)
            np.add.at(Y[i], (b1, il), f[:, None] * y)
        w = self.dcf ** self.wpow                         # wpow 1 = density-weighted, 0 = unweighted (maximum likelihood)
        w = w / w.sum()
        if win is not None:
            w = w * win[None, :]
        return O0.astype(np.float32), O1.astype(np.float32), Y, w.astype(np.float32)

    def fwd_all(self, x):
        return np.stack([self.op.fwd(x if s is None else s[None] * x) for s in self.S])

    def adj_all(self, T, w):
        out = None
        for i, s in enumerate(self.S):
            a = self.op.adj(w[None] * T[i])
            if s is not None:
                a *= np.conj(s)[None]
            out = a if out is None else out + a
        out *= self.M[None]
        return out

    def recon(self, lines, lam_s=1e-3, lam_t=0.0, it=25, win=None, x0=None, verbose=True, cb=None, moco=None, phase=None, prior=None):
        """moco = (list of sparse W_b, support voxel indices): the temporal penalty becomes |x_b - W_b x_{b+1}|^2,
        i.e. coupling along the breathing motion instead of at fixed positions.
        phase = (P unit-modulus (G) or (nb, G), mu): adds mu |Im(conj(P) x)|^2 (image = real amplitude x known smooth phase).
        All penalties are quadratic; CG runs on the real inner product, so the real-linear phase term is fine."""
        O0, O1, Y, w = self.assemble(lines, win)

        def gram(x):
            U = self.fwd_all(x)
            T = O0[None, :, :, None] * U + O1[None, :, :, None] * np.roll(U, -1, axis=1) + np.roll(O1[None, :, :, None] * U, 1, axis=1)
            return self.adj_all(T, w)

        rhs = self.adj_all(Y, w)
        g = gram(rhs)
        sc = np.vdot(rhs, g).real / np.vdot(rhs, rhs).real
        del g

        def normal(x):
            out = gram(x) + np.float32(lam_s * sc) * (x if prior is None else prior * x)   # prior = 1 / relative prior power, (G) or (nb, G)
            if phase is not None:                          # soft phase constraint: penalise the part of x out of phase with P
                P, mu = phase
                out += np.float32(mu * sc) * (P * (1j * (np.conj(P) * x).imag))
            if lam_t > 0 and moco is None:
                out += np.float32(lam_t * sc) * (2 * x - np.roll(x, 1, axis=0) - np.roll(x, -1, axis=0))
            elif lam_t > 0:
                Wm, idx = moco
                nb = self.nb
                xs = x.reshape(nb, -1)[:, idx]
                acc = np.zeros_like(xs)
                for b in range(nb):
                    if isinstance(Wm[b], tuple):               # symmetric half-shift pair
                        e = Wm[b][0] @ xs[b] - Wm[b][1] @ xs[(b + 1) % nb]
                        acc[b] += Wm[b][0].T @ e
                        acc[(b + 1) % nb] -= Wm[b][1].T @ e
                    else:
                        e = xs[b] - Wm[b] @ xs[(b + 1) % nb]
                        acc[b] += e
                        acc[(b + 1) % nb] -= Wm[b].T @ e
                o2 = out.reshape(nb, -1)
                o2[:, idx] += np.float32(lam_t * sc) * acc
            return out

        t0 = time.time()

        def _cb(i, x, r):
            if verbose and (i % 5 == 4 or i == 0):
                print(f'    cg {i+1:3d}: rel resid {r:.4f}  ({time.time()-t0:.0f} s)', flush=True)
            if cb is not None:
                cb(i, x)
        return H.cg(normal, rhs, x0=x0, it=it, cb=_cb)

    def predict(self, x, lines, U=None):
        """Model data (nch, len(lines), 510) at each line's own breathing coordinate."""
        if U is None:
            U = self.fwd_all(x)
        il = lines % H.NILV
        b0, f = self.hat(lines)
        b1 = (b0 + 1) % self.nb
        return (1 - f)[None, :, None] * U[:, b0, il] + f[None, :, None] * U[:, b1, il]

    def chi2(self, x, lines, U=None, bands=BANDS):
        """mean |prediction - data|^2 in thermal-noise units (1 = floor) per channel and readout band."""
        r = self.predict(x, lines, U) - self.v[list(self.chs)][:, lines, H.KILL:]
        return np.array([[(np.abs(r[i][:, a:b]) ** 2).mean() for a, b in bands] for i in range(len(self.chs))])


def crop(x):
    """(…, 176, 120, 160) local -> (…, 100, 100, 100) Tyger frame (SI, AP, LR)."""
    c = [(g - 100) // 2 for g in G]
    y = x[..., c[0]:c[0] + 100, c[1]:c[1] + 100, c[2]:c[2] + 100]
    return H.to_tyger(y) if y.ndim == 3 else np.stack([H.to_tyger(a) for a in y])
