"""Scratch: does the second channel shrink the null space once its sensitivity is fitted to the data?
Static bench (new model, pass-averaged scored-train data, both whitened virtual channels).
Channel 0 sensitivity = 1 by convention, channel 1 = R(r), R on a coarse trilinear basis (14 mm knots).
Alternate: image m from both channels (CG) <-> R from channel 1 given m (CG on the coarse coefficients).
Score: held-out lines / interleaves, thermal units, samples 2-50 / 50-200 / 200-512, each channel."""
import os, sys, time, numpy as np
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY
G = DY.G; BANDS = DY.BANDS
kit = np.load(DY.KIT); D = H.load(); k = D['k']; v = H.prewhiten(D['raw'])[0]
ybar = np.load(os.path.join(H.EXT, 'ybar_train.npy')).astype(np.complex64); cnt = np.load(os.path.join(H.EXT, 'cnt_train.npy'))
dcf = np.load(os.path.join(H.EXT, 'dcf_g.npy')); M = np.load(os.path.join(H.EXT, 'support_g.npy'))
op = H.RxNufft(k, n=G, nt=1)
ss = np.arange(H.NILV, v.shape[1]); ho_l = kit['heldout_lines']; ho_i = ss[np.isin(ss % H.NILV, kit['heldout_interleaves'])]
WP = float(os.environ.get('WPOW', 0.75)); LS = 0.003
w = (dcf ** WP) * cnt[:, None]; w = (w / w.sum()).astype(np.float32)
STEP = 4
def interp_mat(n, step):
    nc = n // step + 2; pos = np.arange(n) / step; i0 = np.floor(pos).astype(int); f = pos - i0
    return sp.csr_matrix((np.concatenate([1 - f, f]), (np.concatenate([np.arange(n)] * 2), np.concatenate([i0, i0 + 1]))), shape=(n, nc)).toarray().astype(np.float32)
B = [interp_mat(n, STEP) for n in G]
up = lambda c: np.einsum('ia,jb,kc,abc->ijk', B[0], B[1], B[2], c, optimize=True)
dn = lambda x: np.einsum('ia,jb,kc,ijk->abc', B[0], B[1], B[2], x, optimize=True)

def score(m, R, tag):
    out = []
    for c, s in ((0, None), (1, R)):
        U = op.fwd(m if s is None else s * m)
        for h in (ho_l, ho_i):
            e = U[h % H.NILV] - v[c, h, H.KILL:]
            out.append([(np.abs(e[:, a:b]) ** 2).mean() for a, b in BANDS])
    print(f'{tag:40s} ch0 lines {out[0][0]:7.0f} {out[0][1]:6.1f} {out[0][2]:5.1f} | ch0 ilv {out[1][0]:7.0f} {out[1][1]:6.1f} {out[1][2]:5.1f} || ch1 lines {out[2][0]:7.0f} {out[2][1]:6.1f} {out[2][2]:5.1f} | ch1 ilv {out[3][0]:7.0f} {out[3][1]:6.1f} {out[3][2]:5.1f}', flush=True)

def fit_m(R, it=20, w1=1.0):
    Ss = [None] if R is None else [None, R]; wc = [1.0, w1]
    def n0(z):
        o = op.adj(w * op.fwd(z))
        if R is not None: o = o + np.float32(w1) * np.conj(R) * op.adj(w * op.fwd(R * z))
        return M * o
    rhs = op.adj(w * ybar[0])
    if R is not None: rhs = rhs + np.float32(w1) * np.conj(R) * op.adj(w * ybar[1])
    rhs = M * rhs; sc = np.vdot(rhs, n0(rhs)).real / np.vdot(rhs, rhs).real
    return H.cg(lambda z: n0(z) + np.float32(LS * sc) * z, rhs, it=it)

def fit_R(m, R0, it=12, lam=1e-2):
    """R = R0 + up(c): least squares for the coarse correction c."""
    res = ybar[1] - op.fwd(R0 * m)
    A_ = lambda c: op.fwd(m * up(c).astype(np.complex64)); At = lambda y: dn(np.conj(m) * op.adj(w * y))
    rhs = At(res); nrm0 = lambda c: At(A_(c)); sc = np.vdot(rhs, nrm0(rhs)).real / np.vdot(rhs, rhs).real
    c = H.cg(lambda c: nrm0(c) + lam * sc * c, rhs, it=it)
    return (R0 + up(c)).astype(np.complex64)

S = np.load(os.path.join(H.EXT, 'sens_g.npy')); R = (S[1] / S[0]).astype(np.complex64)
m = fit_m(None); score(m, R, 'channel 0 only (R = low-res ratio)')
for w1 in (1.0,):
    Rk = R
    for rnd in range(4):
        m = fit_m(Rk, w1=w1); score(m, Rk, f'two channels, round {rnd} (ch1 weight {w1})')
        Rk = fit_R(m, Rk)
    np.save(os.path.join(H.EXT, 'R_fitted.npy'), Rk)
    print('fitted |R| in support: median %.3f, 5-95 pct %.3f-%.3f ; rms change from the low-res ratio %.3f' % (*np.percentile(np.abs(Rk)[M], [50, 5, 95]), np.sqrt((np.abs(Rk - R)[M] ** 2).mean())))
