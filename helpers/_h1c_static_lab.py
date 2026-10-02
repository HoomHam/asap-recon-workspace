"""Scratch: null-space regularisation lab on the STATIC problem (fast). Pass-averaged scored-train data, new model.
Score = held-out lines / held-out interleaves, thermal units per line, samples 2-50 / 50-200 / 200-512.
Knobs: CG iterations, Tikhonov, data-weight exponent, smoothing penalty on the body OUTSIDE the 350 mm box."""
import os, sys, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY
G = DY.G; BANDS = DY.BANDS
kit = np.load(DY.KIT); D = H.load(); k = D['k']; v = H.prewhiten(D['raw'])[0][0]
ybar = np.load(os.path.join(H.EXT, 'ybar_train.npy'))[0].astype(np.complex64); cnt = np.load(os.path.join(H.EXT, 'cnt_train.npy'))
dcf = np.load(os.path.join(H.EXT, 'dcf_g.npy')); M = np.load(os.path.join(H.EXT, 'support_g.npy'))
op = H.RxNufft(k, n=G, nt=1)
ss = np.arange(H.NILV, v.shape[0]); ho_l = kit['heldout_lines']; ho_i = ss[np.isin(ss % H.NILV, kit['heldout_interleaves'])]
inbox = np.zeros(G, bool); inbox[38:138, 10:110, 30:130] = True
wout = (M & ~inbox).astype(np.float32)

def grad_pen(x):
    out = np.zeros_like(x)
    for a in range(3):
        sl0 = [slice(None)] * 3; sl1 = [slice(None)] * 3; sl0[a] = slice(0, -1); sl1[a] = slice(1, None)
        d = (x[tuple(sl1)] - x[tuple(sl0)]) * np.minimum(wout[tuple(sl1)], wout[tuple(sl0)])
        out[tuple(sl1)] += d; out[tuple(sl0)] -= d
    return out

def score(x):
    U = op.fwd(x)
    r = []
    for h in (ho_l, ho_i):
        e = U[h % H.NILV] - v[h, H.KILL:]
        r.append([(np.abs(e[:, a:b]) ** 2).mean() for a, b in BANDS])
    return r

def fit(it=15, lam_s=1e-3, wpow=1.0, lam_o=0.0, tag='', mask=M, its_report=None):
    w = (dcf ** wpow) * cnt[:, None]; w = (w / w.sum()).astype(np.float32)
    rhs = mask * op.adj(w * ybar); n0 = lambda z: mask * op.adj(w * op.fwd(z))
    sc = np.vdot(rhs, n0(rhs)).real / np.vdot(rhs, rhs).real
    def nrm(z):
        o = n0(z) + np.float32(lam_s * sc) * z
        if lam_o > 0: o = o + np.float32(lam_o * sc) * (mask * grad_pen(z))
        return o
    def cb(i, z, r):
        if its_report and (i + 1) in its_report:
            s = score(z); print(f'{tag:46s} it {i+1:3d}: lines {s[0][0]:8.0f} {s[0][1]:6.1f} {s[0][2]:6.1f} | interleaves {s[1][0]:8.0f} {s[1][1]:6.1f} {s[1][2]:6.1f}', flush=True)
    x = H.cg(nrm, rhs, it=it, cb=cb)
    if not its_report:
        s = score(x); print(f'{tag:46s} it {it:3d}: lines {s[0][0]:8.0f} {s[0][1]:6.1f} {s[0][2]:6.1f} | interleaves {s[1][0]:8.0f} {s[1][1]:6.1f} {s[1][2]:6.1f}', flush=True)
    return x

if __name__ == '__main__':
    fit(it=60, tag='baseline (dcf weights, tikhonov 1e-3)', its_report=(5, 8, 12, 15, 20, 30, 60))
    for ls in (1e-2, 3e-2, 1e-1):
        fit(it=60, lam_s=ls, tag=f'tikhonov {ls:g}', its_report=(15, 60))
    for lo in (0.1, 1.0, 10.0):
        fit(it=60, lam_o=lo, tag=f'outside-box smoothing {lo:g}', its_report=(15, 60))
    for p in (0.75, 0.5):
        fit(it=60, wpow=p, tag=f'data-weight exponent {p}', its_report=(15, 30, 60))
