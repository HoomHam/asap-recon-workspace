"""Scratch: trajectory scale/delay by CROSS-VALIDATED data consistency (006KL 1H, end-exp lines).
Fit low-res image on even passes (DCF-weighted CG), predict odd-pass samples; error vs scale / tau."""
import sys, os, json, numpy as np, finufft
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _proton_local as L
KMAX, N, DX, IT = 0.04, 64, 6.0, 15
raw, k, tr = L.load(); lines, _ = L.endexp_lines(raw.shape[2])
dcf = L.pipe_dcf(k, N, DX)
tr_l, te_l = lines[(lines % 832) % 2 == 0], lines[(lines % 832) % 2 == 1]   # split by INTERLEAVE: test k-locations differ from train
def setup(lns):
    il = lns % 832; kk = k[il][:, L.KILL:]; keep = np.linalg.norm(kk - kk[:, :1], axis=-1) <= KMAX
    return il, keep, dcf[il][:, L.KILL:][keep], [raw[L.KILL:, c, :][:, lns].T[keep].astype(np.complex128) for c in range(raw.shape[1])]
TR, TE = setup(tr_l), setup(te_l)
def geom(il, keep, s, tau):
    kk = L.delay_traj(k, tau)[il][:, L.KILL:] if np.any(np.asarray(tau) != 0) else k[il][:, L.KILL:]
    kk = (kk - kk[:, :1]) * s + kk[:, :1]
    return (2 * np.pi * DX * kk[keep]).astype(np.float64)
def cv(s=1.0, tau=0.0):
    kx, kt = geom(TR[0], TR[1], s, tau), geom(TE[0], TE[1], s, tau)
    A = lambda kx, x: finufft.nufft3d2(kx[:, 0], kx[:, 1], kx[:, 2], x, eps=1e-6)
    AH = lambda kx, y: finufft.nufft3d1(kx[:, 0], kx[:, 1], kx[:, 2], y, (N, N, N), eps=1e-6)
    w = TR[2]; num = den = 0.0
    for y, yt in zip(TR[3], TE[3]):
        b = AH(kx, w * y); x = np.zeros((N, N, N), complex); r = b.copy(); p = r.copy(); rr = np.vdot(r, r).real
        for _ in range(IT):
            Ap = AH(kx, w * A(kx, p)); a = rr / np.vdot(p, Ap).real; x += a * p; r -= a * Ap; rn = np.vdot(r, r).real; p = r + (rn / rr) * p; rr = rn
        e = A(kt, x) - yt; num += (TE[2] * np.abs(e) ** 2).sum(); den += (TE[2] * np.abs(yt) ** 2).sum()
    return num / den
res = {}
for s in [0.90, 0.94, 0.97, 1.0, 1.03, 1.06, 1.10]:
    res[f's{s}'] = v = cv(s=s); print(f'scale {s:.2f}: held-out err {v:.5f}', flush=True)
for t in [-1.0, -0.5, 0.0, 0.5, 1.0]:
    res[f't{t}'] = v = cv(tau=t); print(f'tau {t:+.2f}: held-out err {v:.5f}', flush=True)
json.dump(res, open(os.path.join(sys.argv[1], 'traj_cv_ilvsplit.json'), 'w'), indent=1)
