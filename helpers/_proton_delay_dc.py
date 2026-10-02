"""Scratch: gradient-delay estimate by DATA CONSISTENCY (006KL 1H, end-exp lines).
For each tau: weighted-LS low-res fit (CG on A^H W A x = A^H W y, per coil) with |k|<=KMAX,
residual ||W^1/2 (A x - y)|| / ||W^1/2 y||.  Minimum = best delay. Then per-axis coordinate refine."""
import sys, os, json, time
import numpy as np, finufft
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _proton_local as L
KMAX, N, DX, IT = 0.04, 64, 6.0, 20
raw, k, tr = L.load()
lines, _ = L.endexp_lines(raw.shape[2])
lines = lines[::2]                                          # half the end-exp lines: plenty, faster
dcf = L.pipe_dcf(k, N, DX)
il = lines % 832
Y = [raw[L.KILL:, c, :][:, lines].T for c in range(raw.shape[1])]   # (Lines, S)

def resid(tau):
    kk = L.delay_traj(k, tau)[il][:, L.KILL:]
    keep = np.linalg.norm(kk - kk[:, :1], axis=-1) <= KMAX          # select by true radius (offset-free)
    kx = (2 * np.pi * DX * kk[keep]).astype(np.float64)
    w = dcf[il][:, L.KILL:][keep]
    A = lambda x: finufft.nufft3d2(kx[:, 0], kx[:, 1], kx[:, 2], x, eps=1e-6)
    AH = lambda y: finufft.nufft3d1(kx[:, 0], kx[:, 1], kx[:, 2], y, (N, N, N), eps=1e-6)
    num = den = 0.0
    for y in Y:
        y = y[keep].astype(np.complex128)
        b = AH(w * y); x = np.zeros((N, N, N), complex); r = b.copy(); p = r.copy(); rr = np.vdot(r, r).real
        for _ in range(IT):
            Ap = AH(w * A(p)); a = rr / np.vdot(p, Ap).real; x += a * p; r -= a * Ap
            rn = np.vdot(r, r).real; p = r + (rn / rr) * p; rr = rn
        e = A(x) - y
        num += (w * np.abs(e) ** 2).sum(); den += (w * np.abs(y) ** 2).sum()
    return num / den

out = sys.argv[1]; res = {}
t0 = time.time()
for tau in np.arange(-3, 3.01, 0.5):
    res[f'iso {tau:+.2f}'] = r = resid(tau); print(f'iso tau {tau:+.2f} ({tau*5:+.1f} us): resid {r:.5f}  [{time.time()-t0:.0f}s]', flush=True)
iso = min((v, float(kk_.split()[1])) for kk_, v in res.items())[1]
best = np.array([iso] * 3)
for sweep in range(2):
    for a in range(3):
        cand = {}
        for d in (-1, -0.5, -0.25, 0, 0.25, 0.5, 1):
            t = best.copy(); t[a] += d; cand[d] = resid(t)
        dbest = min(cand, key=cand.get); best[a] += dbest
        print(f'axis {a}: tau -> {best.tolist()}  resid {cand[dbest]:.5f}  (scan {", ".join(f"{d:+}:{v:.5f}" for d, v in cand.items())})', flush=True)
res['best_axes'] = best.tolist(); res['resid_best'] = resid(best); res['resid_zero'] = resid(0.0)
print(f'BEST tau per axis (samples) {best.tolist()} = {(best*5).tolist()} us; resid {res["resid_best"]:.5f} vs tau=0 {res["resid_zero"]:.5f}', flush=True)
json.dump(res, open(os.path.join(out, 'delay_dc.json'), 'w'), indent=1)
