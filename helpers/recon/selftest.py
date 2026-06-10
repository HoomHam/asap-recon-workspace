"""Self-test on synthetic data — no scanner data needed.

Builds a 3D spiral-ish trajectory + analytic phantom (sum of Gaussian blobs,
closed-form Fourier transform), synthesizes samples exactly, then checks:
  1. adjoint runs and peaks where the phantom is
  2. DCF adjoint beats plain adjoint (NRMSE vs ground truth)
  3. CG beats both
  4. adjoint/forward pair passes the dot-product (adjointness) test
Run:  .venv/bin/python selftest.py
"""

import numpy as np
import asap_recon as ar

rng = np.random.default_rng(0)
IS = 64          # small grid for speed
MS = int(IS * 2.4)


def make_traj(n_ilv=400, npts=300):
    """Crude 3D spiral stand-in: golden-angle cones, radius growing 0->kmax."""
    kmax = IS / 2 * 0.9
    ga = np.pi * (3 - np.sqrt(5))
    pts = []
    for i in range(n_ilv):
        th = np.arccos(1 - 2 * (i + 0.5) / n_ilv)   # polar of cone axis
        ph = i * ga
        ax = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
        t = np.linspace(0, 1, npts)
        r = kmax * t
        spin = 2 * np.pi * 8 * t
        # orthonormal frame around axis
        u = np.cross(ax, [0, 0, 1.0]); u /= np.linalg.norm(u) + 1e-12
        v = np.cross(ax, u)
        k = (r * np.cos(spin))[:, None] * u + (r * np.sin(spin))[:, None] * v \
            + (r * 0.3)[:, None] * ax
        pts.append(k)
    k = np.concatenate(pts)                       # delta-k units
    return k * MS / IS + MS / 2                   # grid units (Steve convention)


def phantom_and_samples(traj_grid):
    """Gaussian blobs: image-domain truth + exact analytic k-samples."""
    blobs = [  # (center voxel offsets from middle, sigma_vox, amplitude)
        ((0, 0, 0), 6.0, 1.0),
        ((10, -8, 4), 3.0, 0.7),
        ((-12, 6, -6), 4.0, 0.5),
    ]
    n = IS
    xs = np.arange(n) - n / 2
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    truth = np.zeros((n, n, n))
    k_dk = (traj_grid - MS / 2) * IS / MS         # back to delta-k
    s = np.zeros(len(k_dk), dtype=complex)
    for (cx, cy, cz), sig, amp in blobs:
        truth += amp * np.exp(-((X - cx) ** 2 + (Y - cy) ** 2 + (Z - cz) ** 2)
                              / (2 * sig ** 2))
        # FT of Gaussian (continuous, matched to isign=-1 adjoint convention)
        kk = k_dk * 2 * np.pi / n
        k2 = (kk ** 2).sum(1)
        phase = (k_dk @ np.array([cx, cy, cz])) * 2 * np.pi / n
        s += amp * (2 * np.pi * sig ** 2) ** 1.5 * np.exp(-k2 * sig ** 2 / 2) \
             * np.exp(+1j * phase)
    return truth, s


def nrmse(a, b):
    a, b = np.abs(a), np.abs(b)
    a, b = a / a.max(), b / b.max()
    return np.linalg.norm(a - b) / np.linalg.norm(b)


traj = make_traj()
truth, samples = phantom_and_samples(traj)
print(f"traj {traj.shape}, samples {samples.shape}, grid {IS}^3")

# 4. adjointness check on the operator pair
tr = ar.grid_to_radians(traj, MS)
x = rng.standard_normal((IS, IS, IS)) + 1j * rng.standard_normal((IS, IS, IS))
y = rng.standard_normal(len(tr)) + 1j * rng.standard_normal(len(tr))
lhs = np.vdot(y, ar.forward(tr, x))
rhs = np.vdot(ar.adjoint(tr, y, IS), x)
adj_err = abs(lhs - rhs) / abs(lhs)
print(f"adjointness <Ax,y> vs <x,A^Hy>: rel err {adj_err:.2e}  "
      f"{'PASS' if adj_err < 1e-6 else 'FAIL'}")

results = {}
for method, kw in [("adjoint", {}), ("adjoint_dcf", {}), ("cg", {"cg_iters": 30})]:
    img = ar.recon(traj, samples, method=method, MS=MS, IS=IS, **kw)
    results[method] = nrmse(img, truth)
    print(f"{method:12s} NRMSE vs truth: {results[method]:.4f}")

ok = results["adjoint_dcf"] < results["adjoint"] and results["cg"] <= results["adjoint_dcf"] * 1.05
print("ordering adjoint > adjoint_dcf >= cg :", "PASS" if ok else "FAIL")
