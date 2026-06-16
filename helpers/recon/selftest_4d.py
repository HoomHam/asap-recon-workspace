"""Self-test for the 4D operators -- no scanner data. Run: ../.venv/bin/python selftest_4d.py

Checks the three operators that the 4D solver relies on, plus a tiny end-to-end
joint solve:
  1. Finufft4D adjointness  <A x, y> == <x, A^H y>
  2. temporal_diff_op adjointness + circular formula + zero on a constant cine
  3. spatial_wavelet_op adjointness + db4 orthogonality (||Wx|| == ||x||)
  4. recon_4d_joint runs on a tiny synthetic and reduces the data residual
"""

import numpy as np
import asap_recon as ar
import cs_recon_4d as c4

rng = np.random.default_rng(0)
IS = 16
MS = int(IS * 2.4)
B = 3
M = 2000


def adj_err(op, ishape, oshape):
    # <op x, y> == <x, op^H y>  ->  vdot(op(x), y) == vdot(x, op.H(y))
    x = rng.standard_normal(ishape) + 1j * rng.standard_normal(ishape)
    y = rng.standard_normal(oshape) + 1j * rng.standard_normal(oshape)
    lhs = np.vdot(op(x), y)
    rhs = np.vdot(x, op.H(y))
    return abs(lhs - rhs) / abs(lhs)


def ncc(a, b):
    a = np.abs(a).ravel().astype(float); b = np.abs(b).ravel().astype(float)
    a -= a.mean(); b -= b.mean()
    return float((a @ b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-30))


def cone_spiral_traj(n_ilv=120, npts=200):
    """Small golden-angle cone spiral (well-conditioned), grid-index units."""
    kmax = IS / 2 * 0.9
    ga = np.pi * (3 - np.sqrt(5))
    pts = []
    for i in range(n_ilv):
        th = np.arccos(1 - 2 * (i + 0.5) / n_ilv)
        ph = i * ga
        ax = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
        t = np.linspace(0, 1, npts)
        r = kmax * t
        spin = 2 * np.pi * 8 * t
        u = np.cross(ax, [0, 0, 1.0]); u /= np.linalg.norm(u) + 1e-12
        v = np.cross(ax, u)
        k = (r * np.cos(spin))[:, None] * u + (r * np.sin(spin))[:, None] * v \
            + (r * 0.3)[:, None] * ax
        pts.append(k)
    k = np.concatenate(pts)
    return k * MS / IS + MS / 2


# shared spiral-ish trajectory (well-conditioned for the actual recon in test 4)
traj_grid = cone_spiral_traj()
M = traj_grid.shape[0]
traj_rad = ar.grid_to_radians(traj_grid, MS)

passed = True

# 1. Finufft4D adjointness
A4 = c4.Finufft4D(traj_rad, B, IS)
e = adj_err(A4, (B, IS, IS, IS), (B, M))
ok = e < 1e-6
passed &= ok
print(f"1. Finufft4D adjointness    rel err {e:.2e}  {'PASS' if ok else 'FAIL'}")

# 2. temporal diff op
D = c4.temporal_diff_op((B, IS, IS, IS))
e = adj_err(D, (B, IS, IS, IS), (B, IS, IS, IS))
ok = e < 1e-10
x = rng.standard_normal((B, IS, IS, IS)) + 1j * rng.standard_normal((B, IS, IS, IS))
manual = np.roll(x, -1, axis=0) - x
formula_ok = np.allclose(D(x), manual)
xc = np.tile(x[:1], (B, 1, 1, 1))           # constant across bins
zero_ok = np.abs(D(xc)).max() < 1e-10
passed &= ok and formula_ok and zero_ok
print(f"2. temporal D adjointness   rel err {e:.2e}  {'PASS' if ok else 'FAIL'}"
      f"  | circular formula {'PASS' if formula_ok else 'FAIL'}"
      f"  | const->0 {'PASS' if zero_ok else 'FAIL'}")

# 3. spatial per-bin wavelet
W = c4.spatial_wavelet_op((B, IS, IS, IS))
e = adj_err(W, (B, IS, IS, IS), W.oshape)
ortho = abs(np.linalg.norm(W(x)) - np.linalg.norm(x)) / np.linalg.norm(x)
ok = e < 1e-10 and ortho < 1e-6
passed &= ok
print(f"3. wavelet adjointness      rel err {e:.2e}  | orthogonality {ortho:.2e}"
      f"  {'PASS' if ok else 'FAIL'}")

# 4. tiny end-to-end joint solve (shared data, all bins see all samples)
#    truth: one Gaussian blob; synth samples via the 3D forward; M_samp = ones so
#    every bin reconstructs the same image (temporal TV ~ 0) -- just checks the
#    solver runs and the data residual drops.
xs = np.arange(IS) - IS / 2
X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
truth = np.exp(-(X**2 + Y**2 + Z**2) / (2 * 4.0**2)).astype(np.complex128)
y = ar.forward(traj_rad, truth)
M_samp = np.ones((B, M))
xrec, base = c4.recon_4d_joint(traj_grid, y, M_samp, N=IS, MS=MS, lam_s_rel=0.01,
                               lam_t_rel=0.05, max_iter=40, max_power_iter=10,
                               verbose=False)
# all bins share data + uniform membership -> every bin should recover the blob,
# and bins should be near-identical (temporal TV ~ 0). Scale-invariant checks.
corr = ncc(xrec[0], truth)
bin_spread = np.linalg.norm(xrec[0] - xrec[1]) / (np.linalg.norm(xrec[0]) + 1e-30)
ok = np.all(np.isfinite(xrec)) and corr > 0.7 and bin_spread < 0.1
passed &= ok
print(f"4. joint solve recovers blob  ncc(xrec,truth) {corr:.3f}>0.7  "
      f"bin spread {bin_spread:.3f}<0.1  shape {xrec.shape}  {'PASS' if ok else 'FAIL'}")

print("\nALL PASS" if passed else "\nFAILURES ABOVE")
