"""Scratch (h1 challenge): breathing motion model for motion-compensated temporal coupling.
Rank-1 displacement: node b = end-expiration image displaced by a_b * d(r); d from 3D demons between the
end-expiration (nodes 15, 0, 1) and end-inspiration (nodes 7, 8) magnitude images, a_b by a 1D search per node.
estimate(x (16, G) complex) -> d (3, G) voxels, a (16,); warp_mats(d, a, M) -> list of sparse W_b with
(W_b x_{b+1})(r) = x_{b+1}(r + (a_b - a_{b+1}) d(r)) on the support voxels (trilinear)."""
import numpy as np
from scipy import ndimage as ndi
import scipy.sparse as sp


def demons(fixed, moving, iters=80, sigma=2.0):
    import SimpleITK as sitk
    f = sitk.GetImageFromArray(fixed.astype(np.float32)); m = sitk.GetImageFromArray(moving.astype(np.float32))
    m = sitk.HistogramMatching(m, f, 64, 8, True)
    flt = sitk.FastSymmetricForcesDemonsRegistrationFilter()
    flt.SetNumberOfIterations(iters); flt.SetSmoothDisplacementField(True); flt.SetStandardDeviations(sigma)
    field = None
    for shrink in (4, 2, 1):                                   # coarse to fine
        fs = sitk.Shrink(sitk.SmoothingRecursiveGaussian(f, 0.6 * shrink), [shrink] * 3) if shrink > 1 else f
        ms = sitk.Shrink(sitk.SmoothingRecursiveGaussian(m, 0.6 * shrink), [shrink] * 3) if shrink > 1 else m
        if field is not None:
            field = sitk.Resample(field, fs, sitk.Transform(), sitk.sitkLinear, 0.0, sitk.sitkVectorFloat64)
            field = flt.Execute(fs, ms, field)
        else:
            field = flt.Execute(fs, ms)
    d = sitk.GetArrayFromImage(field)                          # (z, y, x, 3) with components (x, y, z) in physical units = voxels
    return np.stack([d[..., 2], d[..., 1], d[..., 0]])          # -> (3, n0, n1, n2) displacement along array axes 0, 1, 2


def warp(img, disp):
    idx = np.meshgrid(*[np.arange(n, dtype=np.float32) for n in img.shape], indexing='ij')
    return ndi.map_coordinates(img, [idx[i] + disp[i] for i in range(3)], order=1, mode='nearest')


def estimate(x, M, verbose=True):
    mag = np.abs(x).astype(np.float32)
    sm = lambda a: ndi.gaussian_filter(a, 0.8)
    E = sm(mag[[15, 0, 1]].mean(0)); I = sm(mag[[7, 8]].mean(0))
    d = demons(I, E) * ndi.gaussian_filter(M.astype(np.float32), 2.0)[None]         # E(r + d) ~ I(r); zero outside the body
    amp = np.sqrt((d ** 2).sum(0)); roi = M & (amp > 0.3 * np.percentile(amp[M], 99))
    a = []
    grid = np.linspace(-0.3, 1.5, 19)
    for b in range(x.shape[0]):
        t = sm(mag[b]); cost = [((warp(E, g_ * d) - t)[roi] ** 2).mean() for g_ in grid]
        i = int(np.argmin(cost)); i = min(max(i, 1), len(grid) - 2)
        c0, c1, c2 = cost[i - 1], cost[i], cost[i + 1]
        a.append(grid[i] + 0.5 * (grid[1] - grid[0]) * (c0 - c2) / (c0 - 2 * c1 + c2 + 1e-30))
    if verbose:
        print('displacement at end inspiration (voxels): 99th pct |d| %.2f, max along axes (LR, AP, SI) %s' % (np.percentile(amp[M], 99), np.round(np.abs(d).reshape(3, -1).max(1), 2)))
        print('node amplitudes a_b:', np.round(a, 2))
    return d, np.array(a)


def _tri(p, lut, shp, nv):
    p0 = np.floor(p).astype(np.int64); f = p - p0
    rows, cols, vals = [], [], []
    for c in range(8):
        o = [(c >> i) & 1 for i in range(3)]
        q = [np.clip(p0[i] + o[i], 0, shp[i] - 1) for i in range(3)]
        w = np.prod([f[i] if o[i] else 1 - f[i] for i in range(3)], axis=0)
        j = lut[np.ravel_multi_index(q, shp)]
        ok = j >= 0
        rows.append(np.arange(nv)[ok]); cols.append(j[ok]); vals.append(w[ok])
    return sp.csr_matrix((np.concatenate(vals).astype(np.float32), (np.concatenate(rows), np.concatenate(cols))), shape=(nv, nv))


def warp_mats(d, a, M, sym=False):
    """Sparse trilinear operators on the support voxels (others never enter): row = target voxel.
    sym=False: W_b with (W_b x_{b+1})(r) = x_{b+1}(r + (a_b - a_{b+1}) d(r)); penalty |x_b - W_b x_{b+1}|^2.
    sym=True:  pairs (Wm_b, Wp_b) evaluating x_b at r - delta/2 and x_{b+1} at r + delta/2, delta = (a_b - a_{b+1}) d:
               both sides get the same interpolation blur, so the penalty no longer taxes fine detail where the shift is fractional."""
    shp = M.shape; idx = np.flatnonzero(M); nv = idx.size
    lut = -np.ones(M.size, np.int64); lut[idx] = np.arange(nv)
    r = np.stack(np.unravel_index(idx, shp)).astype(np.float32)                    # (3, nv)
    dv = d.reshape(3, -1)[:, idx]
    mats = []
    nb = len(a)
    for b in range(nb):
        delta = (a[b] - a[(b + 1) % nb]) * dv
        if sym:
            mats.append((_tri(r - 0.5 * delta, lut, shp, nv), _tri(r + 0.5 * delta, lut, shp, nv)))
        else:
            mats.append(_tri(r + delta, lut, shp, nv))
    return mats, idx
