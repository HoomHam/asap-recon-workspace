"""Scratch (h1 challenge): Marchenko-Pastur PCA denoising across the 16 phases (complex, patch-wise, sliding).
For each patch (p^3 voxels x 16 phases) keep the principal components above the MP noise edge (Veraart 2016
estimator of sigma and rank), average the overlapping patch estimates. Linear subspace per patch, data-driven,
no sparsity transform. mppca(x (16, X, Y, Z) complex, p=5, stride=2) -> denoised, sigma map, rank map."""
import numpy as np


def _mp_rank(s2, M, N):
    """s2: eigenvalues (descending) of X^H X / N for X (N x M), M <= N. Returns (rank kept, sigma^2)."""
    lam = s2
    for p in range(M - 1):
        g = (M - p) / N
        sig2 = lam[p + 1:].mean() if p + 1 < M else lam[-1]
        rng_mp = 4 * sig2 * np.sqrt(g)
        if lam[p + 1] - lam[-1] < rng_mp:                 # the remaining eigenvalues fit inside one MP bulk
            return p + 1, lam[p + 1:].sum() / ((M - p - 1) * (1 - 0))
    return M, 0.0


def mppca(x, p=5, stride=2, mask=None):
    T = x.shape[0]; shp = x.shape[1:]
    out = np.zeros_like(x); wsum = np.zeros(shp, np.float32); sig = np.zeros(shp, np.float32); rk = np.zeros(shp, np.float32)
    N = p ** 3
    for i in range(0, shp[0] - p + 1, stride):
        for j in range(0, shp[1] - p + 1, stride):
            for k in range(0, shp[2] - p + 1, stride):
                if mask is not None and not mask[i + p // 2, j + p // 2, k + p // 2]:
                    continue
                X = x[:, i:i + p, j:j + p, k:k + p].reshape(T, N).T            # (N, T)
                mu = X.mean(0, keepdims=True); Xc = X - mu
                u, s, vh = np.linalg.svd(Xc, full_matrices=False)
                r, s2 = _mp_rank(s ** 2 / N, T, N)
                r = min(r, T - 1)
                Xd = (u[:, :r] * s[:r]) @ vh[:r] + mu
                out[:, i:i + p, j:j + p, k:k + p] += Xd.T.reshape(T, p, p, p)
                wsum[i:i + p, j:j + p, k:k + p] += 1; sig[i:i + p, j:j + p, k:k + p] += np.sqrt(max(s2, 0)); rk[i:i + p, j:j + p, k:k + p] += r
    ok = wsum > 0
    out[:, ok] /= wsum[ok]; out[:, ~ok] = x[:, ~ok]
    sig[ok] /= wsum[ok]; rk[ok] /= wsum[ok]
    return out, sig, rk
