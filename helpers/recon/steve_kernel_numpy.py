"""Steve's cudarecon/cudarenorm, faithfully reimplemented in numpy. CPU-only.

Replicates recon.py line-for-line in vectorized form (single bin, one channel):
  - readout filter      fwt = exp(-((idx % npts)/smoothing)^2)   [recon.py:99]
  - skip |raw| < eps                                              [recon.py:101]
  - trajectory wrap     kidx = idx % nuniquesmp                   [recon.py:103]
  - box cx-2 .. cx+1 per axis, edge-clamped (asymmetric, as his)  [recon.py:118]
  - Gaussian weight     exp(-dsq/0.2)                             [recon.py:124-129]
  - accumulate k += wt*raw*fwt, knorm += wt  (fwt NOT in knorm)   [recon.py:134-136]
  - renorm k/=knorm, zero cells with knorm < 1e-5 ("KLUGE")       [recon.py:59-64]
  - flat C-index ix*MS^2+iy*MS+iz, reshaped order='F'             [results.py:92]
    => returned array axes are (z, y, x), exactly like Steve's kspace3d
  - fftshift(fftn(ifftshift(.))), crop ISLL:ISUL                  [results.py:99,103]

Differences vs GPU run: float64 numpy vs float32 atomics; summation order.
Expect agreement to ~1e-6 relative, not bit-exact.
"""

import numpy as np

EPS = 1.0e-16
KDIST0SQ = 0.2
BXSZ = 2
KNORM_FLOOR = 1e-5


def steve_recon(traj_grid, data, npts, MS=240, IS=100, smoothing=0.0,
                axes="steve"):
    """traj_grid: (Nuniq,3) grid units. data: (M,) complex. npts: samples per
    interleave (for the readout filter). smoothing: Steve's gplb (0 = off).
    axes='steve' returns his (z,y,x) layout (comparable to savedbin*.npy);
    axes='xyz' transposes to match asap_recon/finufft images."""
    data = np.asarray(data, dtype=np.complex128)
    M = len(data)
    idx = np.arange(M)

    rawval = data.copy()
    if smoothing > EPS:
        rawval = rawval * np.exp(-(((idx % npts) / smoothing) ** 2))

    keep = np.abs(rawval) >= EPS
    rawval = rawval[keep]
    kidx = idx[keep] % traj_grid.shape[0]
    kx, ky, kz = (traj_grid[kidx, i] for i in range(3))
    cx, cy, cz = ((k + 0.5).astype(np.int64) for k in (kx, ky, kz))

    MS2 = MS * MS
    k_acc = np.zeros(MS ** 3, dtype=np.complex128)
    knorm = np.zeros(MS ** 3, dtype=np.float64)

    # vectorize his triple loop: 4^3 offsets, edge-clamp exactly like
    # range(max(0,c-2), min(MS,c+2)) — offset valid iff resulting index is
    # inside [max(0,c-2), min(MS,c+2)) for that sample
    offs = range(-BXSZ, BXSZ)  # -2..+1, his asymmetric box
    for dx in offs:
        ix = cx + dx
        vx = (ix >= np.maximum(0, cx - BXSZ)) & (ix < np.minimum(MS, cx + BXSZ)) \
             & (ix >= 0) & (ix < MS)
        for dy in offs:
            iy = cy + dy
            vy = vx & (iy >= 0) & (iy < MS)
            for dz in offs:
                iz = cz + dz
                v = vy & (iz >= 0) & (iz < MS)
                if not v.any():
                    continue
                dsq = (kx[v] - ix[v]) ** 2 + (ky[v] - iy[v]) ** 2 + (kz[v] - iz[v]) ** 2
                wt = np.exp(-dsq / KDIST0SQ)
                flat = ix[v] * MS2 + iy[v] * MS + iz[v]
                k_acc += np.bincount(flat, weights=wt * rawval[v].real, minlength=MS ** 3) \
                    + 1j * np.bincount(flat, weights=wt * rawval[v].imag, minlength=MS ** 3)
                knorm += np.bincount(flat, weights=wt, minlength=MS ** 3)

    # cudarenorm, oldway path incl. the kluge
    nz = knorm > 0
    k_acc[nz] = k_acc[nz] / knorm[nz]
    k_acc[knorm < KNORM_FLOOR] = 0

    kspace3d = np.reshape(k_acc, (MS, MS, MS), order="F")   # (z,y,x) like Steve
    rspace = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(kspace3d)))

    MSc = int(MS / 2 + 0.1)
    ll = int(MSc - IS / 2 + 0.1)
    img = rspace[ll:ll + IS, ll:ll + IS, ll:ll + IS]
    return img if axes == "steve" else np.transpose(img, (2, 1, 0))
