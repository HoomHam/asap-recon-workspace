"""Scratch (h1 challenge): model-based self-navigation. Each line's breathing coordinate is corrected from the data:
Gauss-Newton step along the cyclic coordinate, delta = Re<r, dU/dc> / |dU/dc|^2 with r = line - model prediction -
fixed per-interleave pattern, then inverse-variance smoothing along time (breathing is slow against 17 ms lines).
Lines listed in `exclude` (held-out) do not contribute their own data: they take their neighbours' estimate."""
import numpy as np
from scipy.ndimage import gaussian_filter1d
import _h1c as H


def refine(d, x, train, exclude=None, band=(0, 198), smooth=10.0, step=1.0, U=None):
    U = d.fwd_all(x)[0] if U is None else U                                # (nb, 832, 510)
    lines = d.ss
    il = lines % H.NILV
    b0, f = d.hat(lines); b1 = (b0 + 1) % d.nb
    sl = slice(*band)
    y = d.v[d.chs[0], lines, H.KILL:][:, sl]
    u0 = U[b0, il][:, sl]; u1 = U[b1, il][:, sl]
    r = y - ((1 - f)[:, None] * u0 + f[:, None] * u1)
    dU = (u1 - u0) * d.nb
    istr = np.isin(lines, train)
    rb = np.zeros((H.NILV, r.shape[1]), np.complex128); cnt = np.bincount(il[istr], minlength=H.NILV)
    np.add.at(rb, il[istr], r[istr]); rb[cnt > 0] /= cnt[cnt > 0][:, None]
    r = r - rb[il]
    num = (np.conj(dU) * r).real.sum(1); den = (np.abs(dU) ** 2).sum(1)
    use = np.ones(lines.size, bool)
    if exclude is not None:
        use &= ~np.isin(lines, exclude)
    # time axis = line index (uniform, 17.2 ms); fill a full-length vector so that gaps are handled by the weights
    N = d.L; nn = np.zeros(N); dd = np.zeros(N)
    nn[lines[use]] = num[use]; dd[lines[use]] = den[use]
    delta = gaussian_filter1d(nn, smooth, mode='nearest') / (gaussian_filter1d(dd, smooth, mode='nearest') + 1e-30)
    c_new = d.c.copy()
    c_new[lines] = (d.c[lines] + step * delta[lines]) % 1.0
    info = dict(delta_rms=float(np.sqrt((delta[lines] ** 2).mean())), raw_snr=float(np.sqrt((num ** 2 / (den / 2 + 1e-30)).mean())))
    return c_new, delta, info
