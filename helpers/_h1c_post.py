"""Scratch (h1 challenge): linear post-filters for the SNR arm.
apod(x, tau_ms)   radial k-space weighting exp(-t(|k|)/tau): the readout filter exp(-t/tau) expressed in the image
                  (t(|k|) from the trajectory's mean |k|(sample)); works on (…, n1, n2, n3) complex arrays, 3.5 mm voxels.
"""
import numpy as np
import _h1c as H

_KT = None


def _t_of_k():
    global _KT
    if _KT is None:
        k = H.load()['k']
        kr = np.linalg.norm(k - k[:, :1], axis=-1).mean(0)                 # |k| (start offset removed) per sample
        s_end = int(np.argmax(kr > 0.995 * kr.max()))
        _KT = (np.maximum.accumulate(kr[:s_end + 1]), np.arange(s_end + 1) * H.DWELL)
    return _KT


def apod(x, tau_ms):
    kk, tt = _t_of_k()
    n = x.shape[-3:]
    ax = [np.fft.fftfreq(m, H.DX) for m in n]
    kr = np.sqrt(ax[0][:, None, None] ** 2 + ax[1][None, :, None] ** 2 + ax[2][None, None, :] ** 2)
    w = np.exp(-np.interp(kr, kk, tt, right=tt[-1]) / (tau_ms * 1e-3)).astype(np.float32)
    w /= w.flat[0]
    return np.fft.ifftn(np.fft.fftn(x, axes=(-3, -2, -1)) * w, axes=(-3, -2, -1)).astype(x.dtype)
