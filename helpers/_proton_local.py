"""Scratch (2026-10-01): local FINUFFT recon toolkit for the 006KL 1H dynamic (delay / apodization tests).

load()          raw (512, ch, lines), measured water traj (832, 512, 3), TR
endexp_lines()  line indices in end-expiratory k0 bins (Steve's bin(): bins 0-2 + 12-15 of 16)
delay_traj()    k(n - tau) per axis (tau in samples of 5 us), clamped at the start
pipe_dcf()      Pipe-Menon DCF on the 832 unique interleaves (cycles/FOV units)
recon()         adjoint NUFFT, DCF-weighted, optional readout window, 2-coil SOS
"""
import os, sys
import numpy as np
import mapvbvd, finufft
import sigpy.mri as smr

HERE = os.path.dirname(os.path.abspath(__file__))
DAT = '/Volumes/HoomHamExt/_5t_images_roundtrip/Images/2024-03-12/006KL/Dynamic Proton/meas_MID00174_FID08670_fa_spiral_dyn_fancy_v3_20240130.dat'
TRAJ = os.path.join(HERE, '..', 'data', 'h1', 'traj', 'fa_spiral_dyn_fancy_v3_20240130_water_fov350_gp.npy')
KILL = 2                                   # ADC transient samples


def load():
    tw = mapvbvd.mapVBVD(DAT, quiet=True)
    tw = tw[-1] if isinstance(tw, list) else tw
    tw.image.flagRemoveOS = False
    tw.image.squeeze = True
    tr = float(tw.hdr.MeasYaps[('alTR', '0')]) * 1e-6
    raw = tw.image.unsorted()
    raw = raw if raw.ndim == 3 else raw[:, None, :]
    k = np.load(TRAJ).reshape(832, 512, 3)
    return raw, k, tr


def endexp_lines(nl, nbins=16, keep=(0, 1, 2, 12, 13, 14, 15)):
    sys.path.insert(0, HERE)
    sys.path.insert(0, os.path.join(HERE, '..', '..'))
    from _proton_resp_pneumo import k0_resp
    import raw as rawmod
    _, r, _ = k0_resp(DAT)
    v = (-r).astype(float).copy()
    R = rawmod.raw()
    R.rescale(v)
    b = np.floor(R.bin(v) * nbins).astype(int).clip(0, nbins - 1)
    ss = np.arange(nl) >= 832                # steady state only
    return np.where(np.isin(b, keep) & ss)[0], b


def delay_traj(k, tau):
    """tau: scalar or (3,) in samples; positive = data arrive later than the trajectory says."""
    tau = np.broadcast_to(np.asarray(tau, float), (3,))
    n = np.arange(k.shape[1], dtype=float)
    out = np.empty_like(k)
    for a in range(3):
        for i in range(k.shape[0]):
            out[i, :, a] = np.interp(n - tau[a], n, k[i, :, a])
    return out


def pipe_dcf(k, N, DX, it=30):
    w = smr.pipe_menon_dcf(k.reshape(-1, 3) * N * DX, img_shape=(N, N, N), max_iter=it, show_pbar=False)
    return np.abs(w).reshape(k.shape[:2])


def recon(raw, k, lines, dcf, N, DX, win=None, kmax=None):
    il = lines % k.shape[0]
    kk = k[il][:, KILL:]                                    # (L, S, 3)
    w = dcf[il][:, KILL:].copy()
    if win is not None:
        w = w * win[KILL:][None, :]
    keep = np.ones(kk.shape[:2], bool)
    if kmax is not None:
        keep = np.linalg.norm(kk, axis=-1) <= kmax
    kx = (2 * np.pi * DX * kk[keep]).astype(np.float64)
    sos = 0
    for c in range(raw.shape[1]):
        d = raw[KILL:, c, :][:, lines].T                    # (L, S)
        im = finufft.nufft3d1(kx[:, 0], kx[:, 1], kx[:, 2], (d[keep] * w[keep]).astype(np.complex128),
                              (N, N, N), eps=1e-4, nthreads=0)
        sos = sos + np.abs(im) ** 2
    return np.sqrt(sos)


def sharpness(x, mask=None):
    """Fuderer entropy (lower = sharper) and normalized gradient energy inside mask (higher = sharper)."""
    m = mask if mask is not None else x > 0.15 * np.percentile(x, 99)
    b = x[m] / np.sqrt((x[m] ** 2).sum())
    ent = -(b * np.log(b + 1e-30)).sum()
    g = sum(np.gradient(x, axis=a) ** 2 for a in range(3))
    return ent, g[m].sum() / (x[m] ** 2).sum()
