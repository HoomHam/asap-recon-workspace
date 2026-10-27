"""Scratch (2026-10-26, PCA Step 1): recon-free per-interleave / per-bin / per-breath tables.

Inputs: old-tree AIkill_Dynamic/<study>/d/signal_pneumo.npz (fid_signal = sum|first 8 samples| per LINE of the
dynamic acquisition: numspec spectral lines first, then gas, diss, gas, diss ...), and the navigator arrays from the
b44d (else old-tree) d/recon.mat (nav_diaphragm_z per 26-ilv frame in voxels, sub-voxel parabola fit; nav_time;
nav_volume per gas interleave = min-max normalised, Savgol-smoothed, linearly interpolated z).

Outputs (workspace/outputs/bin_timing_2026-10-26/):
  <study>_ilv.csv     one row per GAS interleave: t, k0 gas / diss (as stored, scan-normalised 0-1), SIGNAL bin
                      fraction (Steve raw.py), DIAPHRAGM bin fraction (raw.bin on nav_volume), nav_volume,
                      nav z (vox, interpolated from the frames), breath index + seconds since the signal EE.
  <study>_bins.csv    per bin x {SIGNAL, DIAPHRAGM}: hard count, kernel w_eff, median t after EE, kernel-weighted
                      means of nav z (vox), nav_volume, k0 gas, k0 diss (weights = exp(-d^2/2) cyclic, Steve's).
  <study>_breaths.csv per breath (signal EE to next EE): t_EE, period, k0 trough (= EE), k0 peak, t_peak - t_EE,
                      k0 at the next EE, peak/trough, next-EE/peak.
usage: python helpers/_bin_timing_export.py <study> ...
"""
import os
import sys
import numpy as np
import scipy.io as sio

sys.path.insert(0, os.path.dirname(__file__))
from _bin_timing import steve_bin, steve_ee_minima, NBINS, TR, EXT

OUT = os.path.join(os.path.dirname(__file__), '..', 'outputs', 'bin_timing_2026-10-26')


def kernel_w(binfrac):
    """(nilv, NBINS) Steve's cyclic Gaussian weight of each interleave in each bin; 0 for bin -1."""
    w = np.zeros((len(binfrac), NBINS))
    ok = binfrac >= 0
    c = binfrac[ok] * NBINS
    for ib in range(NBINS):
        d = np.minimum.reduce([np.abs(ib - c), np.abs(ib + NBINS - c), np.abs(ib - NBINS - c)])
        w[ok, ib] = np.exp(-d * d / 2.0)
    return w


def load_nav(study):
    for tree in ('AIkill_Dynamic_b44d', 'AIkill_Dynamic'):
        p = f'{EXT}/{tree}/{study}/d/recon.mat'
        if os.path.exists(p):
            m = sio.loadmat(p, variable_names=['nav_volume', 'nav_ilvtime', 'nav_time', 'nav_diaphragm_z'])
            return tree, {k: m[k].ravel().astype(float) for k in m if not k.startswith('__')}
    raise FileNotFoundError(study)


def main(study):
    tree, nav = load_nav(study)
    z = np.load(f'{EXT}/AIkill_Dynamic/{study}/d/signal_pneumo.npz')
    fid = z['fid_signal'].astype(float)
    ngas = len(nav['nav_volume'])
    numspec = 2 * ((len(fid) - 2 * ngas) // 2)   # spectral block is even (500 v2 / 780 v3); an odd leftover line is a trailing unpaired gas line
    if numspec < 0 or numspec > 2000:
        raise ValueError(f'{study}: cannot place the spectral block: {len(fid)} lines, {ngas} gas interleaves')
    gas = fid[numspec::2][:ngas]
    dis = fid[numspec + 1::2][:ngas]
    t = nav['nav_ilvtime'] if len(nav['nav_ilvtime']) == ngas else np.arange(ngas) * TR
    # SIGNAL bins (raw.py:470-510)
    ee = steve_ee_minima(gas, t, 50)
    ee_t = t[ee]
    sig = np.full(ngas, -1.0)
    kb = np.full(ngas, -1)
    for j in range(ngas):
        if len(ee_t) and ee_t[0] <= t[j] <= ee_t[-1]:
            k = min(np.searchsorted(ee_t, t[j], side='right') - 1, len(ee_t) - 2)
            sig[j] = (t[j] - ee_t[k]) / (ee_t[k + 1] - ee_t[k])
            kb[j] = k
    # DIAPHRAGM bins (raw.bin on nav_volume)
    dia = steve_bin(nav['nav_volume'])
    # sub-voxel navigator z per interleave: frames -> interleave times (same linear interp as nav_volume)
    fz = nav['nav_diaphragm_z']
    ft = nav['nav_time']
    good = np.isfinite(fz)
    navz = np.interp(t, ft[good], fz[good], left=np.nan, right=np.nan)
    navz[np.isnan(nav['nav_volume'])] = np.nan
    t_ee = np.where(kb >= 0, t - ee_t[np.clip(kb, 0, None)], np.nan)
    os.makedirs(OUT, exist_ok=True)
    hdr = 'ilv,t_s,k0_gas,k0_diss,sigbin_frac,diabin_frac,nav_volume,nav_z_vox,breath,t_since_EE_s'
    np.savetxt(f'{OUT}/{study}_ilv.csv',
               np.column_stack([np.arange(ngas), t, gas, dis, sig, dia, nav['nav_volume'], navz, kb, t_ee]),
               delimiter=',', header=hdr, comments='', fmt='%d,%.4f,%.5f,%.5f,%.5f,%.5f,%.5f,%.3f,%d,%.3f')
    # per-bin tables
    rows = []
    for name, bf in (('SIGNAL', sig), ('DIAPHRAGM', dia)):
        w = kernel_w(bf)
        ib = np.where(bf >= 0, np.floor(bf * NBINS + 0.5).astype(int) % NBINS, -1)
        for b in range(NBINS):
            wb = w[:, b]
            m = ib == b
            ok = np.isfinite(navz)
            def wmean(x, msk):
                ww = wb * msk
                return np.sum(ww * np.nan_to_num(x)) / max(np.sum(ww), 1e-9)
            tm = np.nanmedian(t_ee[m]) if m.any() else np.nan
            rows.append((name, b, int(m.sum()), wb.sum(), tm, wmean(navz, ok), wmean(nav['nav_volume'], ok),
                         wmean(gas, np.ones(ngas)), wmean(dis, np.ones(ngas))))
    with open(f'{OUT}/{study}_bins.csv', 'w') as f:
        f.write('binning,bin,n_hard,w_eff,t_after_EE_med_s,nav_z_vox_wmean,nav_volume_wmean,k0_gas_wmean,k0_diss_wmean\n')
        for r in rows:
            f.write('%s,%d,%d,%.1f,%.3f,%.3f,%.4f,%.5f,%.5f\n' % r)
    # per-breath k0 table
    with open(f'{OUT}/{study}_breaths.csv', 'w') as f:
        f.write('breath,t_EE_s,period_s,k0_trough,k0_peak,t_peak_after_EE_s,k0_next_EE,peak_over_trough,nextEE_over_peak,k0_diss_trough,k0_diss_peak\n')
        for k in range(len(ee) - 1):
            a, b = ee[k], ee[k + 1]
            seg = gas[a:b + 1]
            ip = a + int(np.argmax(seg))
            dseg = dis[a:b + 1]
            f.write('%d,%.3f,%.3f,%.5f,%.5f,%.3f,%.5f,%.3f,%.3f,%.5f,%.5f\n' % (
                k, t[a], t[b] - t[a], gas[a], gas[ip], t[ip] - t[a], gas[b],
                gas[ip] / max(gas[a], 1e-9), gas[b] / max(gas[ip], 1e-9), dis[a], dseg.max()))
    per = np.diff(ee_t)
    pk = np.array([gas[ee[k] + int(np.argmax(gas[ee[k]:ee[k + 1] + 1]))] for k in range(len(ee) - 1)])
    tr = np.maximum(gas[ee[:-1]], 1e-6)
    r = pk / tr
    r = r[np.isfinite(r) & (tr > 1e-4)]
    print(f'{study}: nav from {tree}; lines {len(fid)} = {numspec} spec + 2x{ngas}; {len(ee)} signal EE, '
          f'period {np.median(per):.2f} s; per-breath k0 peak/trough median {np.median(r):.2f} '
          f'(IQR {np.percentile(r,25):.2f}-{np.percentile(r,75):.2f}, n {len(r)}); diss/gas at EE {np.median(dis[ee[:-1]][tr > 1e-4]/tr[tr > 1e-4]):.3f}')


if __name__ == '__main__':
    for s in sys.argv[1:]:
        try:
            main(s)
        except Exception as e:
            print(f'{s}: FAILED {e}')
