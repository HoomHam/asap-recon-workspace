#!/usr/bin/env python3
"""
Image-domain SNR for the dynamic Xe-129 recon outputs (AIkill_Dynamic), gas phase and
total dissolved phase, at end-inspiration and end-expiration.

Method (per session, one binning folder, default DIAPHRAGM `d`):
  noise     8 corner cubes (10^3 voxels each) of the 100^3 volume, minus voxels within
            3 voxels of the bin-averaged lung mask. sigma = std per bin (gas is real-valued,
            zero-mean Gaussian background, so no Rician correction).
  lung mask per bin: gas > 5 sigma, 3D connected components >= 500 voxels.
  bins      end-inspiration = bin with largest lung-mask volume, end-expiration = smallest
            (chosen from the images, recorded in the output; bin index differs per binning).
  gas SNR   mean(gas in that bin's mask) / sigma_gas(bin)
  dissolved same mask and bin:
            split   (rbc_tp_separated 1, or untagged pre-40d23a4 image): D = aRBC + aTP
                    (real, zero-mean background) -> mean(D in mask) / std(D in background)
            unsplit (rbc_tp_separated 0): z complex; sigma = mean of std(Re), std(Im) in
                    background; per-voxel Rician bias correction sqrt(max(|z|^2 - 2 sigma^2, 0));
                    SNR = mean(corrected in mask) / sigma
            no dissolved image (gas-only acquisition): blank
Caveat: gridding correlates noise between voxels; values are per-voxel SNR.

Usage:
  python snr_calc.py                      # all sessions -> outputs/snr_2026-09-13/snr_table.csv
  python snr_calc.py --only 2024-03-06_030DN --binning d
"""
import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
import scipy.io as sio
from scipy import ndimage

DYN = Path('/Volumes/HoomHamExt/AIkill_Dynamic')
OUT = Path(__file__).resolve().parents[1] / 'outputs' / 'snr_2026-09-13'
CORNER = 10
THR = 5.0
MIN_BLOB = 500
EXCL_DILATE = 3


def corner_box(shape, k=CORNER):
    box = np.zeros(shape, dtype=bool)
    for a in (slice(0, k), slice(shape[0] - k, shape[0])):
        for b in (slice(0, k), slice(shape[1] - k, shape[1])):
            for c in (slice(0, k), slice(shape[2] - k, shape[2])):
                box[a, b, c] = True
    return box


def lung_mask(img, sigma):
    m = img > THR * sigma
    lab, n = ndimage.label(m)
    if n == 0:
        return m
    sizes = ndimage.sum(m, lab, range(1, n + 1))
    return np.isin(lab, 1 + np.flatnonzero(sizes >= MIN_BLOB))


def session_snr(key, binning):
    folder = DYN / key / binning
    mat = folder / 'recon.mat'
    if not mat.is_file():
        return dict(key=key, error=f'no {binning}/recon.mat')
    names = {n for n, _, _ in sio.whosmat(mat)}
    want = [v for v in ('gas_phase', 'dissolved_phase_real', 'dissolved_phase_imag',
                        'rbc_tp_separated') if v in names]
    m = sio.loadmat(mat, variable_names=want)
    gas = m['gas_phase'].astype(np.float64)            # (bins, z, y, x)
    nb = gas.shape[0]
    corners = corner_box(gas.shape[1:])

    # background: corners minus a dilated bin-averaged lung mask
    gmean = gas.mean(0)
    mmean = lung_mask(gmean, gmean[corners].std())
    bg = corners & ~ndimage.binary_dilation(mmean, iterations=EXCL_DILATE)

    sig_g = np.array([gas[b][bg].std() for b in range(nb)])
    masks = [lung_mask(gas[b], sig_g[b]) for b in range(nb)]
    vols = np.array([mk.sum() for mk in masks])
    snr_g = np.array([gas[b][masks[b]].mean() / sig_g[b] if vols[b] else np.nan
                      for b in range(nb)])
    b_in, b_ex = int(np.argmax(vols)), int(np.argmin(vols))

    row = dict(key=key, binning=binning, bins=nb, bg_voxels=int(bg.sum()),
               insp_bin=b_in, exp_bin=b_ex, lung_vox_insp=int(vols[b_in]), lung_vox_exp=int(vols[b_ex]),
               gas_snr_insp=snr_g[b_in], gas_snr_exp=snr_g[b_ex],
               gas_max_bin=int(np.nanargmax(snr_g)), gas_min_bin=int(np.nanargmin(snr_g)),
               gas_snr_max=float(np.nanmax(snr_g)), gas_snr_min=float(np.nanmin(snr_g)),
               gas_sigma_insp=sig_g[b_in], gas_sigma_exp=sig_g[b_ex],
               gas_snr_per_bin=' '.join(f'{v:.1f}' for v in snr_g),
               lung_vox_per_bin=' '.join(str(int(v)) for v in vols),
               dp_method='', dp_snr_insp=np.nan, dp_snr_exp=np.nan, dp_snr_per_bin='', error='')

    if 'dissolved_phase_real' not in m:
        row['dp_method'] = 'none (gas-only acquisition)'
        return row
    re_, im_ = m['dissolved_phase_real'].astype(np.float64), m['dissolved_phase_imag'].astype(np.float64)
    tag = int(np.asarray(m['rbc_tp_separated']).ravel()[0]) if 'rbc_tp_separated' in m else -1
    split = tag != 0
    row['dp_method'] = 'split: aRBC+aTP' if split else 'unsplit: Rician-corrected |z|'
    if tag == -1:
        row['dp_method'] += ' (untagged, pre-40d23a4)'
    # Total dissolved SNR, phase-independent, same for split and unsplit (2026-09-13 revision).
    # The stored (real, imag) pair is a linear map M of the complex dissolved image z
    # (identity if unsplit; Steve's RBC/TP inversion if split). Complex noise in z is
    # isotropic, so the background covariance C of (real, imag) = sigma^2 M M^T and
    #   q = v^T C^-1 v = |z|^2 / sigma^2   per voxel   (rotation part of M drops out).
    # Pure noise gives E[q] = 2 (verified 2.0-2.1 on 4 sessions). Rician-corrected
    # magnitude SNR = mean_mask sqrt(max(q - 2, 0)).
    # Why not aRBC+aTP: Steve's split often solves at a spurious R crossing (logged R far
    # from target in 60/86 sessions) and makes aRBC/aTP noise ~-0.99 correlated, so the
    # sum is a projection at an arbitrary angle and loses signal.
    row['dp_method'] = ('split' if split else 'unsplit') + ': whitened |z|, Rician-corrected'
    if tag == -1:
        row['dp_method'] += ' (untagged, pre-40d23a4)'
    snr_d = np.full(nb, np.nan)
    snr_rbc = np.full(nb, np.nan)     # split only, DIAGNOSTIC: real part = aRBC (unreliable)
    snr_tp = np.full(nb, np.nan)      # split only, DIAGNOSTIC: imag part = aTP (unreliable)
    for b in range(nb):
        if not vols[b]:
            continue
        C = np.cov(np.stack([re_[b][bg], im_[b][bg]]))
        Ci = np.linalg.inv(C)
        q = Ci[0, 0] * re_[b] ** 2 + 2 * Ci[0, 1] * re_[b] * im_[b] + Ci[1, 1] * im_[b] ** 2
        snr_d[b] = np.sqrt(np.clip(q - 2, 0, None))[masks[b]].mean()
        if split:
            snr_rbc[b] = re_[b][masks[b]].mean() / re_[b][bg].std()
            snr_tp[b] = im_[b][masks[b]].mean() / im_[b][bg].std()
    row['dp_snr_insp'], row['dp_snr_exp'] = snr_d[b_in], snr_d[b_ex]
    row['dp_snr_per_bin'] = ' '.join(f'{v:.1f}' for v in snr_d)
    if np.any(~np.isnan(snr_d)):
        row['dp_max_bin'], row['dp_min_bin'] = int(np.nanargmax(snr_d)), int(np.nanargmin(snr_d))
        row['dp_snr_max'], row['dp_snr_min'] = float(np.nanmax(snr_d)), float(np.nanmin(snr_d))
    for name, arr in (('rbc', snr_rbc), ('tp', snr_tp)):
        if np.any(~np.isnan(arr)):
            row[f'{name}_snr_max'], row[f'{name}_max_bin'] = float(np.nanmax(arr)), int(np.nanargmax(arr))
            row[f'{name}_snr_min'], row[f'{name}_min_bin'] = float(np.nanmin(arr)), int(np.nanargmin(arr))
            row[f'{name}_snr_ei'] = float(arr[b_in])
            row[f'{name}_snr_per_bin'] = ' '.join(f'{v:.1f}' for v in arr)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--binning', default='d', choices=['s', 'p', 'd'])
    ap.add_argument('--only', help='comma-separated session keys (date_id)')
    ap.add_argument('--out', default=str(OUT / 'snr_table.csv'))
    args = ap.parse_args()

    keys = sorted(p.name for p in DYN.iterdir() if p.is_dir() and not p.name.startswith(('.', '_')))
    if args.only:
        keys = [k for k in keys if k in set(args.only.split(','))]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    fields = ['key', 'binning', 'insp_bin', 'exp_bin', 'gas_snr_insp', 'gas_snr_exp',
              'dp_snr_insp', 'dp_snr_exp',
              'gas_snr_max', 'gas_max_bin', 'gas_snr_min', 'gas_min_bin',
              'dp_snr_max', 'dp_max_bin', 'dp_snr_min', 'dp_min_bin', 'dp_method',
              'rbc_snr_max', 'rbc_max_bin', 'rbc_snr_min', 'rbc_min_bin', 'rbc_snr_ei',
              'tp_snr_max', 'tp_max_bin', 'tp_snr_min', 'tp_min_bin', 'tp_snr_ei',
              'rbc_snr_per_bin', 'tp_snr_per_bin', 'lung_vox_insp', 'lung_vox_exp',
              'gas_sigma_insp', 'gas_sigma_exp', 'bg_voxels', 'bins',
              'gas_snr_per_bin', 'dp_snr_per_bin', 'lung_vox_per_bin', 'error']
    with open(args.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        for i, k in enumerate(keys, 1):
            t0 = time.time()
            try:
                r = session_snr(k, args.binning)
            except Exception as e:                       # keep going; record the failure
                r = dict(key=k, binning=args.binning, error=f'{type(e).__name__}: {e}')
            w.writerow({kk: (f'{v:.2f}' if isinstance(v, float) else v) for kk, v in r.items()})
            f.flush()
            if r.get('error'):
                print(f'[{i}/{len(keys)}] {k}: ERROR {r["error"]}', flush=True)
            else:
                print(f'[{i}/{len(keys)}] {k}: gas insp/exp {r["gas_snr_insp"]:.1f}/{r["gas_snr_exp"]:.1f} '
                      f'(bins {r["insp_bin"]}/{r["exp_bin"]})  dp {r["dp_snr_insp"]:.1f}/{r["dp_snr_exp"]:.1f} '
                      f'[{r["dp_method"]}]  {time.time() - t0:.0f}s', flush=True)
    print(f'CSV -> {args.out}')


if __name__ == '__main__':
    sys.exit(main())
