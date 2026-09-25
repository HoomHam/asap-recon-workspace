#!/usr/bin/env python
"""How much of b's phase to keep?  Scale-resolved static residual of aRBC vs the smoothing width.

For sigma in a list: b_s = phase of b low-passed at sigma (inside the lung, weighted), split the dissolved
image with z * exp(-i (phase_b - phase_b_s)), then measure the STATIC component of aRBC between distant
bins (0 vs 7) in three spatial bands: fine (< 1.5 vox), mid (1.5-6 vox), coarse (> 6 vox), and compare
with |dis| (anatomy floor: what a magnitude map that cannot carry phase shows) and with aTM.
Hooman's worry: a low-passed b is still a fixed field on every bin. If the residual static aRBC in every
band sits at the |dis| floor, the fixed field is physics (coil + B0), not an imprint.

Usage: bphase_sigma.py <run_root> <key> [--sigmas 0,1.5,2.5,4,6,10]
Output: outputs/bphase_fix_2026-09-25/<key>_sigma.csv + fig/<key>_sigma.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import snr_calc as S  # noqa: E402
from bphase_confirm import read_mrd, smooth_phase, resplit, OUT, FIG  # noqa: E402

G6 = lambda v: ndimage.gaussian_filter(v, 6)
G15 = lambda v: ndimage.gaussian_filter(v, 1.5)
BANDS = {'fine<1.5': lambda v: v - G15(v), 'mid1.5-6': lambda v: G15(v) - G6(v), 'coarse>6': G6}


def band_static(V, m):
    return {k: float(np.corrcoef(f(V[0])[m], f(V[7])[m])[0, 1]) for k, f in BANDS.items()}


def run(root, key, sigmas):
    gas, dis, b, meta = read_mrd(Path(root) / key / 'd' / 'output.mrd')
    b = b[0]; phi = np.radians(float(meta['rbc_tp_dphi_deg'])); target = float(meta['rbc_tp_target'])
    nb = gas.shape[0]
    gm = gas.mean(0); corners = S.corner_box(gas.shape[1:])
    bg = corners & ~ndimage.binary_dilation(S.lung_mask(gm, gm[corners].std()), iterations=S.EXCL_DILATE)
    masks = [S.lung_mask(gas[k], gas[k][bg].std()) for k in range(nb)]
    many = np.any(masks, 0); w = ndimage.binary_dilation(many, iterations=3).astype(float)
    z_rot = dis.real * np.exp(1j * phi) + dis.imag
    Z = np.abs(z_rot); m07 = masks[0] & masks[7]
    bu = b / np.maximum(np.abs(b), 1e-12)
    rows = [{'sigma': 'dis_floor', **band_static(Z, m07), 'oow': np.nan, 'tm_' + 'fine<1.5': np.nan}]
    for s in sigmas:
        if s == 0:
            R, T = dis.real, dis.imag; tag = 'container (raw b)'
        elif s == 'poly':
            # global quadratic in place of the measured phase INSIDE the mask too (Steve's outside recipe everywhere)
            idx = np.indices(bu.shape).reshape(3, -1).T.astype(float) / bu.shape[0]
            sel = (w > 0).ravel() & (np.abs(bu).ravel() > 0)
            A = np.c_[np.ones(idx.shape[0]), idx, idx ** 2, idx[:, 0] * idx[:, 1], idx[:, 0] * idx[:, 2], idx[:, 1] * idx[:, 2]]
            ang = np.unwrap(np.angle(bu).ravel()) if False else np.angle(bu).ravel()
            coef = np.linalg.lstsq(A[sel], ang[sel], rcond=None)[0]
            poly = (A @ coef).reshape(bu.shape)
            delta = np.where(w > 0, np.angle(np.exp(1j * (np.angle(bu) - poly))), 0.0)
            R, T = resplit(z_rot * np.exp(-1j * delta)[None], masks, phi, target); tag = 'quadratic everywhere'
        else:
            delta = np.where(w > 0, np.angle(np.exp(1j * (np.angle(bu) - smooth_phase(bu, w, float(s))))), 0.0)
            R, T = resplit(z_rot * np.exp(-1j * delta)[None], masks, phi, target); tag = f'b low-pass σ={s}'
        row = {'sigma': tag, **band_static(R, m07), 'oow': 100 * float(((R[7] < 0) | (T[7] < 0))[masks[7]].mean())}
        row.update({'tm_' + k: v for k, v in band_static(T, m07).items()})
        rows.append(row)
        print(row, flush=True)
    t = pd.DataFrame(rows); t.to_csv(OUT / f'{key}_sigma.csv', index=False)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    x = np.arange(len(t) - 1)
    for j, k in enumerate(BANDS):
        ax[j].bar(x - 0.2, t[k].values[1:], 0.4, label='aRBC', color='tab:red')
        ax[j].bar(x + 0.2, t['tm_' + k].values[1:], 0.4, label='aTM', color='tab:blue')
        ax[j].axhline(t[k].values[0], color='k', ls='--', label='|dis| floor (anatomy)')
        ax[j].set_xticks(x); ax[j].set_xticklabels(t.sigma.values[1:], rotation=30, fontsize=7, ha='right')
        ax[j].set_title(f'static corr(bin0,bin7), band {k}', fontsize=9); ax[j].set_ylim(-0.2, 1)
    ax[0].legend(fontsize=7)
    fig.suptitle(f'{key}: static component of the split maps vs how much of b\'s phase is kept', fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / f'{key}_sigma.png', dpi=130); plt.close(fig)
    print(t.round(2).to_string(index=False))


if __name__ == '__main__':
    a = sys.argv[1:]
    sig = a[a.index('--sigmas') + 1].split(',') if '--sigmas' in a else ['0', '1.5', '2.5', '4', '6', '10', 'poly']
    sig = [s if s == 'poly' else float(s) for s in sig]
    run(a[0], a[1], sig)
