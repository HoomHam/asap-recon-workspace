#!/usr/bin/env python
"""Remove the bin-invariant FINE-SCALE phase from the dissolved image before the RBC/TM split.

Finding (ledger s11 06:00, F79): the voxel-scale texture of aRBC is identical in every bin inside the lung
(corr(bin0, bin7) of the < 1.5-vox component: aRBC 0.84, aTM 0.77, |dis| 0.47, gas 0.33, background 0).
The only full-resolution, frame-invariant factor in the pipeline is calcb's b = conj(F_avg)/|F_avg|,
whose phase is measured voxel by voxel inside the 10-sigma support (a sharp lung-shaped boundary) and
smooth (polynomial) outside. The dissolved image is smooth, so its own fine-scale phase is noise; b's
fine-scale phase is imprinted on every bin and lands in the quadrature (RBC) channel as a static speckle.

Proxy fix from existing maps (no b available offline): in the container's common frame z_rot (TM axis
= 0, all bins consistent), the bin-mean of z_rot carries the static phase; its fine-scale part
delta(x) = angle(mean_b z_rot) - smooth(angle(mean_b z_rot), sigma) is removed from every bin,
z' = z_rot * exp(-i delta), then the split is re-swept to the spectroscopic target and re-applied.
Real fine-scale dissolved phase that is static across bins is removed too - at this SNR it is not
measurable anyway (per-voxel RBC SNR ~2).

Usage: bphase_fix.py <key> [--sigma 2.5] [--write]   (--write stores d/recon_resplit_xecs_bfix.mat on Ext)
Outputs: outputs/bphase_fix_2026-09-25/<key>.json, fig/<key>.png
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import scipy.io as sio
from scipy import ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import snr_calc as S  # noqa: E402
from resplit_rbctp import sweep_phase, split_maps  # noqa: E402

SRC = Path('/Volumes/HoomHamExt/AIkill_Dynamic')
WS = HERE.parents[1]
OUT = WS / 'outputs' / 'bphase_fix_2026-09-25'
FIG = OUT / 'fig'
fine = lambda v, s=1.5: v - ndimage.gaussian_filter(v, s)


def smooth_phase(c, w, sigma):
    num = ndimage.gaussian_filter(c.real * w, sigma) + 1j * ndimage.gaussian_filter(c.imag * w, sigma)
    return np.angle(num / (ndimage.gaussian_filter(w, sigma) + 1e-9))


def run(key, sigma=2.5, write=False):
    p = SRC / key / 'd'
    G = sio.loadmat(p / 'recon.mat', variable_names=['gas_phase'])['gas_phase'].astype(float)
    r = sio.loadmat(p / 'recon_resplit_xecs.mat')
    phi = np.radians(float(r['dphi_new_deg'].ravel()[0])); target = float(r['target'].ravel()[0])
    R0, T0 = r['aRBC'].astype(float), r['aTP'].astype(float)
    nb = G.shape[0]
    gm = G.mean(0); corners = S.corner_box(G.shape[1:])
    bg = corners & ~ndimage.binary_dilation(S.lung_mask(gm, gm[corners].std()), iterations=S.EXCL_DILATE)
    masks = [S.lung_mask(G[b], G[b][bg].std()) for b in range(nb)]
    many = np.any(masks, 0)
    z_rot = R0 * np.exp(1j * phi) + T0                       # common frame, consistent across bins
    zm = z_rot.mean(0)
    w = ndimage.binary_dilation(many, iterations=3).astype(float)
    ang = np.angle(zm)
    delta = np.where(w > 0, np.angle(np.exp(1j * (ang - smooth_phase(zm, w, sigma)))), 0.0)
    z1 = z_rot * np.exp(-1j * delta)[None]
    R1 = np.zeros_like(R0); T1 = np.zeros_like(T0); ph = np.full(nb, np.nan)
    for b in range(nb):
        mk = masks[b]
        if not mk.any():
            continue
        pb, Rr, sR, sT, nv = sweep_phase(z1[b][mk].sum(), phi, target)
        ph[b] = pb
        if np.isnan(pb):
            continue
        R1[b], T1[b] = split_maps(z1[b] * np.exp(1j * pb), phi)
    Z = np.abs(z_rot)
    res = {'key': key, 'sigma': sigma, 'phi_k0_deg': float(np.degrees(phi)),
           'delta_lung_deg': {'median_abs': float(np.degrees(np.median(np.abs(delta[many])))),
                              'p90_abs': float(np.degrees(np.percentile(np.abs(delta[many]), 90)))},
           'resweep_ph_rad': [float(v) for v in ph]}
    m07 = masks[0] & masks[7]; m312 = masks[3] & masks[12]
    for name, R, T in (('before', R0, T0), ('after', R1, T1)):
        d = {}
        for n, V in (('rbc', R), ('tm', T), ('dis', Z)):
            d[f'static_fine_corr07_{n}'] = float(np.corrcoef(fine(V[0])[m07], fine(V[7])[m07])[0, 1])
            d[f'static_fine_corr312_{n}'] = float(np.corrcoef(fine(V[3])[m312], fine(V[12])[m312])[0, 1])
        d['oow_pct_bin7'] = 100 * float(((R[7] < 0) | (T[7] < 0))[masks[7]].mean())
        d['lung_rbc_over_bgsd_bin7'] = float(np.median(np.abs(R[7][masks[7]])) / R[7][bg].std())
        d['corr_rbc_vs_before_bin7'] = float(np.corrcoef(R[7][masks[7]], R0[7][masks[7]])[0, 1])
        res[name] = d
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(exist_ok=True)
    json.dump(res, open(OUT / f'{key}.json', 'w'), indent=1)
    # figure
    ys = np.where(masks[7].any(axis=(0, 2)))[0]; y = int(ys[len(ys) // 2])
    sl = lambda v: v[:, y, :][::-1][15:85, 10:90]
    fig, axs = plt.subplots(3, 5, figsize=(15, 9.3), facecolor='black')
    vr = np.percentile(np.maximum(R0[7], 0)[masks[7]], 99.5)
    for j, b in enumerate((0, 7, 12)):
        axs[0, j].imshow(sl(np.maximum(R0[b], 0)), cmap='gray', vmin=0, vmax=vr); axs[0, j].set_title(f'aRBC before, bin {b}', fontsize=8, color='white')
        axs[1, j].imshow(sl(np.maximum(R1[b], 0)), cmap='gray', vmin=0, vmax=vr); axs[1, j].set_title(f'aRBC after, bin {b}', fontsize=8, color='white')
    axs[0, 3].imshow(sl(np.degrees(delta)), cmap='twilight', vmin=-30, vmax=30); axs[0, 3].set_title(f'removed static fine phase δ [deg] (σ={sigma})', fontsize=8, color='white')
    axs[0, 4].imshow(sl(np.maximum(R0.min(0), 0)), cmap='gray', vmin=0, vmax=vr); axs[0, 4].set_title('aRBC MIN over bins, before', fontsize=8, color='white')
    axs[1, 3].imshow(sl(np.maximum(R1.min(0), 0)), cmap='gray', vmin=0, vmax=vr); axs[1, 3].set_title('aRBC MIN over bins, after', fontsize=8, color='white')
    axs[1, 4].imshow(sl(fine(R0[7])), cmap='gray', vmin=-vr / 3, vmax=vr / 3); axs[1, 4].set_title('fine (<1.5 vox) aRBC bin 7, before', fontsize=8, color='white')
    vt = np.percentile(T0[7][masks[7]], 99.5)
    for j, (img, t) in enumerate([(np.maximum(T0[7], 0), 'aTM before, bin 7'), (np.maximum(T1[7], 0), 'aTM after, bin 7'),
                                  (Z[7], '|dis| bin 7'), (G[7], 'gas bin 7'), (fine(R1[7]), 'fine aRBC bin 7, after')]):
        vm = vt if 'aTM' in t else (np.percentile(img[masks[7]], 99.5) if 'fine' not in t else vr / 3)
        axs[2, j].imshow(sl(img), cmap='gray', vmin=(0 if 'fine' not in t else -vr / 3), vmax=vm); axs[2, j].set_title(t, fontsize=8, color='white')
    for ax in axs.ravel():
        ax.set_facecolor('black'); ax.set_xticks([]); ax.set_yticks([])
    bf, af = res['before'], res['after']
    fig.suptitle(f'{key} · static fine-scale corr(bin0,bin7) aRBC {bf["static_fine_corr07_rbc"]:.2f} → {af["static_fine_corr07_rbc"]:.2f}, '
                 f'aTM {bf["static_fine_corr07_tm"]:.2f} → {af["static_fine_corr07_tm"]:.2f}, |dis| {bf["static_fine_corr07_dis"]:.2f} · '
                 f'δ in lung median {res["delta_lung_deg"]["median_abs"]:.1f}° p90 {res["delta_lung_deg"]["p90_abs"]:.1f}° · '
                 f'out-of-wedge {bf["oow_pct_bin7"]:.1f} → {af["oow_pct_bin7"]:.1f} %', fontsize=8, color='white')
    fig.tight_layout(); fig.savefig(FIG / f'{key}.png', dpi=120, facecolor='black'); plt.close(fig)
    if write:
        sio.savemat(p / 'recon_resplit_xecs_bfix.mat', dict(aRBC=R1.astype(np.float32), aTP=T1.astype(np.float32),
                    delta_deg=np.degrees(delta).astype(np.float32), sigma=sigma, ph_new=ph, dphi_new_deg=np.degrees(phi),
                    target=target, made_by='workspace/helpers/atlas/bphase_fix.py 2026-09-25 (proxy: static fine phase from bin-mean z_rot)'))
    return res


if __name__ == '__main__':
    a = sys.argv[1:]
    sigma = float(a[a.index('--sigma') + 1]) if '--sigma' in a else 2.5
    res = run(a[0], sigma, '--write' in a)
    for k in ('before', 'after'):
        d = res[k]
        print(f'{k:6s} static fine corr07 rbc {d["static_fine_corr07_rbc"]:.2f} tm {d["static_fine_corr07_tm"]:.2f} dis {d["static_fine_corr07_dis"]:.2f} | '
              f'corr312 rbc {d["static_fine_corr312_rbc"]:.2f} | oow {d["oow_pct_bin7"]:.1f}% | lungRBC/bgsd {d["lung_rbc_over_bgsd_bin7"]:.2f} | corr vs before {d["corr_rbc_vs_before_bin7"]:.3f}')
    print('delta in lung:', res['delta_lung_deg'])
