#!/usr/bin/env python
"""Container with calcb phase smoothing (fork 19b7365, sigma 2.5) vs the raw-b container and the offline true-b fix.

Checks per session:
  1. dissolved split: new container aRBC vs offline true-b fix at the same sigma (C44 machinery) -> corr ~ 1
  2. static-texture metric (fine corr bin0/bin7) raw -> new, out-of-wedge raw -> new
  3. gas image real(F*b): change vs raw container (should be ~0: cos of ~2 deg)
  4. frames for the eye: aRBC bins 0/4/7/11 raw container vs new container, coronal + sagittal
Usage: bsmooth_compare.py <key>  (roots fixed below)
Outputs: outputs/bphase_fix_2026-09-25/<key>_bsmooth.json, fig/<key>_bsmooth_frames.png
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
import snr_calc as S  # noqa: E402
from bphase_confirm import read_mrd, smooth_phase, resplit, metrics, OUT, FIG  # noqa: E402

RAW = Path('/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/tyger_bexport_2026-09-25')     # c8366c3: raw b + export
NEW = Path('/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/tyger_bsmooth_2026-09-25')     # 19b7365: smoothed b
fine = lambda v, s=1.5: v - ndimage.gaussian_filter(v, s)


def run(key):
    g0, d0, b0, m0 = read_mrd(RAW / key / 'd' / 'output.mrd')
    g1, d1, b1, m1 = read_mrd(NEW / key / 'd' / 'output.mrd')
    phi0, phi1 = np.radians(float(m0['rbc_tp_dphi_deg'])), np.radians(float(m1['rbc_tp_dphi_deg']))
    nb = g0.shape[0]
    gm = g0.mean(0); corners = S.corner_box(g0.shape[1:])
    bg = corners & ~ndimage.binary_dilation(S.lung_mask(gm, gm[corners].std()), iterations=S.EXCL_DILATE)
    masks = [S.lung_mask(g0[k], g0[k][bg].std()) for k in range(nb)]
    many = np.any(masks, 0); w = ndimage.binary_dilation(many, iterations=3).astype(float)
    Z0 = np.abs(d0.real * np.exp(1j * phi0) + d0.imag); Z1 = np.abs(d1.real * np.exp(1j * phi1) + d1.imag)
    # offline true-b fix at sigma 2.5 from the raw container (C44 recipe)
    bu = b0[0] / np.maximum(np.abs(b0[0]), 1e-12)
    delta = np.where(w > 0, np.angle(np.exp(1j * (np.angle(bu) - smooth_phase(bu, w, 2.5)))), 0.0)
    z_rot0 = d0.real * np.exp(1j * phi0) + d0.imag
    Rf, Tf = resplit(z_rot0 * np.exp(-1j * delta)[None], masks, phi0, float(m0['rbc_tp_target']))
    res = {'key': key, 'dphi_raw': float(np.degrees(phi0)), 'dphi_new': float(np.degrees(phi1)),
           'target_raw': m0['rbc_tp_target'], 'target_new': m1['rbc_tp_target'],
           'raw': metrics(d0.real, d0.imag, Z0, masks, bg), 'new_container': metrics(d1.real, d1.imag, Z1, masks, bg),
           'offline_trueb_s2.5': metrics(Rf, Tf, Z0, masks, bg)}
    mk = masks[7]
    res['corr_aRBC_new_vs_offline_bin7'] = float(np.corrcoef(d1.real[7][mk], Rf[7][mk])[0, 1])
    res['corr_aRBC_new_vs_raw_bin7'] = float(np.corrcoef(d1.real[7][mk], d0.real[7][mk])[0, 1])
    res['corr_aTM_new_vs_raw_bin7'] = float(np.corrcoef(d1.imag[7][mk], d0.imag[7][mk])[0, 1])
    res['dis_mag_relrms_new_vs_raw'] = float(np.sqrt(np.mean((Z1[7][mk] - Z0[7][mk]) ** 2)) / np.sqrt(np.mean(Z0[7][mk] ** 2)))
    res['gas_relrms_new_vs_raw_bin7'] = float(np.sqrt(np.mean((g1[7][mk] - g0[7][mk]) ** 2)) / np.sqrt(np.mean(g0[7][mk] ** 2)))
    res['gas_lungmean_ratio_new_over_raw'] = float(g1[7][mk].mean() / g0[7][mk].mean())
    b1u = b1[0] / np.maximum(np.abs(b1[0]), 1e-12)
    res['b_phase_diff_deg_in_lung'] = {'rms': float(np.degrees(np.sqrt(np.mean(np.angle(bu * np.conj(b1u))[many] ** 2)))),
                                       'p90_abs': float(np.degrees(np.percentile(np.abs(np.angle(bu * np.conj(b1u))[many]), 90)))}
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(exist_ok=True)
    json.dump(res, open(OUT / f'{key}_bsmooth.json', 'w'), indent=1)
    # frames
    ys = np.where(mk.any(axis=(0, 2)))[0]; y = int(ys[len(ys) // 2]); xs = np.where(mk.any(axis=(0, 1)))[0]; x = int(xs[len(xs) // 4])
    bins = (0, 4, 7, 11)
    fig, axs = plt.subplots(4, 4, figsize=(12.5, 12.5), facecolor='black')
    vr = np.percentile(np.maximum(d0.real[7], 0)[mk], 99.5)
    for j, b in enumerate(bins):
        for i, (V, t) in enumerate([(d0.real[b], 'raw b'), (d1.real[b], 'smoothed b')]):
            axs[i, j].imshow(np.maximum(V, 0)[:, y, :][::-1][15:85, 10:90], cmap='gray', vmin=0, vmax=vr)
            axs[i, j].set_title(f'aRBC bin {b} — {t} (cor)', fontsize=8, color='white')
            axs[2 + i, j].imshow(np.maximum(V, 0)[:, :, x][::-1][15:85, 10:90], cmap='gray', vmin=0, vmax=vr)
            axs[2 + i, j].set_title(f'aRBC bin {b} — {t} (sag)', fontsize=8, color='white')
    for ax in axs.ravel():
        ax.set_facecolor('black'); ax.set_xticks([]); ax.set_yticks([])
    r, n = res['raw'], res['new_container']
    fig.suptitle(f'{key} · static fine corr aRBC {r["static_fine_corr07_rbc"]:.2f} → {n["static_fine_corr07_rbc"]:.2f} · '
                 f'negatives {r["oow_pct_bin7"]:.1f} → {n["oow_pct_bin7"]:.1f} % · new container vs offline true-b: corr {res["corr_aRBC_new_vs_offline_bin7"]:.3f} · '
                 f'gas change rms {100 * res["gas_relrms_new_vs_raw_bin7"]:.2f} %', fontsize=8, color='white')
    fig.tight_layout(); fig.savefig(FIG / f'{key}_bsmooth_frames.png', dpi=120, facecolor='black'); plt.close(fig)
    return res


if __name__ == '__main__':
    r = run(sys.argv[1])
    print(json.dumps(r, indent=1))
