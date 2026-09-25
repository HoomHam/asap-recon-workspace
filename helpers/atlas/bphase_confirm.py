#!/usr/bin/env python
"""Direct confirmation of F79 with the exported calcb b (fork c8366c3, item `calcb_b`).

Compares, inside the lung:
  1. the fine-scale (< 1.5 vox) phase of b   vs   the static fine-scale dissolved phase (bin-mean of z_rot)
     -> if F79 is right they are the same field (corr ~ 1, slope ~ 1 in the wedge-visible part)
  2. the proxy fix (C43, delta from the bin mean) vs the TRUE fix (divide out b's fine phase) on the
     static-texture metric and the out-of-wedge fraction
  3. what the true fix leaves: the static fine corr should fall to the |dis| level

Usage: bphase_confirm.py <run_root> <key> [--sigma 2.5]
Outputs: outputs/bphase_fix_2026-09-25/<key>_confirm.json, fig/<key>_confirm.png
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
import mrd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import snr_calc as S  # noqa: E402
from resplit_rbctp import sweep_phase, split_maps  # noqa: E402

WS = HERE.parents[1]
OUT = WS / 'outputs' / 'bphase_fix_2026-09-25'
FIG = OUT / 'fig'
fine = lambda v, s=1.5: v - ndimage.gaussian_filter(v, s)


def read_mrd(path):
    gas = dis = b = None; meta = {}
    with mrd.BinaryMrdReader(str(path)) as r:
        r.read_header()
        for it in r.read_data():
            if isinstance(it, mrd.StreamItem.NdArrayFloat) and it.value.meta.get('gas_phase_image'):
                gas = np.asarray(it.value.data, dtype=np.float64)
            elif isinstance(it, mrd.StreamItem.NdArrayComplexFloat):
                if it.value.meta.get('dissolved_phase_image'):
                    dis = np.asarray(it.value.data, dtype=np.complex128)
                    meta = {k: str(v[0].value) for k, v in it.value.meta.items()}
                elif it.value.meta.get('calcb_b'):
                    b = np.asarray(it.value.data, dtype=np.complex128)
    return gas, dis, b, meta


def smooth_phase(c, w, sigma):
    num = ndimage.gaussian_filter(c.real * w, sigma) + 1j * ndimage.gaussian_filter(c.imag * w, sigma)
    return np.angle(num / (ndimage.gaussian_filter(w, sigma) + 1e-9))


def resplit(z, masks, phi, target):
    nb = z.shape[0]; R = np.zeros(z.shape); T = np.zeros(z.shape)
    for k in range(nb):
        mk = masks[k]
        if not mk.any():
            continue
        p, *_ = sweep_phase(z[k][mk].sum(), phi, target)
        if np.isnan(p):
            continue
        R[k], T[k] = split_maps(z[k] * np.exp(1j * p), phi)
    return R, T


def metrics(R, T, Z, masks, bg):
    m07 = masks[0] & masks[7]
    return {'static_fine_corr07_rbc': float(np.corrcoef(fine(R[0])[m07], fine(R[7])[m07])[0, 1]),
            'static_fine_corr07_tm': float(np.corrcoef(fine(T[0])[m07], fine(T[7])[m07])[0, 1]),
            'static_fine_corr07_dis': float(np.corrcoef(fine(Z[0])[m07], fine(Z[7])[m07])[0, 1]),
            'oow_pct_bin7': 100 * float(((R[7] < 0) | (T[7] < 0))[masks[7]].mean()),
            'lung_rbc_over_bgsd_bin7': float(np.median(np.abs(R[7][masks[7]])) / R[7][bg].std())}


def run(root, key, sigma):
    gas, dis, b, meta = read_mrd(Path(root) / key / 'd' / 'output.mrd')
    assert b is not None, 'no calcb_b item'
    b = b[0]                                                    # nch = 1
    phi = np.radians(float(meta['rbc_tp_dphi_deg'])); target = float(meta['rbc_tp_target'])
    nb = gas.shape[0]
    gm = gas.mean(0); corners = S.corner_box(gas.shape[1:])
    bg = corners & ~ndimage.binary_dilation(S.lung_mask(gm, gm[corners].std()), iterations=S.EXCL_DILATE)
    masks = [S.lung_mask(gas[k], gas[k][bg].std()) for k in range(nb)]
    many = np.any(masks, 0); w = ndimage.binary_dilation(many, iterations=3).astype(float)
    R0, T0 = dis.real, dis.imag
    z_rot = R0 * np.exp(1j * phi) + T0
    Z = np.abs(z_rot)
    # 1. fields
    zm = z_rot.mean(0)
    delta_proxy = np.where(w > 0, np.angle(np.exp(1j * (np.angle(zm) - smooth_phase(zm, w, sigma)))), 0.0)
    bu = b / np.maximum(np.abs(b), 1e-12)
    delta_b = np.where(w > 0, np.angle(np.exp(1j * (np.angle(bu) - smooth_phase(bu, w, sigma)))), 0.0)
    support = np.abs(b) > 0  # |b| = 1 for nch=1 -> use gm-based proxy for the support instead
    supp = gm > 10 * 1.25 * gm[bg].std()
    res = {'key': key, 'sigma': sigma, 'phi_k0_deg': float(np.degrees(phi)),
           'b_abs_unique': [float(np.min(np.abs(b))), float(np.max(np.abs(b)))],
           'corr_deltaproxy_deltab_lung': float(np.corrcoef(delta_proxy[many], delta_b[many])[0, 1]),
           'corr_deltaproxy_deltab_support': float(np.corrcoef(delta_proxy[many & supp], delta_b[many & supp])[0, 1]),
           'corr_deltaproxy_deltab_outside_support': float(np.corrcoef(delta_proxy[many & ~supp], delta_b[many & ~supp])[0, 1]) if (many & ~supp).sum() > 200 else np.nan,
           'slope_proxy_on_b_lung': float(np.polyfit(delta_b[many], delta_proxy[many], 1)[0]),
           'delta_b_lung_deg': {'median_abs': float(np.degrees(np.median(np.abs(delta_b[many])))),
                                'p90_abs': float(np.degrees(np.percentile(np.abs(delta_b[many]), 90)))},
           'delta_proxy_lung_deg': {'median_abs': float(np.degrees(np.median(np.abs(delta_proxy[many])))),
                                    'p90_abs': float(np.degrees(np.percentile(np.abs(delta_proxy[many]), 90)))},
           'fine_phase_b_static_vs_fine_aRBC_corr': float(np.corrcoef(
               (np.sin(delta_b) * Z[7])[masks[7]], fine(R0[7])[masks[7]])[0, 1])}
    # 2. fixes
    res['before'] = metrics(R0, T0, Z, masks, bg)
    Rp, Tp = resplit(z_rot * np.exp(-1j * delta_proxy)[None], masks, phi, target); res['proxy_fix'] = metrics(Rp, Tp, Z, masks, bg)
    Rb, Tb = resplit(z_rot * np.exp(-1j * delta_b)[None], masks, phi, target); res['true_b_fix'] = metrics(Rb, Tb, Z, masks, bg)
    # full smooth-b variant (2steve/06 §4 style): replace b's phase inside the lung by its sigma-smoothed version = same as true_b_fix here
    res['corr_rbc_proxy_vs_trueb_bin7'] = float(np.corrcoef(Rp[7][masks[7]], Rb[7][masks[7]])[0, 1])
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(exist_ok=True)
    json.dump(res, open(OUT / f'{key}_confirm.json', 'w'), indent=1)
    ys = np.where(masks[7].any(axis=(0, 2)))[0]; y = int(ys[len(ys) // 2])
    sl = lambda v: v[:, y, :][::-1][15:85, 10:90]
    fig, axs = plt.subplots(2, 4, figsize=(13, 6.8), facecolor='black')
    vr = np.percentile(np.maximum(R0[7], 0)[masks[7]], 99.5)
    panels = [(np.degrees(delta_b), 'twilight', -30, 30, 'fine phase of b [deg]'),
              (np.degrees(delta_proxy), 'twilight', -30, 30, 'proxy δ from bin-mean z [deg]'),
              (np.degrees(np.angle(bu)), 'twilight', -180, 180, 'phase of b (full)'),
              (supp.astype(float), 'gray', 0, 1, 'calcb support proxy (|F_avg| > 10·noise)'),
              (np.maximum(R0[7], 0), 'gray', 0, vr, 'aRBC bin 7 container'),
              (np.maximum(Rp[7], 0), 'gray', 0, vr, 'aRBC bin 7 proxy fix'),
              (np.maximum(Rb[7], 0), 'gray', 0, vr, 'aRBC bin 7 true-b fix'),
              (fine(R0[7]), 'gray', -vr / 3, vr / 3, 'fine aRBC bin 7 before')]
    for ax, (img, cm, lo, hi, t) in zip(axs.ravel(), panels):
        ax.imshow(sl(img), cmap=cm, vmin=lo, vmax=hi); ax.set_facecolor('black'); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(t, fontsize=8, color='white')
    fig.suptitle(f'{key} · corr(proxy δ, fine phase of b) in lung {res["corr_deltaproxy_deltab_lung"]:.2f} (slope {res["slope_proxy_on_b_lung"]:.2f}) · '
                 f'static fine corr aRBC: before {res["before"]["static_fine_corr07_rbc"]:.2f}, proxy {res["proxy_fix"]["static_fine_corr07_rbc"]:.2f}, '
                 f'true-b {res["true_b_fix"]["static_fine_corr07_rbc"]:.2f} (|dis| {res["before"]["static_fine_corr07_dis"]:.2f})', fontsize=8, color='white')
    fig.tight_layout(); fig.savefig(FIG / f'{key}_confirm.png', dpi=120, facecolor='black'); plt.close(fig)
    return res


if __name__ == '__main__':
    a = sys.argv[1:]
    sigma = float(a[a.index('--sigma') + 1]) if '--sigma' in a else 2.5
    r = run(a[0], a[1], sigma)
    print(json.dumps({k: v for k, v in r.items() if k not in ('key',)}, indent=1))
