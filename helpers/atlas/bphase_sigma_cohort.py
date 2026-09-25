#!/usr/bin/env python
"""Cohort sigma sweep for the calcb phase low-pass: 0 (raw b), 1, 2.5, 5, 10 vox, and no b — on every session
that has a raw-b export (Ext tyger_bexport_2026-09-25/<key>/d/output.mrd, fork c8366c3).

Per session and sigma: split the dissolved image with b's phase low-passed at sigma (offline, C44 recipe),
then the static component of aRBC between bins 0 and 7 in three bands (fine < 1.5, mid 1.5-6, coarse > 6 vox)
against the |dis| floor, the out-of-wedge fraction at EI, and the per-voxel RBC SNR proxy. Writes one CSV and
a cohort figure; the decision is read off the medians and the spread.

Usage: bphase_sigma_cohort.py [--sigmas 0,1,2.5,5,10,nob] [--only KEY ...]
Outputs: outputs/bphase_fix_2026-09-25/cohort_sigma.csv, cohort_sigma.png
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
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
import snr_calc as S  # noqa: E402
from bphase_confirm import read_mrd, smooth_phase, resplit, OUT  # noqa: E402
from bphase_sigma import band_static  # noqa: E402

RAW = Path('/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/tyger_bexport_2026-09-25')
TABLE = HERE.parents[1] / 'outputs' / 'te90_2026-09-25' / 'te90_table.csv'


def one(key, sigmas):
    gas, dis, b, meta = read_mrd(RAW / key / 'd' / 'output.mrd')
    if 'rbc_tp_dphi_deg' not in meta:
        # container gate left this session unsplit (stability refit > 10 deg): the dissolved item is the raw
        # complex image; take angle + target from the offline XeCS re-split row (same fit, offline gate)
        import scipy.io as sio
        r = sio.loadmat(f'/Volumes/HoomHamExt/AIkill_Dynamic/{key}/d/recon_resplit_xecs.mat',
                        variable_names=['dphi_new_deg', 'target'])
        phi = np.radians(float(r['dphi_new_deg'].ravel()[0])); target = float(r['target'].ravel()[0])
        # raw complex -> rotate into the stored-basis convention via one global sweep so that 'sigma 0' is comparable
        from resplit_rbctp import sweep_phase, split_maps
        z0 = np.asarray(dis); dis = np.zeros_like(z0)
        for k in range(gas.shape[0]):
            dis[k] = z0[k]  # placeholder; replaced below once masks exist
        meta = dict(meta, rbc_tp_dphi_deg=str(np.degrees(phi)), rbc_tp_target=str(target), _raw_unsplit='1')
    phi = np.radians(float(meta['rbc_tp_dphi_deg'])); target = float(meta['rbc_tp_target']); nb = gas.shape[0]
    gm = gas.mean(0); corners = S.corner_box(gas.shape[1:])
    bg = corners & ~ndimage.binary_dilation(S.lung_mask(gm, gm[corners].std()), iterations=S.EXCL_DILATE)
    masks = [S.lung_mask(gas[k], gas[k][bg].std()) for k in range(nb)]
    if masks[0].sum() < 500 or masks[7].sum() < 500:
        return []
    many = np.any(masks, 0); w = ndimage.binary_dilation(many, iterations=3).astype(float)
    b_in = int(np.argmax([m.sum() for m in masks]))
    if meta.get('_raw_unsplit'):
        # raw complex image: split it once with the raw b (sigma 0 reference) so all variants start alike
        R0, T0 = resplit(np.asarray(dis), masks, phi, target)
        dis = R0 + 1j * T0
    z = dis.real * np.exp(1j * phi) + dis.imag; Z = np.abs(z); m07 = masks[0] & masks[7]
    bu = b[0] / np.maximum(np.abs(b[0]), 1e-12)
    rows = [{'key': key, 'sigma': 'floor', **band_static(Z, m07), 'phi_k0_deg': float(np.degrees(phi))}]
    for s in sigmas:
        if s == 0:
            R, T = dis.real, dis.imag
        elif s == 'nob':
            R, T = resplit(z * np.conj(bu)[None], masks, phi, target)
        else:
            d = np.where(w > 0, np.angle(np.exp(1j * (np.angle(bu) - smooth_phase(bu, w, float(s))))), 0.0)
            R, T = resplit(z * np.exp(-1j * d)[None], masks, phi, target)
        mk = masks[b_in]
        row = {'key': key, 'sigma': str(s), **band_static(R, m07), 'phi_k0_deg': float(np.degrees(phi)),
               'oow_pct_ei': 100 * float(((R[b_in] < 0) | (T[b_in] < 0))[mk].mean()),
               'negRBC_pct_ei': 100 * float((R[b_in] < 0)[mk].mean()),
               'lung_rbc_over_bgsd': float(np.median(np.abs(R[b_in][mk])) / R[b_in][bg].std()),
               'corr_rbc_vs_raw': float(np.corrcoef(R[b_in][mk], dis.real[b_in][mk])[0, 1])}
        row.update({'tm_' + k: v for k, v in band_static(T, m07).items()})
        rows.append(row)
    return rows


def main():
    a = sys.argv[1:]
    sig = a[a.index('--sigmas') + 1].split(',') if '--sigmas' in a else ['0', '1', '2.5', '5', '10', 'nob']
    sig = [s if s == 'nob' else float(s) for s in sig]
    t = pd.read_csv(TABLE); keys = sorted(t[t.atlas_set.eq('yes')].key)
    if '--only' in a:
        keys = a[a.index('--only') + 1:]
    keys = [k for k in keys if (RAW / k / 'd' / 'output.mrd').exists()]
    rows = []
    for k in keys:
        try:
            r = one(k, sig); rows += r
            print(f'[sigma] {k}: ' + ' | '.join(f'{x["sigma"]}: fine {x["fine<1.5"]:.2f} oow {x.get("oow_pct_ei", np.nan):.1f}' for x in r), flush=True)
        except Exception as e:
            print(f'[sigma] {k}: FAILED {e!r}', flush=True)
    d = pd.DataFrame(rows); d.to_csv(OUT / 'cohort_sigma.csv', index=False)
    # summary
    med = d.groupby('sigma')[['fine<1.5', 'mid1.5-6', 'coarse>6', 'oow_pct_ei', 'negRBC_pct_ei', 'lung_rbc_over_bgsd', 'corr_rbc_vs_raw']].median()
    print('\nmedians over', d.key.nunique(), 'sessions'); print(med.round(3).to_string())
    # excess over floor per session
    fl = d[d.sigma == 'floor'].set_index('key')
    ex = d[d.sigma != 'floor'].copy()
    for band in ('fine<1.5', 'mid1.5-6', 'coarse>6'):
        ex[band + '_excess'] = ex[band].values - fl.loc[ex.key, band].values
    print('\nmedian excess of aRBC static corr over the |dis| floor (>0 = phase-carried static structure):')
    print(ex.groupby('sigma')[[c for c in ex.columns if c.endswith('_excess')]].median().round(3).to_string())
    print('\nfraction of sessions with fine excess > 0.1:'); print(ex.groupby('sigma').apply(lambda g: (g['fine<1.5_excess'] > 0.1).mean()).round(2).to_string())
    order = [s for s in ['0', '1.0', '2.5', '5.0', '10.0', 'nob'] if s in set(ex.sigma)]
    fig, axs = plt.subplots(1, 4, figsize=(17, 4.2))
    for j, (col, title) in enumerate([('fine<1.5_excess', 'fine (<1.5 vox) static excess over |dis| floor'),
                                      ('coarse>6_excess', 'coarse (>6 vox) static excess over floor'),
                                      ('oow_pct_ei', 'out-of-wedge lung voxels at EI [%]'),
                                      ('corr_rbc_vs_raw', 'aRBC corr vs raw-b map (how much changed)')]):
        axs[j].boxplot([ex[ex.sigma == s][col].dropna() for s in order], labels=[{'0': 'raw b', 'nob': 'no b'}.get(s, f'σ={s}') for s in order], showfliers=False)
        axs[j].set_title(title, fontsize=9); axs[j].axhline(0, color='grey', lw=0.6)
    fig.suptitle(f'calcb phase low-pass: cohort sweep, {d.key.nunique()} sessions — bphase_sigma_cohort.py 2026-09-25', fontsize=10)
    fig.tight_layout(); fig.savefig(OUT / 'cohort_sigma.png', dpi=140)


if __name__ == '__main__':
    main()
