"""Scratch (s14): summarise the three pre-fork tests from tune_gas.py runs.
  1 soft-bin stability  <key>_stab_bd{0.5,1,4}, <key>_stab_nb{8,24}  (+ baseline = <key> main run, bd 2 / nb 16)
  2 lower SNR           <key>_snr, <key>_snr_nm{1.5,2,3,4}
  3 registration        <key>_reg (+ proxies in every run)
-> outputs/tune_gas/tests_2026-09-28.md"""
import csv
import os

import numpy as np

O = os.path.expanduser('~/Hooman/Work/Codes/2026_ASAP_Recon/workspace/outputs/tune_gas')
KEYS = ['2024-11-13_025JC', '2024-10-02_011CN']


def load(d):
    p = os.path.join(O, d, 'metrics.csv')
    return list(csv.DictReader(open(p))) if os.path.exists(p) else None


def agg(rows, n, dea, col):
    v = [float(r[col]) for r in rows if r['setting'] == n and r['deapod'] == str(dea) and r.get(col, '') not in ('', None)]
    return float(np.mean(v)) if v else np.nan


out = ['# Pre-fork tests 2026-09-28 (tune_gas.py, 025JC + 011CN)\n',
       'Total = in-lung rms error vs the known object / lung mean, mean of 4 bins. rel = / production (MS 240, LB 300, no deapod) of the SAME run. d = deapod.\n']

# ---------------- 1 stability
out.append('\n## 1. Soft-bin stability of the ranking\n')
cfgs = [('bd 2, nb 16 (prod)', ''), ('bd 0.5', '_stab_bd0.5'), ('bd 1', '_stab_bd1'), ('bd 4', '_stab_bd4'),
        ('nb 8', '_stab_nb8'), ('nb 24', '_stab_nb24')]
cols = [('prod', 1), ('MS200', 0), ('MS200', 1), ('MS160', 1), ('gplb450', 1), ('gplb200', 1)]
out.append('| subject | binning | prod total | noise | alias | ' + ' | '.join(f'{n}{"+d" if d else ""} rel' for n, d in cols) + ' | best |\n')
out.append('|---|---|---|---|---|' + '---|' * len(cols) + '---|\n')
for k in KEYS:
    for lab, suf in cfgs:
        rows = load(k + suf)
        if rows is None:
            out.append(f'| {k[-5:]} | {lab} | (missing) |\n')
            continue
        p = agg(rows, 'prod', 0, 'total')
        cand = {(n, d): agg(rows, n, d, 'total') for n in ('prod', 'gplb200', 'gplb450', 'MS160', 'MS200') for d in (0, 1)}
        best = min(cand, key=lambda x: cand[x] if np.isfinite(cand[x]) else 9)
        out.append(f'| {k[-5:]} | {lab} | {p:.3f} | {agg(rows, "prod", 0, "noise"):.3f} | {agg(rows, "prod", 0, "alias"):.3f} | '
                   + ' | '.join(f'{agg(rows, n, d, "total") / p:.3f}' for n, d in cols)
                   + f' | {best[0]}{"+d" if best[1] else ""} |\n')

# ---------------- 2 SNR
out.append('\n## 2. Lower SNR: does the optimum move?\n')
nms = [(1, '_snr'), (1.5, '_snr_nm1.5'), (2, '_snr_nm2'), (3, '_snr_nm3'), (4, '_snr_nm4')]
gl = ['gplb100', 'gplb150', 'gplb200', 'prod', 'gplb450', 'gplb700']
glv = {'gplb100': 100, 'gplb150': 150, 'gplb200': 200, 'prod': 300, 'gplb450': 450, 'gplb700': 700}
out.append('| subject | noise × | SNR at prod (lung mean / noise sd) | real-data corner SNR equiv | ' + ' | '.join(f'LB {glv[g]}' for g in gl)
           + ' | best LB | penalty of LB 300 | MS200+d rel | FWHM best / prod |\n')
out.append('|---|---|---|---|' + '---|' * len(gl) + '---|---|---|---|\n')
snr_rows = []
for k in KEYS:
    for nm, suf in nms:
        rows = load(k + suf)
        if rows is None:
            out.append(f'| {k[-5:]} | {nm} | (missing) |\n')
            continue
        tot = {g: agg(rows, g, 1, 'total') for g in gl}
        best = min(tot, key=tot.get)
        snr = 1 / agg(rows, 'prod', 0, 'noise')
        rs = agg(rows, 'prod', 0, 'real_snr') / nm            # real-data corner SNR scaled as if the data had this noise
        pen = tot['prod'] / tot[best]
        out.append(f'| {k[-5:]} | {nm:g} | {snr:.1f} | {rs:.1f} | ' + ' | '.join(f'{tot[g]:.3f}' for g in gl)
                   + f' | {glv[best]} | ×{pen:.3f} | {agg(rows, "MS200", 1, "total") / agg(rows, "prod", 0, "total"):.3f} | '
                   f'{agg(rows, best, 1, "fwhm_bin0") / agg(rows, "prod", 1, "fwhm_bin0"):.2f} |\n')
        snr_rows.append((k, nm, snr, glv[best], pen))

# ---------------- 3 registration
out.append('\n## 3. Registration proxies (static object: anatomy identical in every bin)\n')
out.append('hp = 2-vox high-pass of (bin recon + noise) inside lung ⊕ 3 vox. pair corr = mean corr between bins (texture that '
           'repeats bin to bin); anat corr = corr with the noise-free all-data image; incons = rms(bin hp − all-data hp) / lung '
           'mean (what a registration would chase); edge = mean |∇| on the lung boundary / lung mean; anchor SNR = edge / incons.\n\n')
regs = ['prod', 'gplb450', 'gplb700', 'gplb0', 'MS200', 'MS160', 'MS200g450', 'MS160g450', 'MS160g700']
out.append('| subject | noise × | setting | FWHM bin0 | noise | hp pair corr | hp anat corr | incons | edge | anchor SNR |\n'
           '|---|---|---|---|---|---|---|---|---|---|\n')
for k in KEYS:
    for nm, d_ in ((1, '_reg'), (2, '_snr_nm2'), (4, '_snr_nm4')):
        rows = load(k + d_)
        if rows is None:
            continue
        for n in regs:
            if not np.isfinite(agg(rows, n, 0, 'total')):
                continue
            out.append(f'| {k[-5:]} | {nm} | {n} | {agg(rows, n, 0, "fwhm_bin0"):.2f} | {agg(rows, n, 0, "noise"):.3f} | '
                       f'{agg(rows, n, 0, "hp_pair_corr"):.3f} | {agg(rows, n, 0, "hp_anat_corr"):.3f} | {agg(rows, n, 0, "hp_incons"):.4f} | '
                       f'{agg(rows, n, 0, "edge"):.4f} | {agg(rows, n, 0, "anchor_snr"):.2f} |\n')

open(os.path.join(O, 'tests_2026-09-28.md'), 'w').write(''.join(out))
print(''.join(out))
