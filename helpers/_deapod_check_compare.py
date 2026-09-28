"""Scratch (s14): Tyger rerun with fork ee3c91f (deapod) vs the offline tree (b44 / A) vs b44, per variable.
usage: python _deapod_check_compare.py KEY ... -> prints + outputs/tune_gas/deapod_check_2026-09-28.md"""
import os
import sys

import numpy as np
import scipy.io as sio

TY = '/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/tyger_deapod_check_2026-09-28'
OFF = '/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_b44d'
B44 = '/Volumes/HoomHamExt/AIkill_Dynamic_b44'
IM = ('gas_phase', 'gas_phase_magnitude', 'dissolved_phase_real', 'dissolved_phase_imag', 'dissolved_phase_magnitude')
lines = ['# Deapod check 2026-09-28: Tyger ee3c91f vs offline tree (b44 / A)\n',
         '| key | variable | relerr Tyger vs offline | relerr Tyger vs b44 (undivided) | max abs diff / max |\n|---|---|---|---|---|\n']
for key in sys.argv[1:]:
    t = sio.loadmat(os.path.join(TY, key, 'd', 'recon.mat'))
    o = sio.loadmat(os.path.join(OFF, key, 'd', 'recon.mat'))
    b = sio.loadmat(os.path.join(B44, key, 'd', 'recon.mat'))
    print(f'== {key}: Tyger deapod flag {np.ravel(t.get("deapod", [-1]))[0]}, offline flag {np.ravel(o.get("deapod", [-1]))[0]}')
    for k in IM:
        if k not in t:
            continue
        re_o = np.linalg.norm(t[k] - o[k]) / np.linalg.norm(o[k])
        re_b = np.linalg.norm(t[k] - b[k]) / np.linalg.norm(b[k])
        mx = np.abs(t[k] - o[k]).max() / np.abs(o[k]).max()
        print(f'  {k:26s} vs offline {re_o:.2e}   vs b44 {re_b:.2e}   max {mx:.1e}')
        lines.append(f'| {key[-5:]} | {k} | {re_o:.2e} | {re_b:.2e} | {mx:.1e} |\n')
    for k in ('calcb_b',):
        re = np.linalg.norm(t[k] - b[k]) / np.linalg.norm(b[k])
        print(f'  {k:26s} vs b44 {re:.2e}')
        lines.append(f'| {key[-5:]} | {k} (unchanged by design) | — | {re:.2e} | |\n')
    nv = np.nanmax(np.abs(np.ravel(t['nav_volume']) - np.ravel(b['nav_volume'])))
    ph_t, ph_b = str(np.ravel(t.get('rbc_tp_ph_rad', ['']))[0]), str(np.ravel(b.get('rbc_tp_ph_rad', ['']))[0])
    print(f'  nav_volume max |Δ| {nv:.2e};  rbc_tp_ph identical: {ph_t == ph_b}')
    if ph_t != ph_b:
        print('   tyger', ph_t, '\n   b44  ', ph_b)
    lines.append(f'| {key[-5:]} | nav_volume max abs diff | {nv:.2e} | | |\n| {key[-5:]} | rbc_tp_ph_rad identical | {ph_t == ph_b} | | |\n')
open(os.path.expanduser('~/Hooman/Work/Codes/2026_ASAP_Recon/workspace/outputs/tune_gas/deapod_check_2026-09-28.md'), 'w').write(''.join(lines))
