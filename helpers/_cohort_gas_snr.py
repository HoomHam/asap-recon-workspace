"""Scratch (s14): real-data gas corner SNR (lung mean / far-corner sd of the production bin image, median over bins)
for every b44 production session, in the same units as tune_gas.py's real_snr. -> outputs/tune_gas/cohort_gas_snr.csv"""
import csv, glob, os, sys
import numpy as np
from scipy.io import loadmat
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pf_test import lung_masks
rows = []
for p in sorted(glob.glob('/Volumes/HoomHamExt/AIkill_Dynamic_b44/*/d/recon.mat')):
    key = p.split('/')[-3]
    try:
        g = np.asarray(loadmat(p, variable_names=['gas_phase'])['gas_phase'], dtype=np.float64)
    except Exception as e:
        print(key, 'FAIL', e); continue
    m = lung_masks(np.abs(g).mean(0))
    snr = [g[b][m['lung']].mean() / g[b][m['corner']].std() for b in range(g.shape[0])]
    rows.append(dict(key=key, snr_med=float(np.median(snr)), snr_min=float(np.min(snr)), snr_max=float(np.max(snr)), lung_vox=int(m['lung'].sum())))
    print(f"{key}: median {rows[-1]['snr_med']:.1f} (min {rows[-1]['snr_min']:.1f}, max {rows[-1]['snr_max']:.1f})", flush=True)
out = os.path.expanduser('~/Hooman/Work/Codes/2026_ASAP_Recon/workspace/outputs/tune_gas/cohort_gas_snr.csv')
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
s = np.array([r['snr_med'] for r in rows])
print(f'n={len(s)} median {np.median(s):.1f}  quartiles {np.percentile(s,25):.1f} / {np.percentile(s,75):.1f}  '
      f'<13: {np.sum(s<13)}  <10: {np.sum(s<10)}  <7: {np.sum(s<7)}  <5: {np.sum(s<5)}')
