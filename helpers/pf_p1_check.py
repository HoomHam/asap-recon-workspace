"""Legacy dat/PF vs dat/P1 gpdyn: is PF the oldway=0 Hermitian branch of Steve's cudarenorm?

oldway=0 averages each k cell with the conjugate of its mirror and zeroes the DC cell, so the full MS^3
real image sums to zero. Relative to P1 (oldway=1) the PF background is then shifted by -S_bin / MS^3,
S_bin = total object signal of the bin. The test regresses the measured far-bg shift (PF - P1) on that
prediction across the 16 bins: slope ~1 (MS = 240) = DC zeroed = PF is the Hermitian branch.

Also reports geometry (xcorr shift, lung corr), bg sigma ratio and noise corr, gpdyn mtimes, and whether
the Bin_*.mat navigator files of the two runs agree (diaphragm bin jitter).

usage: python pf_p1_check.py [study ...]      (default: every WorkVault Analysis/*/dat/{PF,P1})
out:   workspace/outputs/pf_p1_check/{summary.csv, <study>.txt}
"""
import csv
import glob
import os
import sys
from datetime import datetime

import numpy as np
from numpy.fft import fftn, ifftn
from scipy import ndimage as ndi
from scipy.io import loadmat

ROOT = '/Volumes/WorkVault/Work/Analysis/'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'outputs', 'pf_p1_check')
MS = 240


def masks(mean):
    lung = mean > 0.2 * np.percentile(mean, 99.9)
    lung = ndi.binary_opening(lung, iterations=1)
    lab, nl = ndi.label(lung)
    lung = lab == (1 + np.argmax(ndi.sum(lung, lab, range(1, nl + 1))))
    far = ~ndi.binary_dilation(lung, iterations=16)
    edge = np.zeros_like(far)
    s = max(4, mean.shape[0] // 8)
    edge[:s] = edge[-s:] = True
    edge[:, :s] = edge[:, -s:] = True
    edge[:, :, :s] = edge[:, :, -s:] = True
    return lung, far & edge


def bin_files(study):
    rows = []
    for f in sorted(glob.glob(os.path.join(ROOT, study, 'dat', 'PF', 'Bin_*.mat'))):
        g = f.replace('/PF/', '/P1/')
        if not os.path.exists(g):
            rows.append(f'{os.path.basename(f)}: PF only')
            continue
        ma, mc = loadmat(f), loadmat(g)
        if 'b' not in ma or 'b' not in mc:  # older exports hold t_vector / V_vector only
            rows.append(f'{os.path.basename(f)}: no b key ({[k for k in ma if not k.startswith("__")]})')
            continue
        a, c = ma['b'].ravel().astype(float), mc['b'].ravel().astype(float)
        if a.shape != c.shape:
            rows.append(f'{os.path.basename(f)}: shape {a.shape} vs {c.shape}')
            continue
        if np.allclose(a, c, equal_nan=True):
            rows.append(f'{os.path.basename(f)}: identical')
            continue
        ok = (a >= 0) & (c >= 0)
        d = ((a - c + 0.5) % 1) - 0.5
        q1, q3 = np.percentile(d[ok], [25, 75])
        rows.append(f'{os.path.basename(f)}: DIFFER circ median {np.median(d[ok]):+.4f} IQR {q1:+.4f}..{q3:+.4f} (cycle units)')
    return rows


def check(study):
    d = os.path.join(ROOT, study, 'dat')
    pf = loadmat(os.path.join(d, 'PF', 'gpdyn.mat'))['gpdyn']
    p1 = loadmat(os.path.join(d, 'P1', 'gpdyn.mat'))['gpdyn']
    mt = [datetime.fromtimestamp(os.path.getmtime(os.path.join(d, s, 'gpdyn.mat'))).strftime('%Y-%m-%d %H:%M')
          for s in ('PF', 'P1')]
    if pf.shape != p1.shape:
        return dict(study=study, shape=f'{pf.shape} vs {p1.shape}', mtime_pf=mt[0], mtime_p1=mt[1]), ['shape mismatch']
    nb = pf.shape[0]
    lung, corner = masks(0.5 * (pf.mean(0) + p1.mean(0)))
    xc = np.real(ifftn(fftn(pf.mean(0)) * np.conj(fftn(p1.mean(0)))))
    n = pf.shape[1]
    sh = np.array(np.unravel_index(np.argmax(xc), xc.shape))
    sh = np.where(sh > n // 2, sh - n, sh)
    lines, pred, meas, sdr, bgc, lc = [], [], [], [], [], []
    lines.append('bin | pred dbg = -S/240^3 | meas dbg | bg sd PF/P1 | bg corr | lung corr')
    for b in range(nb):
        bg1 = p1[b][corner].mean()
        S = p1[b].sum() - bg1 * p1[b].size
        pred.append(-S / MS**3)
        meas.append(pf[b][corner].mean() - bg1)
        sdr.append(pf[b][corner].std() / p1[b][corner].std())
        bgc.append(np.corrcoef(pf[b][corner], p1[b][corner])[0, 1])
        lc.append(np.corrcoef(pf[b][lung], p1[b][lung])[0, 1])
        lines.append(f'{b:2d} | {pred[-1]:9.2f} | {meas[-1]:9.2f} | {sdr[-1]:.3f} | {bgc[-1]:+.3f} | {lc[-1]:.4f}')
    pred, meas = np.array(pred), np.array(meas)
    slope, icpt = np.linalg.lstsq(np.vstack([pred, np.ones_like(pred)]).T, meas, rcond=None)[0]
    row = dict(study=study, shape=str(pf.shape), mtime_pf=mt[0], mtime_p1=mt[1],
               shift=' '.join(map(str, sh)), lung_corr=round(float(np.mean(lc)), 4),
               dc_slope=round(float(slope), 3), dc_icpt=round(float(icpt), 2),
               dc_corr=round(float(np.corrcoef(pred, meas)[0, 1]), 3),
               dc_ratio_med=round(float(np.median(meas / pred)), 3),
               bg_sd_ratio=round(float(np.mean(sdr)), 3), bg_corr=round(float(np.mean(bgc)), 3),
               lung_vox=int(lung.sum()))
    lines.append(f'\nDC fit: meas = {slope:.3f}*pred + {icpt:.2f}; corr {row["dc_corr"]}; median ratio {row["dc_ratio_med"]}')
    lines += ['', *bin_files(study)]
    return row, lines


def main():
    os.makedirs(OUT, exist_ok=True)
    studies = sys.argv[1:] or sorted(p.split('/')[-4]  # <ROOT>/<study>/dat/PF/gpdyn.mat
                                     for p in glob.glob(ROOT + '*/dat/PF/gpdyn.mat')
                                     if os.path.exists(p.replace('/PF/', '/P1/')))
    rows = []
    for s in studies:
        try:
            row, lines = check(s)
        except Exception as e:  # keep the cohort run going; record the failure
            row, lines = dict(study=s, shape=f'ERROR {type(e).__name__}: {e}'), [repr(e)]
        rows.append(row)
        with open(os.path.join(OUT, f'{s}.txt'), 'w') as f:
            f.write(f'{s}\n{row}\n\n' + '\n'.join(lines) + '\n')
        print(row, flush=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(OUT, 'summary.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__':
    main()
