"""Scratch: QC contact sheet for the 2026-10-21 leftover recons (saturation, ADC b-sets, BigMac/WorkVault 2022,
early-2022 entry-1/2 calibrations; chunk-split early sessions shown as their *_merged versions). Per output: time-averaged gas magnitude, axial + coronal through the lung centroid,
sampling pattern + calcb fallbacks from tyger.log, crude SNR (lung p95 / corner-cube noise sd).

usage: python3 helpers/_contact_2026_10_21.py [method d|s]
writes: outputs/leftover_recon_2026-10-21/contact_<method>.png + qc_<method>.csv
"""
import csv
import pathlib
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

ROOT = pathlib.Path('/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon')
OUT = pathlib.Path(__file__).resolve().parents[1] / 'outputs' / 'leftover_recon_2026-10-21'
KEYS = {  # key -> tree
    **{k: 'AIkill_Dynamic_2022' for k in [
        '2022-02-02_AW', '2022-02-02_unk-MID16_merged', '2022-02-03_000ID', '2022-02-03_000ID_MID143',
        '2022-02-04_027CB', '2022-02-09_unk-MID38_merged', '2022-02-09_unk-MID74_merged',
        '2022-02-15_unk-MID159_merged', '2022-02-24_000HH_b0', '2022-02-24_000HH_b', '2022-02-25_054SK_merged',
        '2022-03-09_015MG_merged', '2022-03-10_006RL_A_merged', '2022-03-10_006RL_B_merged',
        '2022-03-25_016ID', '2022-03-30_020WC', '2022-04-13_000LL', '2022-04-19_000HH', '2022-04-20_017AK_b0',
        '2022-04-20_017AK_b', '2022-04-21_000HH_b0', '2022-04-21_000HH_b', '2022-04-22_000HH_b0', '2022-04-22_000HH_b',
        '2022-04-25_000LL_b0', '2022-04-25_000LL_b', '2022-04-26_000ID_b0', '2022-04-26_000ID_b',
        '2022-04-29_018BI_b0', '2022-04-29_018BI_b', '2022-05-02_020ST_b0', '2022-05-02_020ST_b',
        '2022-06-02_008LL', '2022-06-02_008LL_dptraj', '2022-06-27_000LL', '2022-06-28_019WR', '2022-07-05_021JM', '2022-07-06_000HH',
        '2022-07-07_001SW', '2022-07-08_000HH', '2022-07-18_025VP', '2022-07-20_023DB', '2022-07-21_014JW',
        '2022-08-12_000LL', '2022-10-31_000LL', '2022-11-04_000HH', '2022-11-07_000HH', '2023-01-24_000LL',
        '2023-11-09_000HH_b0', '2023-11-09_000HH_b']},
    '2023-02-17_000LL': 'AIkill_Dynamic_b44d', '2023-03-02_000LL': 'AIkill_Dynamic_b44d',
}


def load(key, tree, m):
    d = ROOT / tree / key / m
    if not (d / 'recon.mat').exists():
        return None, d
    r = sio.loadmat(d / 'recon.mat', variable_names=['gas_phase_magnitude', 'gas_phase'])
    g = r.get('gas_phase_magnitude', r.get('gas_phase'))
    return np.abs(g).mean(0), d


def main():
    m = sys.argv[1] if len(sys.argv) > 1 else 'd'
    OUT.mkdir(parents=True, exist_ok=True)
    keys = list(KEYS)
    n = len(keys)
    ncol = 8
    nrow = int(np.ceil(n / ncol))
    fig, ax = plt.subplots(2 * nrow, ncol, figsize=(2.1 * ncol, 2.3 * 2 * nrow))
    rows = []
    for j, key in enumerate(keys):
        a0, a1 = ax[2 * (j // ncol), j % ncol], ax[2 * (j // ncol) + 1, j % ncol]
        g, d = load(key, KEYS[key], m)
        log = (d / 'tyger.log').read_text(errors='ignore') if (d / 'tyger.log').exists() else ''
        pat = (re.findall(r'identified sampling pattern (\S+)', log) or ['?'])[0]
        nfb = len(re.findall(r'flat phase fallback', log))
        if g is None:
            a0.set_title(f'{key}\nMISSING {m}', fontsize=6, color='r')
            rows.append(dict(key=key, method=m, status='missing', pattern=pat))
            continue
        c = np.concatenate([g[:10, :10, :10].ravel(), g[:10, -10:, :10].ravel(), g[-10:, :10, -10:].ravel(),
                            g[-10:, -10:, -10:].ravel()])
        noise = c.std() if c.std() > 0 else 1.0
        mask = g > np.median(c) + 5 * noise
        snr = np.percentile(g[mask], 95) / noise if mask.sum() > 500 else 0.0
        ctr = [int(np.round(x)) for x in np.argwhere(mask).mean(0)] if mask.sum() > 500 else [50, 50, 50]
        v = np.percentile(g, 99.7)
        a0.imshow(g[ctr[0]], cmap='gray', vmax=v)
        a1.imshow(g[:, ctr[1]], cmap='gray', vmax=v)
        a0.set_title(f'{key}\n{pat} snr~{snr:.0f}{" fb" + str(nfb) if nfb else ""}', fontsize=6)
        rows.append(dict(key=key, method=m, status='ok', pattern=pat, calcb_fallbacks=nfb, snr_crude=round(snr, 1),
                         lung_vox=int(mask.sum())))
    for a in ax.flat:
        a.axis('off')
    plt.tight_layout()
    fig.savefig(OUT / f'contact_{m}.png', dpi=80)
    with open(OUT / f'qc_{m}.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['key', 'method', 'status', 'pattern', 'calcb_fallbacks', 'snr_crude',
                                           'lung_vox'])
        w.writeheader()
        w.writerows(rows)
    print(OUT / f'contact_{m}.png')
    for r in rows:
        print(r)


if __name__ == '__main__':
    main()
