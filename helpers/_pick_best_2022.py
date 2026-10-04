"""Scratch: per 2022 session keep the better of two Tyger recons (all coils vs weak coil dropped) under the
canonical folder name AIkill_Dynamic_2022/<date>_<id>; the loser keeps a suffix (_allcoil / _autodrop).
Metric = median over 16 bins of gas SNR (lung-mask mean / far-background sd, same as the 2026-10-03 table).
Hooman 2026-10-04: "rerun all 8 with all coils, keep the better". HH 2022-10-26 folders are NOT renamed (paths
already given to XeCS). Usage: _pick_best_2022.py [--apply]  (default dry run). Log -> _logs/_pick_best.tsv
"""
import pathlib
import sys

import numpy as np
import scipy.io as sio
from scipy.ndimage import binary_dilation, gaussian_filter

O = pathlib.Path('/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_2022')
# canonical key -> (current canonical content label, alternative folder suffix, alternative label)
PAIRS = {k: ('autodrop', '_allcoil', 'allcoil') for k in
         ['2022-09-27_021JM', '2022-09-28_017AK', '2022-09-30_029CK', '2022-10-06_019WR', '2022-10-10_001SW',
          '2022-10-11_011BA', '2022-10-19_023DB', '2022-10-27_000HH']}
PAIRS['2022-09-21_027GB'] = ('autodrop', '_8coil', 'allcoil')
PAIRS['2022-11-09_030DN'] = ('allcoil', '_autodrop', 'autodrop')
PAIRS['2022-11-08_018BI'] = ('allcoil', '_autodrop', 'autodrop')


def snr(folder):
    g = np.abs(sio.loadmat(folder / 'd' / 'recon.mat', variable_names=['gas_phase'])['gas_phase'])
    m = gaussian_filter(g.mean(0), 1.5)
    lung = m > np.median(m) + 0.35 * (np.percentile(m, 99.5) - np.median(m))
    far = ~binary_dilation(lung, iterations=6)
    return float(np.median([g[b][lung].mean() / g[b][far].std() for b in range(g.shape[0])]))


def main():
    apply = '--apply' in sys.argv
    log = []
    for key, (cur_lab, alt_suf, alt_lab) in PAIRS.items():
        cur, alt = O / key, O / (key + alt_suf)
        if not (cur / 'd' / 'recon.mat').exists() or not (alt / 'd' / 'recon.mat').exists():
            print(f'{key}: missing one side, skip'); continue
        sc, sa = snr(cur), snr(alt)
        win = alt_lab if sa > sc else cur_lab
        print(f'{key}: {cur_lab} {sc:.1f} vs {alt_lab} {sa:.1f} -> keep {win}' + ('' if sa > sc else ' (no change)'))
        log.append(f'{key}\t{cur_lab}\t{sc:.2f}\t{alt_lab}\t{sa:.2f}\t{win}')
        if sa > sc and apply:
            tmp = O / (key + '__swap')
            cur.rename(tmp)
            alt.rename(cur)
            tmp.rename(O / (key + '_' + cur_lab))
            (cur / 'CHOSEN.txt').write_text(f'kept {alt_lab} (median gas SNR {sa:.2f}) over {cur_lab} ({sc:.2f}); '
                                             f'loser in {key}_{cur_lab}\n')
        elif apply:
            (cur / 'CHOSEN.txt').write_text(f'kept {cur_lab} (median gas SNR {sc:.2f}) over {alt_lab} ({sa:.2f}); '
                                             f'loser in {alt.name}\n')
    if apply:
        (O / '_logs' / '_pick_best.tsv').write_text('key\tcanon_label\tsnr\talt_label\tsnr\tkept\n' + '\n'.join(log) + '\n')
        print('applied; log', O / '_logs' / '_pick_best.tsv')


if __name__ == '__main__':
    main()
