#!/usr/bin/env python
"""Is the fine structure in the aRBC maps signal, noise, or a fixed artifact?  Three reproducibility tests.

Run 2026-09-25 inline (ledger session 11); saved here verbatim for re-running. All tests use the
XeCS re-split maps (Ext AIkill_Dynamic/<key>/d/recon_resplit_xecs.mat) and the atlas set
(outputs/te90_2026-09-25/te90_table.csv atlas_set == yes). "hp" = map minus Gaussian(6 vox) = structure
finer than ~6 vox, envelope removed. Masks as snr_calc (gas > 5 sigma at the bin).

  mirror   : DIAPHRAGM bins k and 15-k hold the same lung volume but independent interleaves.
             corr(bin k, bin 15-k) inside the common lung mask, pairs (3,12),(4,11),(5,10) averaged.
             Noise does not reproduce; signal and position-locked artifacts do.
  subjects : different subjects, insp bin, same scanner voxels (lung overlap >= 50 %).
             A scanner/sequence-fixed pattern reproduces; physiology does not.
  dates    : same subject, different dates, insp bin, rigid 3D shift found on the gas image.
             Physiology reproduces; a per-session b-reference error does not.
Outputs: outputs/split_trust_2026-09-25/{mirror_bin_repro,cross_subject_fixedvoxel,cross_date_same_subject}.csv
Usage: rbc_repro_tests.py mirror|subjects|dates
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
from scipy import ndimage

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import snr_calc as S  # noqa: E402

SRC = Path('/Volumes/HoomHamExt/AIkill_Dynamic')
WS = HERE.parents[1]
OUT = WS / 'outputs' / 'split_trust_2026-09-25'
TABLE = WS / 'outputs' / 'te90_2026-09-25' / 'te90_table.csv'
hp = lambda v: v - ndimage.gaussian_filter(v, 6)


def atlas_keys():
    t = pd.read_csv(TABLE)
    return sorted(t[t.atlas_set.eq('yes')].key)


def load(key, bins=None):
    p = SRC / key / 'd'
    gas = sio.loadmat(p / 'recon.mat', variable_names=['gas_phase'])['gas_phase'].astype(float)
    r = sio.loadmat(p / 'recon_resplit_xecs.mat', variable_names=['aRBC', 'aTP', 'dphi_new_deg', 'insp_bin'])
    phi = np.radians(r['dphi_new_deg'].ravel()[0]); bi = int(r['insp_bin'].ravel()[0])
    R, T = r['aRBC'].astype(float), r['aTP'].astype(float)
    Z = np.abs(R * np.exp(1j * phi) + T)
    corners = S.corner_box(gas.shape[1:]); gm = gas.mean(0)
    bg = corners & ~ndimage.binary_dilation(S.lung_mask(gm, gm[corners].std()), iterations=S.EXCL_DILATE)
    mask = lambda b: S.lung_mask(gas[b], gas[b][bg].std())
    return dict(gas=gas, R=R, T=T, Z=Z, phi=phi, bi=bi, mask=mask)


def roll3(v, s):
    return np.roll(np.roll(np.roll(v, s[0], 0), s[1], 1), s[2], 2)


def test_mirror():
    rows = []
    for k in atlas_keys():
        d = load(k); out = {'key': k}
        acc = {n: [] for n in ('gas', 'dis', 'rbc', 'tm', 'rbc_hp', 'tm_hp')}
        for a, b in ((3, 12), (4, 11), (5, 10)):
            m = d['mask'](a) & d['mask'](b)
            if m.sum() < 500:
                continue
            for n, A, B in (('gas', d['gas'][a], d['gas'][b]), ('dis', d['Z'][a], d['Z'][b]),
                            ('rbc', d['R'][a], d['R'][b]), ('tm', d['T'][a], d['T'][b]),
                            ('rbc_hp', hp(d['R'][a]), hp(d['R'][b])), ('tm_hp', hp(d['T'][a]), hp(d['T'][b]))):
                acc[n].append(np.corrcoef(A[m], B[m])[0, 1])
        out.update({n: float(np.mean(v)) if v else np.nan for n, v in acc.items()})
        rows.append(out)
    t = pd.DataFrame(rows); t.to_csv(OUT / 'mirror_bin_repro.csv', index=False)
    print(t[['gas', 'dis', 'tm', 'rbc', 'tm_hp', 'rbc_hp']].median().round(2).to_dict())


def _insp_maps(k):
    d = load(k); b = d['bi']
    return dict(mk=d['mask'](b), gas=d['gas'][b], dis=d['Z'][b], rbc=d['R'][b], tm=d['T'][b],
                gas_hp=hp(d['gas'][b]), dis_hp=hp(d['Z'][b]), rbc_hp=hp(d['R'][b]), tm_hp=hp(d['T'][b]))


def test_subjects():
    tr = pd.read_csv(OUT / 'split_trust.csv')
    keys = sorted(tr[tr.grade.isin(['A', 'B'])].key)
    cache = {k: _insp_maps(k) for k in keys}
    rows = []
    for a, b in itertools.combinations(keys, 2):
        if a[11:] == b[11:]:
            continue
        m = cache[a]['mk'] & cache[b]['mk']
        if m.sum() < 0.5 * min(cache[a]['mk'].sum(), cache[b]['mk'].sum()) or m.sum() < 3000:
            continue
        row = {'a': a, 'b': b, 'n': int(m.sum())}
        for n in ('gas', 'dis', 'tm', 'rbc', 'gas_hp', 'dis_hp', 'tm_hp', 'rbc_hp'):
            row[n] = np.corrcoef(cache[a][n][m], cache[b][n][m])[0, 1]
        rows.append(row)
    d = pd.DataFrame(rows); d.to_csv(OUT / 'cross_subject_fixedvoxel.csv', index=False)
    print(len(d), 'pairs', d[['gas', 'dis', 'tm', 'rbc', 'gas_hp', 'dis_hp', 'tm_hp', 'rbc_hp']].median().round(3).to_dict())


def test_dates():
    keys = atlas_keys()
    by = pd.Series(keys).groupby(pd.Series(keys).str[11:]).apply(sorted)
    pairs = [(a, b) for ks in by if len(ks) > 1 for a, b in itertools.combinations(ks, 2)]
    cache = {k: _insp_maps(k) for k in sorted(set(sum(map(list, pairs), [])))}
    rows = []
    for a, b in pairs:
        A, B = cache[a], cache[b]
        best, bs = -9, (0, 0, 0)
        grid = [(dz, dy, dx) for dz in range(-10, 11, 2) for dy in range(-8, 9, 2) for dx in range(-8, 9, 2)]
        for s in grid:
            m = A['mk'] & roll3(B['mk'], s)
            if m.sum() < 2000:
                continue
            c = np.corrcoef(A['gas'][m], roll3(B['gas'], s)[m])[0, 1]
            if c > best:
                best, bs = c, s
        for s in [(bs[0] + i, bs[1] + j, bs[2] + l) for i in (-1, 0, 1) for j in (-1, 0, 1) for l in (-1, 0, 1)]:
            m = A['mk'] & roll3(B['mk'], s)
            c = np.corrcoef(A['gas'][m], roll3(B['gas'], s)[m])[0, 1]
            if c > best:
                best, bs = c, s
        m = A['mk'] & roll3(B['mk'], bs)
        row = {'a': a, 'b': b, 'shift': bs, 'overlap': int(m.sum()), 'gas_align': round(best, 3)}
        for n in ('dis', 'tm', 'rbc', 'gas_hp', 'dis_hp', 'tm_hp', 'rbc_hp'):
            row[n] = round(np.corrcoef(A[n][m], roll3(B[n], bs)[m])[0, 1], 3)
        rows.append(row)
    d = pd.DataFrame(rows); d.to_csv(OUT / 'cross_date_same_subject.csv', index=False)
    print(d.to_string(index=False)); print(d[['gas_align', 'dis', 'tm', 'rbc', 'tm_hp', 'rbc_hp']].median().round(2).to_dict())


if __name__ == '__main__':
    {'mirror': test_mirror, 'subjects': test_subjects, 'dates': test_dates}[sys.argv[1]]()
