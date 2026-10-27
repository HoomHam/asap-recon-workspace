"""Scratch (2026-10-26): verify the b44d_sig s/ rerun against the b44d d/ production tree per session.

Checks: recon.mat keys; calcb_b identity (same pooled phase reference -> same absolute scale); per-bin masked
gas sums for s/ vs the SIGNAL rows of outputs/bin_timing_2026-10-26/<key>_bins.csv (k0 order); real/magnitude
ratio per bin (cos(dphi) loss diagnostic); corner-noise mean/std of gas_phase (zero-mean check); dissolved
sum(real+imag) vs |sum z| vs sum|z| per bin.
usage: python helpers/_sig_vs_d_check.py <key> ...
"""
import os
import sys
import numpy as np
import scipy.io as sio

EXT = '/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon'
BT = os.path.join(os.path.dirname(__file__), '..', 'outputs', 'bin_timing_2026-10-26')


def load(p):
    m = sio.loadmat(p)
    return {k: v for k, v in m.items() if not k.startswith('__')}


def main(key):
    s = load(f'{EXT}/AIkill_Dynamic_b44d_sig/{key}/s/recon.mat')
    d = load(f'{EXT}/AIkill_Dynamic_b44d/{key}/d/recon.mat')
    print(f'== {key}')
    print('  s keys:', sorted(s))
    print('  d-only keys:', sorted(set(d) - set(s)), ' s-only:', sorted(set(s) - set(d)))
    if 'calcb_b' in s and 'calcb_b' in d:
        bs, bd = s['calcb_b'], d['calcb_b']
        diff = np.abs(bs - bd)
        print(f'  calcb_b: shape {bs.shape} max|diff| {diff.max():.3e} (|b| max {np.abs(bd).max():.3f}); '
              f'identical={np.array_equal(bs, bd)}')
    gs, gd = s['gas_phase'], d['gas_phase']
    print(f'  gas_phase shape {gs.shape} dtype {gs.dtype}; s min/max {gs.min():.3e}/{gs.max():.3e}; '
          f'd min/max {gd.min():.3e}/{gd.max():.3e}')
    # corner noise (zero-mean check), both trees
    for nm, g in (('s', gs), ('d', gd)):
        c = g[:, :8, :8, :8]
        print(f'  {nm} corner gas_phase: mean {c.mean():.2f} std {c.std():.2f}  (mean/std {c.mean() / c.std():.3f})')
    sig = gs[:, :8, :8, :8].std()
    # pooled scale: bin-mean image of s vs d (both weighted means of the same raw -> same pooled lung sum)
    thr = 10 * sig
    ps, pd_ = gs.mean(0), gd.mean(0)
    ms, md = ps > thr, pd_ > thr
    print(f'  pooled bin-mean image: s sum>{thr:.0f} {ps[ms].sum():.4e} ({ms.sum()} vox) | '
          f'd sum {pd_[md].sum():.4e} ({md.sum()} vox) | ratio s/d {ps[ms].sum() / pd_[md].sum():.4f}')
    # per-bin masked sums (per-bin mask, 10 sigma) and real/mag ratio
    gm = s.get('gas_phase_magnitude')
    rows = []
    for b in range(gs.shape[0]):
        m = gs[b] > thr
        tot = gs[b][m].sum()
        r = gs[b][m].sum() / gm[b][m].sum() if gm is not None else np.nan
        rows.append((b, int(m.sum()), tot, r))
    bins = None
    p = f'{BT}/{key}_bins.csv'
    if os.path.exists(p):
        import csv
        bins = {int(r['bin']): r for r in csv.DictReader(open(p)) if r['binning'] == 'SIGNAL'}
    print('  bin  nvox   sum_gas(real)   real/mag   | bins.csv k0_gas_wmean  t_after_EE  n_hard')
    for b, n, tot, r in rows:
        extra = ''
        if bins and b in bins:
            extra = f'| {float(bins[b]["k0_gas_wmean"]):.4f}  {float(bins[b]["t_after_EE_med_s"]):.3f}  {bins[b]["n_hard"]}'
        print(f'  {b:3d} {n:6d} {tot:14.4e} {r:9.4f}   {extra}')
    tots = np.array([r[2] for r in rows])
    if bins:
        k0 = np.array([float(bins[b]['k0_gas_wmean']) for b in range(len(rows))])
        rho = np.corrcoef(tots, k0)[0, 1]
        print(f'  corr(sum_gas per bin, k0_gas_wmean) = {rho:.3f}; argmin sum_gas = bin {tots.argmin()}, argmax = bin {tots.argmax()}')
    # dissolved: three global estimators per bin, lung mask from gas of the same bin
    if 'dissolved_phase_real' in s:
        dr, di = s['dissolved_phase_real'], s['dissolved_phase_imag']
        print(f'  rbc_tp_separated = {int(np.asarray(s["rbc_tp_separated"]).ravel()[0])}')
        print('  bin   sum(aRBC)    sum(aTP)    sum(a+a)    |sum z|(dphi 53)   sum|z|   RBC/TP')
        dphi = np.radians(53.0)
        for b in range(gs.shape[0]):
            m = gs[b] > thr
            a, t = dr[b][m].sum(), di[b][m].sum()
            z = dr[b][m] * np.exp(1j * dphi) + di[b][m]
            print(f'  {b:3d} {a:11.4e} {t:11.4e} {a + t:11.4e} {np.abs(z.sum()):11.4e} {np.abs(z).sum():11.4e} {a / t:7.3f}')


if __name__ == '__main__':
    for k in sys.argv[1:]:
        try:
            main(k)
        except Exception as e:
            print(f'{k}: FAILED {e!r}')
