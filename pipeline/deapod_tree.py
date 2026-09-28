#!/usr/bin/env python
"""Offline deapodised copy of the b44 production tree: every image in recon.mat divided by the gridding kernel's
rolloff A(r) -- exactly what fork ee3c91f does inside dyn_recon (after the RBC/TP split, not in the navigator), so
recon.mat here == what a Tyger rerun with ee3c91f would give (checked on 025JC / 011CN, run_deapod_check_2026-09-28.sh).

Divided: gas_phase, gas_phase_magnitude, dissolved_phase_real, dissolved_phase_imag, dissolved_phase_magnitude.
Unchanged: calcb_b, nav_*, rbc_tp_* (split phase / mask are computed before the division in the fork too).
Added: deapod = 1, deapod_source, made_by.   MS / IS read from each session's input.mrd header.

usage: deapod_tree.py [KEY ...]          (default: every session in the b44 tree; existing outputs skipped)
"""
import json
import os
import sys
import time

import numpy as np
import scipy.io as sio
import mrd

SRC = '/Volumes/HoomHamExt/AIkill_Dynamic_b44'
RAW = '/Volumes/HoomHamExt/AIkill_Dynamic'
DST = '/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_b44d'
IMAGES = ('gas_phase', 'gas_phase_magnitude', 'dissolved_phase_real', 'dissolved_phase_imag', 'dissolved_phase_magnitude')
KDIST0SQ = 0.2   # recon.py cudarecon kernel; results.KDIST0SQ in the fork


def rolloff(MS, IS):
    """results.rolloff (fork ee3c91f): A on the IS crop of the MS grid, Steve's crop indices."""
    c = int(MS / 2 + .1)
    ll = int(c - IS / 2 + .1)
    r = np.arange(ll, ll + IS) - c
    r2 = r[:, None, None] ** 2 + r[None, :, None] ** 2 + r[None, None, :] ** 2
    return np.exp(-2 * np.pi ** 2 * (KDIST0SQ / 2.0) * r2 / MS ** 2)


def ms_is(key):
    r = mrd.BinaryMrdReader(os.path.join(RAW, key, 'd', 'input.mrd'))
    h = r.read_header()
    ul = {p.name: p.value for p in h.user_parameters.user_parameter_long}
    return int(ul.get('MS', 240)), int(ul.get('IS', 100))


def one(key):
    src = os.path.join(SRC, key, 'd', 'recon.mat')
    dst_dir = os.path.join(DST, key, 'd')
    dst = os.path.join(dst_dir, 'recon.mat')
    if os.path.exists(dst):
        return 'exists, skip'
    MS, IS = ms_is(key)
    A = rolloff(MS, IS).astype(np.float32)[None]
    m = {k: v for k, v in sio.loadmat(src).items() if not k.startswith('__')}
    if int(np.ravel(m.get('deapod', [0]))[0]) == 1:
        raise RuntimeError('source already deapodised')
    done = []
    for k in IMAGES:
        if k in m:
            assert m[k].shape[-3:] == (IS, IS, IS), (k, m[k].shape, IS)
            m[k] = (m[k] / A).astype(np.float32)
            done.append(k)
    m['deapod'] = 1
    m['deapod_source'] = f'{src} / A(r), kdist0sq {KDIST0SQ}, MS {MS}, IS {IS} (== fork ee3c91f dyn_recon)'
    m['made_by'] = 'workspace/pipeline/deapod_tree.py 2026-09-28'
    os.makedirs(dst_dir, exist_ok=True)
    sio.savemat(dst + '.tmp.mat', m, do_compression=False)
    os.replace(dst + '.tmp.mat', dst)
    return f'MS {MS} IS {IS} min A {A.min():.3f}; divided {len(done)}: {",".join(done)}'


def main():
    keys = sys.argv[1:] or sorted(k for k in os.listdir(SRC) if os.path.exists(os.path.join(SRC, k, 'd', 'recon.mat')))
    os.makedirs(DST, exist_ok=True)
    log = []
    for k in keys:
        t0 = time.time()
        try:
            msg = one(k)
        except Exception as e:
            msg = 'FAILED ' + repr(e)
        log.append({'key': k, 'result': msg})
        print(f'[deapod] {k}: {msg} ({time.time() - t0:.0f} s)', flush=True)
    prov = os.path.join(DST, 'PROVENANCE.json')
    old = json.load(open(prov)) if os.path.exists(prov) else {'runs': []}
    old['what'] = ('b44 production recon.mat with every image divided by the gridding-kernel rolloff; identical to a '
                   'Tyger rerun with fork ee3c91f (deapod). Source tree ' + SRC + '. ASAP facts F99 / card C55.')
    old['runs'].append({'at': time.strftime('%Y-%m-%d %H:%M:%S'), 'sessions': log})
    json.dump(old, open(prov, 'w'), indent=1)


if __name__ == '__main__':
    main()
