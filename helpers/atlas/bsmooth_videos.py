#!/usr/bin/env python
"""Atlas-style bin videos, raw b vs smoothed b, for the sessions re-run with fork 19b7365.

Rows per orientation: 'raw b' (AIkill_Dynamic recon.mat + recon_resplit_xecs.mat, the atlas maps) and
'smooth b' (tyger_bsmooth_2026-09-25 output.mrd: gas real(F*b) + container split aRBC + i*aTP).
Same layout/normalisation as rbc_tm_atlas.py (C39): 10 coronal/sagittal/axial slices at shared fractions,
gas-picked window, 16 bins x 5 repeats at 5 fps.
Usage: bsmooth_videos.py <key> [...]
Output: outputs/bphase_fix_2026-09-25/videos/<ID>_{gp,rbc,tm}.mp4
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scipy.io as sio

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
import rbc_tm_atlas as A  # noqa: E402
from bphase_confirm import read_mrd  # noqa: E402

SRC = Path('/Volumes/HoomHamExt/AIkill_Dynamic')
NEW = Path('/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/tyger_bsmooth_2026-09-25')
OUTV = A.WS / 'outputs' / 'bphase_fix_2026-09-25' / 'videos'


def raw_session(key):
    s = A.load_session(key)
    s.date = 'raw b'
    return s


def new_session(key):
    gas, dis, b, meta = read_mrd(NEW / key / 'd' / 'output.mrd')
    vols = {'gp': gas.astype(np.float32), 'rbc': np.maximum(dis.real, 0).astype(np.float32),
            'tm': np.maximum(dis.imag, 0).astype(np.float32)}
    m = {'dphi': float(meta['rbc_tp_dphi_deg']), 'F_lump': float(meta.get('rbc_tp_F_lump', 'nan')),
         'insp_bin': -1, 'ratio_scalar': float(meta.get('rbc_tp_ratio_scalar', 'nan')), 'carrier': '', 'gate': ''}
    s = A.Session(key, vols, m)
    s.date = 'smooth b'
    return s


def main(keys):
    A.VID = OUTV
    OUTV.mkdir(parents=True, exist_ok=True)
    for key in keys:
        sid = key[11:]
        idg = A.IDGroup([raw_session(key), new_session(key)])
        for kind in A.KINDS:
            A.render_video(f'{sid}', idg, kind)
        print('[video]', key, '->', OUTV, flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
