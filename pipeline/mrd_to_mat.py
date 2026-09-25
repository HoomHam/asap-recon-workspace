#!/usr/bin/env python
"""Minimal output.mrd -> recon.mat for the production (b44) tree. No figures, no gifs (post_process.py does those).

Writes: gas_phase (real(F*b)), gas_phase_magnitude (|F*b|, nch=1), dissolved_phase_real/imag (= aRBC / aTP when
rbc_tp_separated == 1, else the unsplit complex image), dissolved_phase_magnitude, calcb_b (complex, nch=1),
rbc_tp_* meta as scalars/strings, navigator arrays if present.
Usage: mrd_to_mat.py <session_dir_with_d> [...]   (e.g. /Volumes/HoomHamExt/AIkill_Dynamic_b44/2024-09-10_045VS)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scipy.io as sio
import mrd

NAV_KEYS = ('nav_coronal', 'nav_diaphragm_z', 'nav_time', 'nav_volume', 'nav_ilvtime', 'nav_volmeastime')


def convert(sess_dir):
    d = Path(sess_dir) / 'd'
    out = {}
    with mrd.BinaryMrdReader(str(d / 'output.mrd')) as r:
        r.read_header()
        for it in r.read_data():
            if isinstance(it, mrd.StreamItem.NdArrayFloat):
                m = it.value.meta
                if m.get('gas_phase_image'):
                    out['gas_phase'] = np.asarray(it.value.data, dtype=np.float32)
                elif m.get('gas_phase_magnitude'):
                    out['gas_phase_magnitude'] = np.asarray(it.value.data, dtype=np.float32)
                else:
                    for k in NAV_KEYS:
                        if m.get(k):
                            out[k] = np.asarray(it.value.data, dtype=np.float32)
            elif isinstance(it, mrd.StreamItem.NdArrayComplexFloat):
                m = it.value.meta
                if m.get('dissolved_phase_image'):
                    z = np.asarray(it.value.data, dtype=np.complex64)
                    out['dissolved_phase_real'] = z.real.astype(np.float32)
                    out['dissolved_phase_imag'] = z.imag.astype(np.float32)
                    out['dissolved_phase_magnitude'] = np.abs(z).astype(np.float32)
                    for k, v in m.items():
                        if k.startswith('rbc_tp'):
                            out[k] = str(v[0].value)
                elif m.get('calcb_b'):
                    out['calcb_b'] = np.asarray(it.value.data, dtype=np.complex64)
                elif m.get('gas_phase_complex'):
                    pass  # 128 MB, not needed offline
    out['made_by'] = 'workspace/pipeline/mrd_to_mat.py 2026-09-25 (fork calcb_phase_sigma 4.4)'
    sio.savemat(str(d / 'recon.mat'), out, do_compression=False)
    return {k: (v.shape if hasattr(v, 'shape') else v) for k, v in out.items() if k.startswith(('gas_phase', 'dissolved_phase_real', 'rbc_tp_separated', 'rbc_tp_dphi', 'calcb'))}


if __name__ == '__main__':
    for s in sys.argv[1:]:
        try:
            print('[mat]', Path(s).name, convert(s), flush=True)
        except Exception as e:
            print('[mat]', Path(s).name, 'FAILED', repr(e), flush=True)
