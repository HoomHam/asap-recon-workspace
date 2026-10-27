#!/usr/bin/env python
"""Minimal output.mrd -> recon.mat for the production (b44) tree. No figures, no gifs (post_process.py does those).

Writes: gas_phase (real(F*b)), gas_phase_magnitude (|F*b|, nch=1), dissolved_phase_real/imag (= aRBC / aTP when
rbc_tp_separated == 1, else the unsplit complex image), dissolved_phase_magnitude, calcb_b (complex, nch=1),
rbc_tp_* meta as scalars/strings, navigator arrays if present.
Usage: mrd_to_mat.py <session_dir_with_d> [...]   (e.g. /Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_b44/2024-09-10_045VS)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scipy.io as sio
import mrd

NAV_KEYS = ('nav_coronal', 'nav_diaphragm_z', 'nav_time', 'nav_volume', 'nav_ilvtime', 'nav_volmeastime',
            'ilvbin', 'ilvtime', 'ee_times', 'bin_weff', 'bin_mode', 'mix_tau_I', 'mix_n_insp', 'mix_tB_med', 'mix_T_med',
            'pneumo_clock_mapped', 'pneumo_clock_offset_s', 'pneumo_clock_T_cal_s', 'pneumo_clock_ilv_true_ms', 'pneumo_clock_chk_a', 'pneumo_clock_chk_b', 'pneumo_clock_chk_n_matched', 'pneumo_clock_chk_n_ee', 'pneumo_clock_chk_resid_sd_s', 'pneumo_clock_chk_ok', 'pneumo_time_aligned', 'pneumo_vol_aligned')


def convert(sess_dir):
    # <session>/d by default; a path that itself holds output.mrd (e.g. the SIGNAL tree's <session>/s) is used as is
    d = Path(sess_dir) if (Path(sess_dir) / 'output.mrd').is_file() else Path(sess_dir) / 'd'
    out = {}
    # keep keys an earlier recon.mat already carries (post_process.py: fid_signal, pneumo_*, diaphragm_pos)
    if (d / 'recon.mat').is_file():
        out.update({k: v for k, v in sio.loadmat(str(d / 'recon.mat')).items() if not k.startswith('__')})
    with mrd.BinaryMrdReader(str(d / 'output.mrd')) as r:
        r.read_header()
        for it in r.read_data():
            if isinstance(it, mrd.StreamItem.NdArrayFloat):
                m = it.value.meta
                if m.get('gas_phase_image'):
                    out['gas_phase'] = np.asarray(it.value.data, dtype=np.float32)
                    out['deapod'] = int(m['deapod'][0].value) if m.get('deapod') else 0   # fork ee3c91f+: rolloff divided
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
    out['made_by'] = 'workspace/pipeline/mrd_to_mat.py (fork calcb_phase_sigma 4.4; deapod flag from output meta)'
    sio.savemat(str(d / 'recon.mat'), out, do_compression=False)
    return {k: (v.shape if hasattr(v, 'shape') else v) for k, v in out.items() if k.startswith(('gas_phase', 'dissolved_phase_real', 'rbc_tp_separated', 'rbc_tp_dphi', 'calcb'))}


if __name__ == '__main__':
    for s in sys.argv[1:]:
        try:
            print('[mat]', Path(s).name, convert(s), flush=True)
        except Exception as e:
            print('[mat]', Path(s).name, 'FAILED', repr(e), flush=True)
