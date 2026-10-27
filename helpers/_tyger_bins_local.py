"""Scratch (2026-10-26, PCA Step 1): reproduce Tyger's SIGNAL bin membership OFFLINE, exactly, from <s>/input.mrd.

Runs the root (fork, diaphragm-recon) raw.load_from_arr on the stored input.mrd — the same code the 262b16f image runs
before gridding — on the CPU (no GPU needed for loading / binning). Exports what the fork does not: per-interleave
bin fraction (g_raw.ilvbin[SIGNAL]), interleave time, excluded flag (auto low-signal frame ranges, results.dyn_recon
mask rule), Tyger's k0 signal (volmeasvol, min-max normalised) and the EE times.

Outputs (workspace/outputs/bin_timing_2026-10-26/):
  <key>_ilvbin_tyger.csv  ilv, t_s, ilvbin_frac, excluded, k0_tyger_norm
  <key>_ee_tyger.csv      EE times (s)
  summary line: EE count, exclude ranges (s), n excluded interleaves, n binned

usage (conda base python, from the repo root so gtypes/raw/specfit import):
  cd <root> && python3 workspace/helpers/_tyger_bins_local.py <key> ...
"""
import os
import re
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, ROOT)
import mrd                                            # noqa: E402
from gtypes import gvar, imgtype, bintype             # noqa: E402
from raw import traj, raw                             # noqa: E402

EXT = '/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon'
OUT = os.path.join(ROOT, 'workspace', 'outputs', 'bin_timing_2026-10-26')


def load(input_mrd):
    """tyger_recon.reconstruct_from_mrd lines 183-261, up to (not including) calcb."""
    g, g_traj, g_raw = gvar(), traj(), raw()
    ref_acq_arr = dyn_acq_arr = pneumo_arr = gp = dp = None
    pneumo_tmap = None
    with mrd.BinaryMrdReader(input_mrd) as reader:
        header = reader.read_header()
        for item in reader.read_data():
            if isinstance(item, (mrd.StreamItem.NdArrayDouble, mrd.StreamItem.NdArrayFloat)):
                if item.value.meta.get('gas_phase_trajectory'):
                    gp = item.value.data
                elif item.value.meta.get('dissolved_phase_trajectory'):
                    dp = item.value.data
                elif item.value.meta.get('pneumotach'):
                    pneumo_arr = item.value.data
                elif item.value.meta.get('pneumo_gas_true_s'):
                    pneumo_tmap = item.value.data
            elif isinstance(item, mrd.StreamItem.NdArrayComplexFloat):
                if item.value.meta.get('reference_acquisition'):
                    ref_acq_arr = item.value.data
                elif item.value.meta.get('dynamic_acquisition'):
                    dyn_acq_arr = item.value.data
    user_long, user_double, user_string = {}, {}, {}
    if header.user_parameters is not None:
        user_long = {p.name: p.value for p in header.user_parameters.user_parameter_long}
        user_double = {p.name: p.value for p in header.user_parameters.user_parameter_double}
        user_string = {p.name: p.value for p in header.user_parameters.user_parameter_string}
    g.MS = int(user_long.get('MS', g.MS))
    g.IS = int(user_long.get('IS', g.IS))
    g.nbins = int(user_long.get('nbins', g.nbins))
    killpts = int(user_long.get('killpts', 2))
    meta = {'TR': user_double.get('TR', 0.0), 'TE': user_double.get('TE', 0.0),
            'DPoff': user_double.get('DPoff', 0.0), 'dtdyn': user_double.get('dtdyn', 0.0),
            'dtspec': user_double.get('dtspec', 0.0), 'numspec': int(user_long.get('numspec', 0)),
            'sigminN': int(user_long.get('sigminN', 50))}
    if int(user_long.get('pneumoclock', 0)) and pneumo_tmap is not None:
        meta['pneumo_tmap'] = np.asarray(pneumo_tmap, dtype='float64')
        meta['pneumo_offset_s'] = float(user_double.get('pneumo_offset_s', 0.0))
    g_traj.killpts = killpts
    g_traj.load_traj_from_array(gp, dp, int(user_long.get('nusimg', 32)))
    dyn_acq_arr = dyn_acq_arr[:, killpts:, :]
    if ref_acq_arr is not None:
        ref_acq_arr = ref_acq_arr[:, killpts:, :]
    g_raw.load_from_arr(g_traj, ref_acq_arr, dyn_acq_arr, pneumo_arr, 'mrd_siemens', meta)
    exstr = user_string.get('excluderanges', '')
    if exstr:
        a = list(map(int, re.findall(r'\d+', exstr)))
        manual = [] if len(a) % 2 else [range(a[j], a[j + 1]) for j in range(0, len(a), 2) if a[j + 1] > a[j]]
        g_raw.excluderanges += manual
    return g, g_raw, user_string.get('binning', 'SIGNAL')


def main(key, sub=None, tree=None):
    # input.mrd source: env TREE/SUB, else the b44d_sig s/ copy, else the old tree d/ (same raw, same EE/exclusions:
    # load_from_arr computes the SIGNAL EE list and the auto exclusions whatever the header's binning string says)
    tree = tree or os.environ.get('TREE'); sub = sub or os.environ.get('SUB')
    cands = ([(tree, sub)] if tree and sub else []) + [('AIkill_Dynamic_b44d_sig', 's'), ('AIkill_Dynamic', 'd')]
    for tr, sb in cands:
        p = f'{EXT}/{tr}/{key}/{sb}/input.mrd'
        if os.path.exists(p):
            break
    else:
        raise FileNotFoundError(f'{key}: no input.mrd in {cands}')
    g, g_raw, binning = load(p)
    bt = bintype.SIGNAL
    bins = np.asarray(g_raw.ilvbin[bt], dtype=float)
    t = np.asarray(g_raw.ilvtime, dtype=float)
    # results.dyn_recon mask rule (results.py:292-296): interleave excluded if iilv*TR inside any exclude range
    excl = np.zeros(len(t), dtype=int)
    for iilv in range(g_raw.ntotalilvs):
        for el in g_raw.excluderanges:
            if iilv * g_raw.TR >= el.start and iilv * g_raw.TR <= el.stop:
                excl[iilv] = 1
    k0 = np.asarray(g_raw.volmeasvol[bt], dtype=float)
    ee = np.asarray(g_raw.volmeasEEtime[bt], dtype=float)
    os.makedirs(OUT, exist_ok=True)
    np.savetxt(f'{OUT}/{key}_ilvbin_tyger.csv', np.column_stack([np.arange(len(t)), t, bins, excl, k0]),
               delimiter=',', header='ilv,t_s,ilvbin_frac,excluded,k0_tyger_norm', comments='',
               fmt='%d,%.4f,%.6f,%d,%.6f')
    np.savetxt(f'{OUT}/{key}_ee_tyger.csv', ee, delimiter=',', header='t_EE_s', comments='', fmt='%.4f')
    print(f'{key}: binning={binning} nilv {len(t)} TR {g_raw.TR:.4f} EE {len(ee)} | exclude ranges (s) '
          f'{[(r.start, r.stop) for r in g_raw.excluderanges]} -> {excl.sum()} ilv excluded, '
          f'{int(np.sum(bins >= 0))} binned, {int(np.sum((bins >= 0) & (excl == 0)))} binned & kept')


if __name__ == '__main__':
    for k in sys.argv[1:]:
        try:
            main(k)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f'{k}: FAILED {e!r}')
