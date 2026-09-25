#!/usr/bin/env python3
"""
resplit_rbctp.py -- offline RBC/TP re-split of Steve's DIAPHRAGM recons with the corrected
basis angle (2steve notes 04 + 05, facts F59/F60). No Tyger re-run, root code untouched.

Why: results.py:324-325 sets the RBC/TP basis angle to 2*pi*df*TEeff, which is the phase
difference extrapolated to t = 0 (RF intercept, 9-32 deg), not the one at the image's k0
(t = TE + killpts*dt, measured 88-103 deg). The 2x2 split is linear and invertible, so the
saved aRBC/aTP (recon.mat dissolved_phase_real/imag) give back the complex dissolved image
per bin, which is then re-split with the measured angle and the note-04 stop criterion
(masked ratio, min |R - target| among phases with both masked sums positive).

Stages (run in this order; all idempotent):
  fit     [keys|--all] [-j N]  rerun Steve's spectral fit (raw.load_from_arr, CPU, unchanged
                               code) per session in parallel subprocesses -> fits/<key>.json
  fit-one <key>                one session, in-process (used by `fit`)
  split   [keys|--all] [--dry] re-split from recon.mat + fits/<key>.json ->
                               <Ext>/<key>/d/recon_resplit.mat (never touches recon.mat) and
                               rows/<key>.json (metrics)
  summary                      rows/*.json -> resplit_summary.csv
  fig     <key> [...]          before/after figure -> fig/<key>_before_after.png
  cohort                       cohort scatter (angles, corr, TP<0) -> fig/cohort_before_after.png

Inputs : /Volumes/HoomHamExt/AIkill_Dynamic/<key>/d/{input.mrd, recon.mat, tyger.log}
Outputs: workspace/outputs/resplit_2026-09-24/{fits,rows,logs,fig}/, resplit_summary.csv
         <Ext>/<key>/d/recon_resplit.mat   (aRBC, aTP float32 (nbins,z,y,x) + scalars)

Gate (session skipped for the .mat, still summarised): spectral fit missing; df outside
[DF_MIN, DF_MAX] Hz; fRBC further than FRBC_TOL from the DPoff-expected offset; |sin dphi_new|
< SIN_MIN. Defaults are loose (250-400 Hz, +-150 Hz); --strict uses the note-05 numbers
(300-380 Hz, +-100 Hz). Both verdicts are written to the summary.

Python: conda base (/opt/homebrew/Caskroom/miniforge/base/bin/python, has `mrd`).
"""
import argparse
import csv
import json
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]                        # 2026_ASAP_Recon (Steve's code, read-only)
OUT = HERE.parent / 'outputs' / 'resplit_2026-09-24'
DYN = Path('/Volumes/HoomHamExt/AIkill_Dynamic')
PY = sys.executable
GAMMA_HZ_PER_PPM = 11.777e6 * 1.494 * 1e-6    # 129Xe at the Avanto (F54): 17.6 Hz/ppm
RBC_PPM = 218.0                               # RBC resonance used for the expected fRBC
SIN_MIN = 0.3                                 # |sin dphi| gate (noise gain 1/|sin| > 3.3)
RATIO_MIN, RATIO_MAX = 0.02, 3.0              # spectral RBC:TP target must be a physical ratio
COND = ((0.7, 'good'), (SIN_MIN, 'marginal'), (0.0, 'poor'))   # |sin dphi_new| classes
STAB_MAX_DEG = 10.0                           # xecs: max k0-angle change when one more leading sample is dropped (F69)
GATES = dict(loose=dict(DF_MIN=250.0, DF_MAX=400.0, FRBC_TOL=150.0),
             strict=dict(DF_MIN=300.0, DF_MAX=380.0, FRBC_TOL=100.0))
FIT = 'steve'                                 # 'steve' (raw.py fit rerun) | 'xecs' (calspec prior fit csv)
SFX = ''                                      # output suffix, '' or '_xecs' (set by main)
XECS_CSV_DEFAULTS = [HERE.parent / 'notes' / 'calspec_package_2026-09-16' / 'resplit_inputs_2026-09-24.csv',
                     ROOT.parent / '2026_XeCS_Recon' / 'workspace' / 'outputs' / 'calspec' / 'resplit_inputs_2026-09-24.csv']
_XECS = None


def xecs_rows(path=None):
    """Per-session split inputs from the XeCS prior-constrained time-domain fit (F65–F67):
    dphi_k0_deg (M3 lumped membrane, RBC − mem, at TE + killpts*dt_img), ratio_lumped (a_RBC/|mem1+mem2|,
    the split's self-consistent target), ratio_scalar, F_lump, df_hz, valid/reason, ..."""
    global _XECS
    if _XECS is None:
        p = Path(path) if path else next((q for q in XECS_CSV_DEFAULTS if q.is_file()), None)
        if p is None or not p.is_file():
            raise SystemExit(f'--fit xecs: no csv at {path or XECS_CSV_DEFAULTS}')
        _XECS = {}
        for r in csv.DictReader(open(p)):
            k = r.get('session') or r.get('key')
            if r.get('primary', '1') in ('0', 'False', 'false') and k in _XECS:
                continue                         # keep the primary twix row per session
            _XECS[k] = r
        print(f'xecs fit table: {p} ({len(_XECS)} sessions)', file=sys.stderr)
    return _XECS


SOLVED_RE = re.compile(r'DP bin (\d+): RBC/TP phase solved at ph=([0-9.]+) rad \(R=([-0-9.eE+]+) '
                       r'crossed target RBCTPratio=([-0-9.eE+]+)\)')


def wrap_deg(d):
    return (d + 180.0) % 360.0 - 180.0


def all_split_keys():
    """Sessions whose Tyger run split at least one bin (the 86 of note 04 + merged reruns)."""
    keys = []
    for d in sorted(DYN.iterdir()):
        log = d / 'd' / 'tyger.log'
        if (d / 'd' / 'recon.mat').is_file() and log.is_file() and \
                'RBC/TP phase solved' in log.read_text(errors='replace'):
            keys.append(d.name)
    return keys


# ----------------------------------------------------------------------------- fit stage
def fit_one(key):
    """Rerun Steve's spectral fit on d/input.mrd exactly as tyger_recon.py feeds it
    (killpts trimmed, header meta in SI) and return the parameters the split needs."""
    import warnings
    warnings.filterwarnings('ignore')
    import mrd
    sys.path.insert(0, str(ROOT))
    from gtypes import gvar            # noqa: F401  (raw imports it)
    from raw import traj, raw
    p = DYN / key / 'd' / 'input.mrd'
    gt, gr = traj(), raw()
    ref = dyn = pne = gtr = dtr = None
    with mrd.BinaryMrdReader(str(p)) as r:
        h = r.read_header()
        for it in r.read_data():
            if isinstance(it, (mrd.StreamItem.NdArrayDouble, mrd.StreamItem.NdArrayFloat)):
                m = it.value.meta
                if m.get('gas_phase_trajectory'):
                    gtr = it.value.data
                elif m.get('dissolved_phase_trajectory'):
                    dtr = it.value.data
                elif m.get('pneumotach'):
                    pne = it.value.data
            elif isinstance(it, mrd.StreamItem.NdArrayComplexFloat):
                m = it.value.meta
                if m.get('reference_acquisition'):
                    ref = it.value.data
                elif m.get('dynamic_acquisition'):
                    dyn = it.value.data
    ul = {q.name: q.value for q in h.user_parameters.user_parameter_long}
    ud = {q.name: q.value for q in h.user_parameters.user_parameter_double}
    kp = int(ul.get('killpts', 2))
    gt.killpts = kp
    gt.load_traj_from_array(gtr, dtr, int(ul.get('nusimg', 32)))
    dyn = dyn[:, kp:, :]
    ref = ref[:, kp:, :] if ref is not None else None
    meta = {'TR': ud['TR'], 'TE': ud['TE'], 'DPoff': ud['DPoff'], 'dtdyn': ud['dtdyn'],
            'dtspec': ud['dtspec'], 'numspec': int(ul['numspec'])}
    t0 = time.time()
    gr.load_from_arr(gt, ref, dyn, pne, 'mrd_siemens', meta)
    d = dict(key=key, nch=int(gr.nch), ilvperTR=int(gr.ilvperTR), TE=float(gr.TE), TR=float(gr.TR),
             DPoff=float(ud['DPoff']), dtdyn=float(ud['dtdyn']), dtspec=float(ud['dtspec']),
             numspec=int(ul['numspec']), killpts=kp, spectBW=float(gt.spectBW), trajBW=float(gt.BW),
             fit_seconds=round(time.time() - t0, 1), fit_ok=bool(len(gr.fRBC)))
    if d['fit_ok']:
        d.update(fRBC=float(gr.fRBC[0]), fTP=float(gr.fTP[0]), ratio=float(gr.RBCTPratio[0]),
                 RBCphase=[float(v) for v in gr.RBCphase], TPphase=[float(v) for v in gr.TPphase],
                 TEphase=[float(v) for v in gr.TEphase], deltaphase=float(gr.deltaphase),
                 TEeff=float(gr.TEeff), TEeff2=float(gr.TEeff2))
    return d


def cmd_fit(keys, jobs):
    (OUT / 'fits').mkdir(parents=True, exist_ok=True)
    (OUT / 'logs').mkdir(parents=True, exist_ok=True)
    todo = [k for k in keys if not (OUT / 'fits' / f'{k}.json').is_file()]
    print(f'fit: {len(keys)} sessions, {len(todo)} to do, {jobs} parallel', flush=True)

    def run(k):
        t0 = time.time()
        with open(OUT / 'logs' / f'{k}.fit.log', 'w') as lg:
            r = subprocess.run([PY, __file__, 'fit-one', k], stdout=subprocess.PIPE, stderr=lg, text=True)
        ok = r.returncode == 0 and (OUT / 'fits' / f'{k}.json').is_file()
        print(f'  {k}: {"ok" if ok else "FAILED rc=%d" % r.returncode} ({time.time() - t0:.0f} s)', flush=True)
        return ok

    with ThreadPoolExecutor(max_workers=jobs) as ex:
        res = list(ex.map(run, todo))
    print(f'fit done: {sum(res)}/{len(todo)} ok')


def cmd_fit_one(key):
    d = fit_one(key)
    (OUT / 'fits').mkdir(parents=True, exist_ok=True)
    json.dump(d, open(OUT / 'fits' / f'{key}.json', 'w'), indent=1)
    print(json.dumps({k: v for k, v in d.items() if not isinstance(v, list)}))


# --------------------------------------------------------------------------- split stage
def basis_angles(f):
    """(dphi_old_deg, dphi_new_deg) from a fit dict. Old = what results.py:324-325 used.
    New = fitted phase difference at the first spectral sample (t = TE + killpts*dtspec),
    moved to the image k0 (t = TE + killpts*dtdyn), differenced as complex numbers."""
    df = f['fRBC'] - f['fTP']
    old = np.degrees(2 * np.pi * df * f['TEeff'])
    meas = np.angle(np.exp(1j * (f['RBCphase'][0] - f['TPphase'][0])))
    corr = 2 * np.pi * df * f['killpts'] * (f['dtdyn'] - f['dtspec'])
    new = np.degrees(meas + corr)
    return float(wrap_deg(old)), float(wrap_deg(new)), float(np.degrees(corr))


def gate(f, dphi_new, mode):
    g = GATES[mode]
    reasons = []
    if not f.get('fit_ok'):
        return ['no_fit']
    df = f['fRBC'] - f['fTP']
    if not (g['DF_MIN'] <= df <= g['DF_MAX']):
        reasons.append(f'df={df:.0f}Hz')
    if not (RATIO_MIN <= f['ratio'] <= RATIO_MAX):
        reasons.append(f'ratio={f["ratio"]:.3f}')
    f_exp = (RBC_PPM - f['DPoff']) * GAMMA_HZ_PER_PPM
    if abs(f['fRBC'] - f_exp) > g['FRBC_TOL']:
        reasons.append(f'fRBC={f["fRBC"]:.0f}Hz(exp{f_exp:.0f})')
    if abs(np.sin(np.radians(dphi_new))) < SIN_MIN:
        reasons.append(f'|sin dphi|={abs(np.sin(np.radians(dphi_new))):.2f}')
    return reasons


def parse_solved(key):
    """bin -> (ph, R_logged, target) for bins the Tyger run actually split."""
    txt = (DYN / key / 'd' / 'tyger.log').read_text(errors='replace')
    return {int(b): (float(ph), float(R), float(t)) for b, ph, R, t in SOLVED_RE.findall(txt)}


def sweep_phase(S, dphi_rad, target, step=0.001):
    """Global phase from the masked complex sum S only (the split is linear, so the masked
    sums of aRBC/aTP at any ph follow from S*exp(i ph)). Basis: phiTP = 0, phiRBC = dphi.
    Returns (ph, R, sumRBC, sumTP, n_valid) at the ph minimising |R - target| among phases
    where both masked sums are positive; ph = nan if no such phase exists."""
    ph = np.arange(0.0, 2 * np.pi, step)
    s = S * np.exp(1j * ph)
    sR = s.imag / np.sin(dphi_rad)
    sT = s.real - sR * np.cos(dphi_rad)
    ok = (sR > 0) & (sT > 0)
    if not ok.any():
        return np.nan, np.nan, np.nan, np.nan, 0
    R = np.full(ph.shape, np.inf)
    R[ok] = sR[ok] / sT[ok]
    i = int(np.argmin(np.abs(R - target)))
    return float(ph[i]), float(R[i]), float(sR[i]), float(sT[i]), int(ok.sum())


def split_maps(z, dphi_rad):
    """z = aRBC*exp(i dphi) + aTP  ->  (aRBC, aTP)."""
    aRBC = z.imag / np.sin(dphi_rad)
    aTP = z.real - aRBC * np.cos(dphi_rad)
    return aRBC, aTP


def resplit_session(key, mode, write=True):
    import scipy.io as sio
    from scipy import ndimage
    sys.path.insert(0, str(HERE))
    import snr_calc as S

    fj = OUT / 'fits' / f'{key}.json'
    if not fj.is_file():
        return dict(key=key, status='no_fit_json')
    f = json.load(open(fj))
    row = dict(key=key, fit=FIT, nch=f['nch'], DPoff=f['DPoff'], TE_ms=f['TE'] * 1e3, fit_ok=f['fit_ok'])
    if not f['fit_ok']:
        # Steve's fit is still needed to invert the OLD split (its basis angles produced the saved aRBC/aTP)
        row.update(status='gated', gate_loose='no_fit', gate_strict='no_fit')
        return _save_row(row)
    df = f['fRBC'] - f['fTP']
    dphi_old, dphi_new, kp_corr = basis_angles(f)
    # alternative estimate for comparison: Steve's own line (intercept + trim terms) with TE added,
    # i.e. 2*pi*df*(TEeff + TE); differs from dphi_new by the first sample's residual off the line
    dphi_line = float(wrap_deg(np.degrees(2 * np.pi * df * (f['TEeff'] + f['TE']))))
    target = f['ratio']
    row.update(fRBC=f['fRBC'], fTP=f['fTP'], df=df, TEeff_us=f['TEeff'] * 1e6,
               dphi_old=dphi_old, dphi_steve_k0=dphi_new, dphi_line=dphi_line,
               line_minus_new=float(wrap_deg(dphi_line - dphi_new)), killpts_corr_deg=kp_corr,
               gain_old=1 / abs(np.sin(np.radians(dphi_old))), ratio_steve=target,
               fit_consistent=abs(wrap_deg(dphi_line - dphi_new)) <= 20.0)
    gl, gs = gate(f, dphi_new, 'loose'), gate(f, dphi_new, 'strict')
    row.update(gate_loose=';'.join(gl) or 'pass', gate_strict=';'.join(gs) or 'pass')
    if FIT == 'xecs':
        # angle + target from the XeCS prior-constrained fit; Steve's numbers stay as reference columns
        x = xecs_rows().get(key) or xecs_rows().get(key.replace('_merged', ''))   # merged inputs share the cal block
        if x is None:
            row.update(status='no_xecs_row', xecs_reason='session has no usable cal block in the XeCS table')
            return _save_row(row)
        fx = lambda k: float(x[k]) if x.get(k) not in (None, '', 'nan', 'NaN') else np.nan
        dphi_new, target = fx('dphi_k0_deg'), fx('ratio_lumped')
        row.update(dphi_new=dphi_new, target=target, df_xecs=fx('df_hz'), ratio_scalar=fx('ratio_scalar'),
                   F_lump=fx('F_lump'), ratio_m2=fx('ratio_m2'), dphi_k0_m2=fx('dphi_k0_m2_deg'),
                   dphi_rep_sd=fx('dphi_rep_sd_deg'), snr_diss=fx('snr_diss'), carrier=x.get('carrier', ''),
                   carrier_ppm=fx('carrier_ppm_above_gas'), model_used=x.get('model_used', ''),
                   xecs_valid=x.get('valid', ''), xecs_reason=x.get('reason', ''),
                   steve_minus_xecs=float(wrap_deg(row['dphi_steve_k0'] - dphi_new)) if np.isfinite(dphi_new) else np.nan)
        row.update(stab_dphi=fx('stab_dphi_deg'), stab_F=fx('stab_F_ratio'), xecs_stable=x.get('stable', ''))
        why = []
        if str(x.get('valid', '')).lower() in ('0', 'false', 'no'):
            why.append(f'xecs invalid: {x.get("reason", "")}')
        if np.isfinite(row['stab_dphi']) and row['stab_dphi'] > STAB_MAX_DEG:
            why.append(f'unstable: dphi moves {row["stab_dphi"]:.1f} deg when one more sample is dropped')
        if not np.isfinite(dphi_new) or not np.isfinite(target):
            why.append('xecs dphi/ratio missing')
        elif abs(np.sin(np.radians(dphi_new))) < SIN_MIN:
            why.append(f'|sin dphi|={abs(np.sin(np.radians(dphi_new))):.2f}')
        if np.isfinite(target) and not (RATIO_MIN <= target <= RATIO_MAX * 2):
            why.append(f'ratio={target:.3f}')
        row['gate_xecs'] = ';'.join(why) or 'pass'
        if why:
            row['status'] = 'gated'
            return _save_row(row)
    else:
        row.update(dphi_new=dphi_new, target=target)
        if (mode == 'strict' and gs) or (mode == 'loose' and gl):
            row['status'] = 'gated'
            return _save_row(row)
    sin_new = abs(np.sin(np.radians(dphi_new)))
    row.update(gain_new=1 / sin_new, cond_new=next(name for thr, name in COND if sin_new >= thr))

    m = sio.loadmat(DYN / key / 'd' / 'recon.mat',
                    variable_names=['gas_phase', 'dissolved_phase_real', 'dissolved_phase_imag'])
    gas = m['gas_phase'].astype(np.float64)
    re_, im_ = m['dissolved_phase_real'].astype(np.float64), m['dissolved_phase_imag'].astype(np.float64)
    nb = gas.shape[0]
    solved = parse_solved(key)
    # masks and background exactly as snr_calc (note 04 numbers are comparable)
    corners = S.corner_box(gas.shape[1:])
    gmean = gas.mean(0)
    bg = corners & ~ndimage.binary_dilation(S.lung_mask(gmean, gmean[corners].std()), iterations=S.EXCL_DILATE)
    masks = [S.lung_mask(gas[b], gas[b][bg].std()) for b in range(nb)]
    vols = np.array([mk.sum() for mk in masks])
    b_in = int(np.argmax(vols))

    phR_old, phT_old = 2 * np.pi * f['fRBC'] * f['TEeff'], 2 * np.pi * f['fTP'] * f['TEeff']
    dnew = np.radians(dphi_new)
    f = dict(f, ratio=target)                    # sweep target: Steve's area ratio or the XeCS lumped ratio
    aRBC = np.zeros_like(re_, dtype=np.float32)
    aTP = np.zeros_like(re_, dtype=np.float32)
    per = dict(ph=np.full(nb, np.nan), R=np.full(nb, np.nan), nvalid=np.zeros(nb, int),
               was_split=np.zeros(nb, bool), corr_old=np.full(nb, np.nan), corr_new=np.full(nb, np.nan),
               negTP_old=np.full(nb, np.nan), negTP_new=np.full(nb, np.nan),
               negRBC_old=np.full(nb, np.nan), negRBC_new=np.full(nb, np.nan),
               lungRBC_old=np.full(nb, np.nan), lungTP_old=np.full(nb, np.nan),
               lungRBC_new=np.full(nb, np.nan), lungTP_new=np.full(nb, np.nan),
               ratio_lung_new=np.full(nb, np.nan), lung_vox=vols)
    for b in range(nb):
        mk = masks[b]
        if b in solved:
            z = re_[b] * np.exp(1j * phR_old) + im_[b] * np.exp(1j * phT_old)   # rotated complex image
            per['was_split'][b] = True
            per['corr_old'][b] = np.corrcoef(re_[b][bg], im_[b][bg])[0, 1]
            if mk.any():
                per['negTP_old'][b] = 100 * (im_[b][mk] < 0).mean()
                per['negRBC_old'][b] = 100 * (re_[b][mk] < 0).mean()
                per['lungRBC_old'][b] = re_[b][mk].mean()
                per['lungTP_old'][b] = im_[b][mk].mean()
        else:
            z = re_[b] + 1j * im_[b]                                            # stored unsplit
        if not mk.any():
            continue
        ph, R, sR, sT, nv = sweep_phase(z[mk].sum(), dnew, f['ratio'])
        per['ph'][b], per['R'][b], per['nvalid'][b] = ph, R, nv
        if np.isnan(ph):
            continue
        r_, t_ = split_maps(z * np.exp(1j * ph), dnew)
        aRBC[b], aTP[b] = r_.astype(np.float32), t_.astype(np.float32)
        per['corr_new'][b] = np.corrcoef(r_[bg], t_[bg])[0, 1]
        per['negTP_new'][b] = 100 * (t_[mk] < 0).mean()
        per['negRBC_new'][b] = 100 * (r_[mk] < 0).mean()
        per['lungRBC_new'][b] = r_[mk].mean()
        per['lungTP_new'][b] = t_[mk].mean()
        per['ratio_lung_new'][b] = r_[mk].mean() / t_[mk].mean() if t_[mk].mean() else np.nan

    row.update(status='resplit', nbins=nb, bins_split_old=int(per['was_split'].sum()),
               bins_resplit=int(np.isfinite(per['ph']).sum()), insp_bin=b_in,
               corr_old_insp=per['corr_old'][b_in], corr_new_insp=per['corr_new'][b_in],
               corr_old_med=np.nanmedian(per['corr_old']), corr_new_med=np.nanmedian(per['corr_new']),
               negTP_old_insp=per['negTP_old'][b_in], negTP_new_insp=per['negTP_new'][b_in],
               negTP_old_med=np.nanmedian(per['negTP_old']), negTP_new_med=np.nanmedian(per['negTP_new']),
               negRBC_old_insp=per['negRBC_old'][b_in], negRBC_new_insp=per['negRBC_new'][b_in],
               R_new_insp=per['R'][b_in], R_new_med=np.nanmedian(per['R']),
               ratio_lung_new_insp=per['ratio_lung_new'][b_in],
               ph_new_insp=per['ph'][b_in], ph_old_insp=solved.get(b_in, (np.nan,) * 3)[0],
               lungRBC_old_insp=per['lungRBC_old'][b_in], lungTP_old_insp=per['lungTP_old'][b_in],
               lungRBC_new_insp=per['lungRBC_new'][b_in], lungTP_new_insp=per['lungTP_new'][b_in])
    if write:
        out = DYN / key / 'd' / f'recon_resplit{SFX}.mat'
        extra = {}
        if FIT == 'xecs':
            extra = dict(ratio_scalar=row['ratio_scalar'], F_lump=row['F_lump'], dphi_steve_k0_deg=row['dphi_steve_k0'],
                         ratio_steve=row['ratio_steve'], model_used=row['model_used'], snr_diss=row['snr_diss'],
                         carrier=row['carrier'], dphi_rep_sd_deg=row['dphi_rep_sd'],
                         scale_note='aRBC/aTP maps are on the LUMPED scale (TP = |mem1+mem2|); divide the RBC/TP '
                                    'ratio by F_lump for the literature scalar scale a_RBC/(a1+a2)')
        sio.savemat(out, dict(aRBC=aRBC, aTP=aTP, ph_new=per['ph'], R_new=per['R'], nvalid=per['nvalid'],
                              was_split_old=per['was_split'], lung_vox=vols, insp_bin=b_in,
                              dphi_old_deg=dphi_old, dphi_new_deg=dphi_new, gain_new=row['gain_new'],
                              target=target, fRBC=f['fRBC'], fTP=f['fTP'], df=df, TEeff=f['TEeff'],
                              TE=f['TE'], DPoff=f['DPoff'], corr_bg_old=per['corr_old'], corr_bg_new=per['corr_new'],
                              negTP_lung_old_pct=per['negTP_old'], negTP_lung_new_pct=per['negTP_new'],
                              gate_loose=row['gate_loose'], gate_strict=row['gate_strict'], fit_source=FIT,
                              basis='z = aRBC*exp(1j*dphi_new) + aTP after global rotation ph_new; '
                                    'source = recon.mat dissolved_phase_real/imag inverted with dphi_old',
                              made_by='workspace/helpers/resplit_rbctp.py 2026-09-24', **extra))
        row['mat'] = str(out)
    return _save_row(row)


def _save_row(row):
    (OUT / f'rows{SFX}').mkdir(parents=True, exist_ok=True)
    json.dump({k: (None if isinstance(v, float) and np.isnan(v) else v) for k, v in row.items()},
              open(OUT / f'rows{SFX}' / f'{row["key"]}.json', 'w'), indent=1, default=float)
    return row


def cmd_split(keys, mode, write):
    for k in keys:
        t0 = time.time()
        r = resplit_session(k, mode, write)
        if r.get('status') == 'resplit':
            print(f"  {k}: dphi {r['dphi_old']:6.1f} -> {r['dphi_new']:6.1f} deg  corr(bg) {r['corr_old_insp']:+.2f} -> "
                  f"{r['corr_new_insp']:+.2f}  TP<0 lung {r['negTP_old_insp']:4.0f}% -> {r['negTP_new_insp']:4.0f}%  "
                  f"R {r['R_new_insp']:.3f} (target {r['target']:.3f})  {time.time() - t0:.0f} s", flush=True)
        else:
            print(f"  {k}: {r.get('status')} [{r.get('gate_xecs', r.get('gate_loose', ''))} {r.get('xecs_reason', '')}]", flush=True)


def cmd_summary():
    rows = [json.load(open(p)) for p in sorted((OUT / f'rows{SFX}').glob('*.json'))]
    fields = ['key', 'fit', 'status', 'gate_loose', 'gate_strict', 'gate_xecs', 'xecs_valid', 'xecs_reason', 'nch', 'DPoff',
              'TE_ms', 'fRBC', 'fTP', 'df', 'df_xecs', 'target', 'ratio_steve', 'ratio_scalar', 'F_lump', 'ratio_m2',
              'TEeff_us', 'dphi_old', 'dphi_new', 'dphi_steve_k0', 'steve_minus_xecs', 'dphi_k0_m2', 'dphi_rep_sd',
              'stab_dphi', 'stab_F', 'xecs_stable',
              'dphi_line', 'line_minus_new', 'fit_consistent', 'killpts_corr_deg', 'model_used', 'snr_diss', 'carrier',
              'carrier_ppm', 'gain_old', 'gain_new', 'cond_new', 'nbins',
              'bins_split_old', 'bins_resplit', 'insp_bin', 'ph_old_insp', 'ph_new_insp', 'R_new_insp', 'R_new_med',
              'ratio_lung_new_insp', 'corr_old_insp', 'corr_new_insp', 'corr_old_med', 'corr_new_med',
              'negTP_old_insp', 'negTP_new_insp', 'negTP_old_med', 'negTP_new_med', 'negRBC_old_insp',
              'negRBC_new_insp', 'lungRBC_old_insp', 'lungTP_old_insp', 'lungRBC_new_insp', 'lungTP_new_insp', 'mat']
    p = OUT / f'resplit_summary{SFX}.csv'
    with open(p, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow({k: (f'{v:.4g}' if isinstance(v, float) else v) for k, v in r.items()})
    n = len(rows)
    ok = [r for r in rows if r.get('status') == 'resplit']
    print(f'summary -> {p}: {n} sessions, {len(ok)} resplit, {n - len(ok)} gated/missing')
    if ok:
        co = np.array([r['corr_old_insp'] for r in ok], float)
        cn = np.array([r['corr_new_insp'] for r in ok], float)
        to = np.array([r['negTP_old_insp'] for r in ok], float)
        tn = np.array([r['negTP_new_insp'] for r in ok], float)
        dn = np.array([r['dphi_new'] for r in ok], float)
        dl = np.array([r['dphi_line'] for r in ok], float)
        print(f'  dphi_new: median {np.nanmedian(dn):.1f} deg, range {np.nanmin(dn):.1f}..{np.nanmax(dn):.1f}; '
              f'dphi_line (2pi*df*(TEeff+TE)): median {np.nanmedian(dl):.1f}, median |new-line| {np.nanmedian(np.abs(wrap_deg(dn - dl))):.1f}')
        print(f'  corr(bg) insp bin: median {np.nanmedian(co):+.2f} -> {np.nanmedian(cn):+.2f}')
        print(f'  lung TP<0 % insp bin: median {np.nanmedian(to):.0f} -> {np.nanmedian(tn):.0f}; '
              f'sessions with TP<0 > 50 %: {int((to > 50).sum())} -> {int((tn > 50).sum())}')
        for name in ('good', 'marginal', 'poor'):
            sel = [r for r in ok if r.get('cond_new') == name]
            if sel:
                print(f'  conditioning {name:8s}: {len(sel):2d} sessions  (|sin dphi_new| '
                      f'{min(1 / r["gain_new"] for r in sel):.2f}..{max(1 / r["gain_new"] for r in sel):.2f})')
        inc = [r['key'] for r in ok if not r.get('fit_consistent', True)]
        print(f'  first-sample vs line estimate differ > 20 deg (fit_consistent=False): {len(inc)} {inc}')
        if FIT == 'xecs':
            d = np.array([r.get('steve_minus_xecs') for r in ok if r.get('steve_minus_xecs') is not None], float)
            fl = np.array([r.get('F_lump') for r in ok if r.get('F_lump') is not None], float)
            rr = np.array([r['target'] / r['ratio_steve'] for r in ok if r.get('ratio_steve')], float)
            print(f'  Steve k0 angle − XeCS k0 angle: median {np.nanmedian(d):+.1f} deg, MAD {np.nanmedian(np.abs(d - np.nanmedian(d))):.1f}, '
                  f'|diff|>20: {int((np.abs(d) > 20).sum())}')
            print(f'  target ratio XeCS lumped / Steve area: median {np.nanmedian(rr):.2f}; F_lump median {np.nanmedian(fl):.2f} '
                  f'range {np.nanmin(fl):.2f}..{np.nanmax(fl):.2f}')
    for r in rows:
        if r.get('status') != 'resplit':
            print(f"  {r.get('status')}: {r['key']}  loose[{r.get('gate_loose')}] strict[{r.get('gate_strict')}] "
                  f"xecs[{r.get('gate_xecs', '-')}] {r.get('xecs_reason', '')}")


# ------------------------------------------------------------------------------- figure
def cmd_fig(keys):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import scipy.io as sio
    from scipy import ndimage
    sys.path.insert(0, str(HERE))
    import snr_calc as S
    (OUT / 'fig').mkdir(parents=True, exist_ok=True)
    for key in keys:
        rp = OUT / f'rows{SFX}' / f'{key}.json'
        if not rp.is_file() or json.load(open(rp)).get('status') != 'resplit':
            print(f'{key}: no resplit row'); continue
        row = json.load(open(rp))
        f = json.load(open(OUT / 'fits' / f'{key}.json'))
        m = sio.loadmat(DYN / key / 'd' / 'recon.mat',
                        variable_names=['gas_phase', 'dissolved_phase_real', 'dissolved_phase_imag'])
        n = sio.loadmat(DYN / key / 'd' / f'recon_resplit{SFX}.mat', variable_names=['aRBC', 'aTP', 'ph_new'])
        gas = m['gas_phase'].astype(float)
        Ro, To = m['dissolved_phase_real'].astype(float), m['dissolved_phase_imag'].astype(float)
        Rn, Tn = n['aRBC'].astype(float), n['aTP'].astype(float)
        b = row['insp_bin']
        corners = S.corner_box(gas.shape[1:]); gm = gas.mean(0)
        bg = corners & ~ndimage.binary_dilation(S.lung_mask(gm, gm[corners].std()), iterations=S.EXCL_DILATE)
        mk = S.lung_mask(gas[b], gas[b][bg].std())
        y = int(np.argmax(mk.sum(axis=(0, 2))))
        phR_old, phT_old = 2 * np.pi * f['fRBC'] * f['TEeff'], 2 * np.pi * f['fTP'] * f['TEeff']
        z = Ro[b] * np.exp(1j * phR_old) + To[b] * np.exp(1j * phT_old)
        mag = np.abs(z)
        # R(ph) curve with the new basis
        ph = np.arange(0, 2 * np.pi, 0.001)
        s = z[mk].sum() * np.exp(1j * ph)
        dn = np.radians(row['dphi_new'])
        sR = s.imag / np.sin(dn); sT = s.real - sR * np.cos(dn)
        Rc = np.where((sR > 0) & (sT > 0), sR / np.where(sT != 0, sT, np.nan), np.nan)

        fig, ax = plt.subplots(2, 4, figsize=(17, 8.2))
        lim_o = np.percentile(np.abs(np.concatenate([Ro[b][:, y, :].ravel(), To[b][:, y, :].ravel()])), 99.5)
        lim_n = np.percentile(np.abs(np.concatenate([Rn[b][:, y, :].ravel(), Tn[b][:, y, :].ravel()])), 99.5)
        panels = [(0, 0, Ro[b], lim_o, f'OLD aRBC  (dphi={row["dphi_old"]:.1f} deg, gain x{row["gain_old"]:.1f})'),
                  (0, 1, To[b], lim_o, f'OLD aTP   lung mean {row["lungTP_old_insp"]:.0f}, TP<0 {row["negTP_old_insp"]:.0f} %'),
                  (1, 0, Rn[b], lim_n, f'NEW aRBC  (dphi={row["dphi_new"]:.1f} deg, gain x{row["gain_new"]:.2f})'),
                  (1, 1, Tn[b], lim_n, f'NEW aTP   lung mean {row["lungTP_new_insp"]:.0f}, TP<0 {row["negTP_new_insp"]:.0f} %')]
        for i, j, img, lim, ttl in panels:
            a = ax[i, j]
            im = a.imshow(img[:, y, :], cmap='RdBu_r', vmin=-lim, vmax=lim)
            a.contour(mk[:, y, :], levels=[0.5], colors='k', linewidths=0.5)
            a.set_title(ttl, fontsize=9); a.set_xticks([]); a.set_yticks([])
            plt.colorbar(im, ax=a, fraction=0.046)
        a = ax[0, 2]; a.imshow(mag[:, y, :], cmap='gray'); a.contour(mk[:, y, :], levels=[0.5], colors='y', linewidths=0.5)
        a.set_title('|dissolved| (basis-independent)', fontsize=9); a.set_xticks([]); a.set_yticks([])
        a = ax[1, 2]; a.imshow(gas[b][:, y, :], cmap='gray'); a.set_title(f'gas bin {b} (mask source)', fontsize=9)
        a.set_xticks([]); a.set_yticks([])
        a = ax[0, 3]
        a.scatter(Ro[b][bg], To[b][bg], s=1, alpha=.3, label=f'old corr {row["corr_old_insp"]:+.3f}')
        a.scatter(Rn[b][bg], Tn[b][bg], s=1, alpha=.3, label=f'new corr {row["corr_new_insp"]:+.3f}')
        a.set_aspect('equal', 'datalim'); a.legend(fontsize=8); a.set_title('background aRBC vs aTP', fontsize=9)
        a.set_xlabel('aRBC'); a.set_ylabel('aTP')
        a = ax[1, 3]
        a.plot(ph, Rc, lw=1); a.axhline(row['target'], color='k', ls='--', lw=.8, label=f'target {row["target"]:.3f}')
        a.axvline(row['ph_new_insp'], color='tab:red', lw=.8, label=f'chosen ph {row["ph_new_insp"]:.2f} (R={row["R_new_insp"]:.3f})')
        a.set_yscale('log'); a.set_xlabel('global phase ph (rad)'); a.set_ylabel('masked R = sum aRBC / sum aTP')
        a.set_title('note-04 stop: min |R - target|, both sums > 0', fontsize=9); a.legend(fontsize=7)
        src = 'XeCS prior fit, lumped membrane' if FIT == 'xecs' else "Steve's fit"
        fig.suptitle(f'{key}  bin {b} (insp), coronal y={y}: RBC/TP re-split with the k0 basis angle from {src} '
                     f'(F59)  df={row["df"]:.0f} Hz, DPoff {row["DPoff"]:.0f} ppm, target {row["target"]:.3f}', fontsize=11)
        fig.tight_layout()
        p = OUT / 'fig' / f'{key}_before_after{SFX}.png'
        fig.savefig(p, dpi=110); plt.close(fig)
        print(f'fig -> {p}')


def cmd_cohort():
    """One-page cohort figure from resplit_summary.csv: basis angles old vs new, background
    correlation old vs new, lung TP<0 % old vs new."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = [r for r in csv.DictReader(open(OUT / f'resplit_summary{SFX}.csv')) if r['status'] == 'resplit']
    g = lambda k: np.array([float(r[k]) for r in rows])
    cond = np.array([r['cond_new'] for r in rows])
    col = np.where(cond == 'good', 'tab:green', 'tab:orange')
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    ax[0].scatter(g('dphi_old'), g('dphi_new'), c=col, s=14)
    ax[0].scatter(g('dphi_old'), g('dphi_line'), facecolors='none', edgecolors='gray', s=14, label='2πΔf·(TEeff+TE) (line)')
    ax[0].axhline(90, color='k', ls=':', lw=.8); ax[0].axhline(0, color='k', lw=.5); ax[0].axvline(0, color='k', lw=.5)
    ax[0].set_xlabel('Δφ used by Steve\'s split (deg)'); ax[0].set_ylabel('Δφ at k0 (deg)')
    ax[0].set_title(f'basis angle, {len(rows)} sessions (green |sin|≥0.7, orange marginal)', fontsize=9); ax[0].legend(fontsize=7)
    ax[1].scatter(g('corr_old_insp'), g('corr_new_insp'), c=col, s=14)
    ax[1].plot([-1, 1], [-1, 1], 'k--', lw=.6); ax[1].set_xlim(-1.02, 1.02); ax[1].set_ylim(-1.02, 1.02)
    ax[1].set_xlabel('old corr(aRBC, aTP) background'); ax[1].set_ylabel('new corr'); ax[1].set_title('noise anti-correlation (= −cos Δφ)', fontsize=9)
    ax[2].scatter(g('negTP_old_insp'), g('negTP_new_insp'), c=col, s=14)
    ax[2].plot([0, 100], [0, 100], 'k--', lw=.6); ax[2].set_xlabel('old % lung voxels with aTP < 0'); ax[2].set_ylabel('new %')
    ax[2].set_title('TP negative in lung (insp bin)', fontsize=9)
    fig.suptitle(f'Offline RBC/TP re-split of the DIAPHRAGM recons with the k0 basis angle (F59), fit = {FIT} — resplit_rbctp.py 2026-09-24')
    fig.tight_layout(); p = OUT / 'fig' / f'cohort_before_after{SFX}.png'; fig.savefig(p, dpi=120); print(f'fig -> {p}')


# ----------------------------------------------------------------------------------- cli
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('stage', choices=['fit', 'fit-one', 'split', 'summary', 'fig', 'cohort'])
    ap.add_argument('keys', nargs='*')
    ap.add_argument('--all', action='store_true', help='every session whose Tyger log split a bin')
    ap.add_argument('-j', '--jobs', type=int, default=6)
    ap.add_argument('--strict', action='store_true', help='note-05 gate numbers instead of the loose defaults')
    ap.add_argument('--dry', action='store_true', help='split: metrics only, no recon_resplit.mat')
    ap.add_argument('--fit', choices=['steve', 'xecs'], default='steve',
                    help='source of the k0 angle + ratio target: Steve raw.py refit (default) or the XeCS csv')
    ap.add_argument('--xecs-csv', default=None, help='path of the XeCS resplit_inputs csv (default: package copy)')
    a = ap.parse_args()
    global FIT, SFX
    FIT, SFX = a.fit, ('_xecs' if a.fit == 'xecs' else '')
    if a.fit == 'xecs':
        xecs_rows(a.xecs_csv)
    keys = all_split_keys() if a.all else a.keys
    if a.stage == 'fit':
        cmd_fit(keys, a.jobs)
    elif a.stage == 'fit-one':
        cmd_fit_one(keys[0])
    elif a.stage == 'split':
        cmd_split(keys, 'strict' if a.strict else 'loose', not a.dry)
    elif a.stage == 'summary':
        cmd_summary()
    elif a.stage == 'fig':
        cmd_fig(keys)
    elif a.stage == 'cohort':
        cmd_cohort()


if __name__ == '__main__':
    main()
