"""One-shot patch (2026-10-26, Hooman: "apply your mix approach for three subjects, call it mix-binning"):
adds a 4th respiratory binning mode MIX to the fork (branch diaphragm-recon) + two kernel knobs, all default-off.

MIX (E'ad in the decision draft, workspace/notes/Binning_Decision_Draft_2026-10-26.md §19):
  EE = Steve's SIGNAL k0 minima; per session tau_I = max(0.375*T_med, 1.25*tB_med) (tB = per-breath k0 peak time);
  n_insp = round(16*tau_I/T_med) clamped [4, 10]; bins 0..n_insp-1 on the absolute clock after EE; bins n_insp..15 on
  the navigator-volume RANK pooled over kept breaths (NaN nav -> -1; nav inside auto-excluded ranges set NaN first);
  breaths with T < 1.5*tau_I -> -1.
Kernel knobs (MRD header): binsigma (double, bins; default 1.0 == Steve's bindist0sq 2), binpou (long; 1 = each
interleave's bin weights normalised to sum 1). Export: ilvbin, ilvtime, ee_times, mix_tau_I, mix_n_insp, bin_weff
as NdArrays (alongside nav_*).
Run from the repo root. Idempotent-ish: every replacement asserts exactly one match.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def patch(path, pairs):
    p = ROOT / path
    s = p.read_text()
    for old, new in pairs:
        n = s.count(old)
        assert n == 1, f'{path}: expected 1 match, got {n} for:\n{old[:120]}'
        s = s.replace(old, new)
    p.write_text(s)
    print(f'patched {path} ({len(pairs)} edits)')


# ---------------------------------------------------------------- gtypes.py
patch('gtypes.py', [
    ("    DIAPHRAGM = 1\n    SIGNAL = 2\n",
     "    DIAPHRAGM = 1\n    SIGNAL = 2\n    MIX = 3          # inspiration on the absolute clock after EE, exhale on navigator rank (fork, 2026-10-26)\n"),
    ("    dplb = 40\n    freqfilter = 0\n",
     "    dplb = 40\n    freqfilter = 0\n"
     "    binsigma = 1.0   # soft-bin kernel width in bins: exp(-d^2 / (2 sigma^2)); 1.0 == Steve's bindist0sq = 2\n"
     "    binpou = 0       # 1: normalise each interleave's bin weights to sum 1 (partition of unity)\n"),
])

# ---------------------------------------------------------------- raw.py (4 bintypes)
patch('raw.py', [
    ("        self.ilvvol = [[], [], []]   # NORMALIZED lung volume at time of each interleave\n        self.ilvbin = [[], [], []]\n",
     "        self.ilvvol = [[], [], [], []]   # NORMALIZED lung volume at time of each interleave (one list per bintype)\n        self.ilvbin = [[], [], [], []]\n"),
    ("        self.volmeastime = [[], [], []]\n        self.volmeasvol = [[], [], []] # UNNORMALIZED\n        self.volmeasEEtime = [[], [], []]             # identified time of each EE\n",
     "        self.volmeastime = [[], [], [], []]\n        self.volmeasvol = [[], [], [], []] # UNNORMALIZED\n        self.volmeasEEtime = [[], [], [], []]             # identified time of each EE\n"),
])

# ---------------------------------------------------------------- recon.py (kernel width + partition of unity)
patch('recon.py', [
    ("def cudarecon(params, trajx, trajy, trajz, binarr, raw, k, knorm):",
     "def cudarecon(params, trajx, trajy, trajz, binarr, raw, k, knorm, binnorm):"),
    ("    smoothing = params[7]\n    binwt = 1.0\n",
     "    smoothing = params[7]\n"
     "    # params[8] = bin kernel width^2 scale (2 sigma^2; Steve's bindist0sq = 2 when absent), params[9] = 1 -> multiply\n"
     "    # by binnorm[ilv] (host-computed 1 / sum_b w_b so each interleave's weights sum to 1: partition of unity)\n"
     "    bd0sq = bindist0sq\n"
     "    pou = 0\n"
     "    if(params.shape[0] > 9):\n"
     "        bd0sq = params[8]\n"
     "        pou = params[9]\n"
     "    binwt = 1.0\n"),
    ("            binwt = exp(-bindist * bindist / bindist0sq)\n",
     "            binwt = exp(-bindist * bindist / bd0sq)\n"
     "            if(pou > 0):\n"
     "                binwt = binwt * binnorm[ilvnum]\n"),
])

# ---------------------------------------------------------------- results.py (params + call sites)
res = (ROOT / 'results.py').read_text()
# every cudarecon launch gets the extra device array; three single-image launches pass a dummy, dyn_recon the real one
n_calls = res.count('d_binarr, d_rawdata, self.d_kspace, self.d_kspacenorm)')
assert n_calls == 4, n_calls
res = res.replace('d_binarr, d_rawdata, self.d_kspace, self.d_kspacenorm)',
                  'd_binarr, d_rawdata, self.d_kspace, self.d_kspacenorm, d_binnorm)')
# single-image launches: define d_binnorm next to d_binarr (three sites use "d_binarr = cuda.to_device(np.zeros(1)) if(usegpu) else np.zeros(1)")
n_dummy = res.count('d_binarr = cuda.to_device(np.zeros(1)) if(usegpu) else np.zeros(1)\n')
assert n_dummy == 3, n_dummy
res = res.replace('d_binarr = cuda.to_device(np.zeros(1)) if(usegpu) else np.zeros(1)\n',
                  'd_binarr = cuda.to_device(np.zeros(1)) if(usegpu) else np.zeros(1)\n'
                  '        d_binnorm = d_binarr\n')
# dyn_recon: kernel params + binnorm from g.binsigma / g.binpou
old = ("        reconparams = np.array([g.MS, g.get_MScenter(), g_traj.nuniquesmp, g_raw.npts, g.nbins, 0, 0, 0])\n"
       "        binarr = np.ascontiguousarray(bins * g.nbins)\n")
new = ("        reconparams = np.array([g.MS, g.get_MScenter(), g_traj.nuniquesmp, g_raw.npts, g.nbins, 0, 0, 0])\n"
       "        binarr = np.ascontiguousarray(bins * g.nbins)\n"
       "        # soft-bin kernel: width g.binsigma (bins; 1.0 == Steve's bindist0sq 2) and optional partition of unity\n"
       "        # (g.binpou: each interleave's weights over the nbins bins sum to 1). Header knobs binsigma / binpou.\n"
       "        bd0sq = 2.0 * float(getattr(g, 'binsigma', 1.0)) ** 2\n"
       "        pou = int(getattr(g, 'binpou', 0))\n"
       "        binnorm = np.ones(len(binarr))\n"
       "        if(pou):\n"
       "            ibs = np.arange(g.nbins)[None, :]\n"
       "            bc = binarr[:, None]\n"
       "            bd = np.minimum(np.minimum(np.abs(ibs - bc), np.abs(ibs + g.nbins - bc)), np.abs(ibs - g.nbins - bc))\n"
       "            wsum = np.sum(np.exp(-bd * bd / bd0sq), axis=1)\n"
       "            binnorm = np.where((binarr >= 0) & (wsum > 0), 1.0 / np.maximum(wsum, 1e-12), 1.0)\n"
       "        binnorm = np.ascontiguousarray(binnorm)\n"
       "        d_binnorm = cuda.to_device(binnorm) if(usegpu) else binnorm\n"
       "        print(f'dyn_recon: bin kernel sigma {np.sqrt(bd0sq / 2):.2f} bins (bindist0sq {bd0sq:.2f}), partition of unity {pou}', file=stderr)\n")
assert res.count(old) == 1
res = res.replace(old, new)
old = ("                    reconparams = np.array([g.MS, g.get_MScenter(), g_traj.nuniquesmp, g_raw.npts, g.nbins, ibin, 0, \\\n"
       "                            g.gplb if(itype == imgtype.GPDYN) else g.dplb])\n")
new = ("                    reconparams = np.array([g.MS, g.get_MScenter(), g_traj.nuniquesmp, g_raw.npts, g.nbins, ibin, 0, \\\n"
       "                            g.gplb if(itype == imgtype.GPDYN) else g.dplb, bd0sq, pou])\n")
assert res.count(old) == 1, 'dyn_recon per-bin reconparams'
res = res.replace(old, new)
(ROOT / 'results.py').write_text(res)
print('patched results.py')

# ---------------------------------------------------------------- tyger_recon.py
patch('tyger_recon.py', [
    # header knobs
    ("    killpts = int(user_long.get('killpts', 2))\n",
     "    killpts = int(user_long.get('killpts', 2))\n"
     "    g.binsigma = float(user_double.get('binsigma', g.binsigma))\n"
     "    g.binpou = int(user_long.get('binpou', g.binpou))\n"),
    # navigator also for MIX
    ("    if binning == 'DIAPHRAGM' and not len(g_raw.ilvbin[bintype.DIAPHRAGM]):\n",
     "    if binning in ('DIAPHRAGM', 'MIX') and not len(g_raw.ilvbin[bintype.DIAPHRAGM]):\n"),
    # bin choice
    ("    if binning == 'DIAPHRAGM' and len(g_raw.ilvbin[bintype.DIAPHRAGM]):\n        bt = bintype.DIAPHRAGM\n",
     "    binexp = {}\n"
     "    if binning == 'MIX' and len(g_raw.ilvbin[bintype.DIAPHRAGM]):\n"
     "        g_raw.ilvbin[bintype.MIX], binexp = _mix_binning(g, g_raw)\n"
     "        bt = bintype.MIX\n"
     "    elif binning == 'DIAPHRAGM' and len(g_raw.ilvbin[bintype.DIAPHRAGM]):\n        bt = bintype.DIAPHRAGM\n"),
    ("    g_res.dyn_recon(g, g_raw, g_traj, g_raw.ilvbin[bt], g.usegpu)\n\n    print('writing reconstructed images to output mrd', file=stderr)\n    _write_results_to_mrd(g_res, header, output, nav=nav)\n",
     "    g_res.dyn_recon(g, g_raw, g_traj, g_raw.ilvbin[bt], g.usegpu)\n\n"
     "    # export the bins actually gridded (+ EE list, kernel weight per bin) so a reader can reproduce per-bin averages\n"
     "    binexp.update(_bin_export(g, g_raw, bt))\n"
     "    if nav is None:\n        nav = {}\n"
     "    nav.update(binexp)\n\n"
     "    print('writing reconstructed images to output mrd', file=stderr)\n    _write_results_to_mrd(g_res, header, output, nav=nav if len(nav) else None)\n"),
    # helpers before _write_results_to_mrd
    ("def _write_results_to_mrd(g_res, header, output, nav=None):",
     '''def _mix_binning(g, g_raw, split_min=0.375, tB_factor=1.25, n_min=4, n_max=10):
    """MIX binning (fork 2026-10-26): inspiration on the absolute clock after the SIGNAL EE, exhale on the navigator
    volume rank. Per session: tau_I = max(split_min * T_med, tB_factor * tB_med) with tB = time of the per-breath k0
    peak; n_insp = round(nbins * tau_I / T_med) clamped [n_min, n_max]. Bins 0..n_insp-1: phi = (n_insp/nbins) *
    (t - t_EE) / tau_I; bins n_insp..nbins-1: phi = n_insp/nbins + (1 - n_insp/nbins) * (rank + 0.5) / N, rank =
    descending navigator volume among all kept exhale interleaves (pooled over breaths). Breaths with T < 1.5 tau_I,
    interleaves outside the EE list, NaN navigator and auto-excluded interleaves -> -1 (navigator volume inside an
    auto-excluded range is set NaN first; the DIAPHRAGM navigator interpolates across such gaps)."""
    nb = g.nbins
    t = np.asarray(g_raw.ilvtime, dtype='float64')
    ee = np.asarray(g_raw.volmeasEEtime[bintype.SIGNAL], dtype='float64')
    k0 = np.asarray(g_raw.volmeasvol[bintype.SIGNAL], dtype='float64')
    nav = np.asarray(g_raw.ilvvol[bintype.DIAPHRAGM], dtype='float64').copy()
    bins = np.full(len(t), -1.0)
    out = {}
    if ee.size < 3:
        print('MIX: fewer than 3 SIGNAL EE minima -- no MIX bins', file=stderr)
        return bins, out
    excl = np.zeros(len(t), dtype=bool)
    for el in g_raw.excluderanges:
        excl |= (t >= el.start) & (t <= el.stop)
    nav[excl] = np.nan
    k = np.searchsorted(ee, t, side='right') - 1
    inside = (k >= 0) & (k < ee.size - 1)
    kk = np.where(inside, k, 0)
    tee = np.where(inside, t - ee[kk], np.nan)
    per = np.where(inside, ee[np.minimum(kk + 1, ee.size - 1)] - ee[kk], np.nan)
    T_med = float(np.median(np.diff(ee)))
    tB = []
    for j in range(ee.size - 1):
        m = inside & (k == j)
        if np.sum(m) >= 10:
            tB.append(float(tee[m][np.argmax(k0[m])]))
    tB_med = float(np.median(tB)) if len(tB) else 0.0
    tau = max(split_min * T_med, tB_factor * tB_med)
    n_insp = int(min(max(round(nb * tau / T_med), n_min), n_max))
    fi = n_insp / nb
    keep = inside & (per >= 1.5 * tau) & ~excl
    insp = keep & (tee < tau)
    bins[insp] = fi * tee[insp] / tau
    exh = keep & (tee >= tau) & np.isfinite(nav)
    if np.any(exh):
        v = nav[exh]
        r = np.empty(len(v))
        r[np.argsort(-v)] = np.arange(len(v))
        bins[exh] = fi + (1.0 - fi) * (r + 0.5) / len(v)
    bins = np.where(bins >= 0, np.clip(bins, 0.0, 1.0 - 1e-9), -1.0)
    print(f'MIX: T_med {T_med:.2f} s, tB_med {tB_med:.2f} s ({len(tB)} breaths), tau_I {tau:.2f} s, n_insp {n_insp}, '
          f'{int(np.sum(insp))} inspiratory + {int(np.sum(exh))} exhale interleaves binned, '
          f'{int(np.sum(~keep & inside))} dropped (short breath / excluded), {int(np.sum(exh | insp))} of {len(t)} total',
          file=stderr)
    out['mix_tau_I'] = np.array([tau], dtype='float32')
    out['mix_n_insp'] = np.array([n_insp], dtype='float32')
    out['mix_tB_med'] = np.array([tB_med], dtype='float32')
    out['mix_T_med'] = np.array([T_med], dtype='float32')
    return bins, out


def _bin_export(g, g_raw, bt):
    """Per-interleave bin coordinate as gridded, interleave times, SIGNAL EE list, effective kernel weight per bin."""
    from recon import bindist0sq
    bins = np.asarray(g_raw.ilvbin[bt], dtype='float64')
    nb = g.nbins
    bd0sq = 2.0 * float(getattr(g, 'binsigma', 1.0)) ** 2 if hasattr(g, 'binsigma') else bindist0sq
    ok = bins >= 0
    c = bins * nb
    w = np.zeros((len(bins), nb))
    for ib in range(nb):
        d = np.minimum(np.minimum(np.abs(ib - c), np.abs(ib + nb - c)), np.abs(ib - nb - c))
        w[:, ib] = np.where(ok, np.exp(-d * d / bd0sq), 0.0)
    if int(getattr(g, 'binpou', 0)):
        w = w / np.maximum(w.sum(1, keepdims=True), 1e-12)
    return {'ilvbin': bins.astype('float32'), 'ilvtime': np.asarray(g_raw.ilvtime, dtype='float32'),
            'ee_times': np.asarray(g_raw.volmeasEEtime[bintype.SIGNAL], dtype='float32'),
            'bin_weff': w.sum(0).astype('float32'),
            'bin_mode': np.array([int(bt)], dtype='float32')}


def _write_results_to_mrd(g_res, header, output, nav=None):'''),
])

# ---------------------------------------------------------------- workspace converter + driver + readers
patch('workspace/pipeline/convert_mixed_adc.py', [
    ("    p.add_argument('--binning', default='SIGNAL', choices=['SIGNAL', 'PNEUMOTACH', 'DIAPHRAGM'])",
     "    p.add_argument('--binning', default='SIGNAL', choices=['SIGNAL', 'PNEUMOTACH', 'DIAPHRAGM', 'MIX'])\n"
     "    p.add_argument('--binsigma', type=float, default=None, help='soft-bin kernel width in bins (fork MIX patch; 1.0 = Steve)')\n"
     "    p.add_argument('--binpou', type=int, default=None, help='1: partition-of-unity bin weights (fork MIX patch)')"),
    ("_EXTRA_LONG, _EXTRA_STR = {}, {}\n", "_EXTRA_LONG, _EXTRA_STR, _EXTRA_DBL = {}, {}, {}\n"),
    ("    up.user_parameter_string = list(up.user_parameter_string) + [\n        C.mrd.UserParameterStringType(name=k, value=str(v)) for k, v in _EXTRA_STR.items()]\n    if _EXTRA_LONG or _EXTRA_STR:\n        print(f'convert_mixed_adc: header extras {_EXTRA_LONG} {_EXTRA_STR}', flush=True)",
     "    up.user_parameter_string = list(up.user_parameter_string) + [\n        C.mrd.UserParameterStringType(name=k, value=str(v)) for k, v in _EXTRA_STR.items()]\n"
     "    up.user_parameter_double = list(up.user_parameter_double) + [\n        C.mrd.UserParameterDoubleType(name=k, value=float(v)) for k, v in _EXTRA_DBL.items()]\n"
     "    if _EXTRA_LONG or _EXTRA_STR or _EXTRA_DBL:\n        print(f'convert_mixed_adc: header extras {_EXTRA_LONG} {_EXTRA_STR} {_EXTRA_DBL}', flush=True)"),
    ("        _EXTRA_LONG['sigminN'] = a.sigmin_n\n",
     "        _EXTRA_LONG['sigminN'] = a.sigmin_n\n"
     "    if a.binsigma is not None:\n        _EXTRA_DBL['binsigma'] = a.binsigma\n"
     "    if a.binpou is not None:\n        _EXTRA_LONG['binpou'] = a.binpou\n"),
])
patch('workspace/pipeline/dyn_recon.py', [
    ("BINNING_LETTER = {'SIGNAL': 's', 'PNEUMOTACH': 'p', 'DIAPHRAGM': 'd'}",
     "BINNING_LETTER = {'SIGNAL': 's', 'PNEUMOTACH': 'p', 'DIAPHRAGM': 'd', 'MIX': 'm'}"),
    ("    if ds.get('sigmin_n') is not None:\n        cmd += ['--sigmin-n', str(ds['sigmin_n'])]\n",
     "    if ds.get('sigmin_n') is not None:\n        cmd += ['--sigmin-n', str(ds['sigmin_n'])]\n"
     "    if ds.get('binsigma') is not None:\n        cmd += ['--binsigma', str(ds['binsigma'])]\n"
     "    if ds.get('binpou') is not None:\n        cmd += ['--binpou', str(ds['binpou'])]\n"),
    ("    ap.add_argument('--sigmin-n', type=int, default=None,\n",
     "    ap.add_argument('--binsigma', type=float, default=None, help='soft-bin kernel width in bins (fork MIX patch; default Steve 1.0)')\n"
     "    ap.add_argument('--binpou', type=int, default=None, help='1: partition-of-unity bin weights (fork MIX patch)')\n"
     "    ap.add_argument('--sigmin-n', type=int, default=None,\n"),
    ("                merge=merge, ref=args.ref, excluderanges=args.excluderanges, sigmin_n=args.sigmin_n)",
     "                merge=merge, ref=args.ref, excluderanges=args.excluderanges, sigmin_n=args.sigmin_n,\n"
     "                binsigma=args.binsigma, binpou=args.binpou)"),
])
NEWKEYS = "'ilvbin', 'ilvtime', 'ee_times', 'bin_weff', 'bin_mode', 'mix_tau_I', 'mix_n_insp', 'mix_tB_med', 'mix_T_med'"
patch('workspace/pipeline/mrd_to_mat.py', [
    ("NAV_KEYS = ('nav_coronal', 'nav_diaphragm_z', 'nav_time', 'nav_volume', 'nav_ilvtime', 'nav_volmeastime')",
     f"NAV_KEYS = ('nav_coronal', 'nav_diaphragm_z', 'nav_time', 'nav_volume', 'nav_ilvtime', 'nav_volmeastime',\n            {NEWKEYS})"),
])
patch('workspace/pipeline/post_process.py', [
    ("NAV_KEYS = ('nav_coronal', 'nav_diaphragm_z', 'nav_time',\n            'nav_volume', 'nav_ilvtime', 'nav_volmeastime')",
     f"NAV_KEYS = ('nav_coronal', 'nav_diaphragm_z', 'nav_time',\n            'nav_volume', 'nav_ilvtime', 'nav_volmeastime',\n            {NEWKEYS})"),
])
print('ALL PATCHES APPLIED')
