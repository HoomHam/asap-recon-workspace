"""Pneumotach clock alignment patch (applied 2026-10-27 on Hooman's order). Run from the repo root once.

Finding (PCA pneumo_clock_fit.csv, XeCS dps_offset F432/F433, draft §20/§23): the converter zeroes the pneumotach
clock at the file's first packet and raw.py assumes packet 0 == gas interleave 0. In truth the pneumotach log starts
at the record start (first cal FID; 26.7 s before gas interleave 0, plus a ~5.8 s term inside 2024-07-29..09-10), and
raw.py's ilvtime is a SYNTHETIC uniform iilv*TR grid (17.200 ms) while the true mean gas-interleave spacing on the
MDH clock is 17.37 ms. Every old-tree p/ recon is mis-binned (phase concentration 0.36 vs 0.95 on the true clock).

Fix (no free parameters): the converter reads the MDH timestamps (2.5 ms units) of every line and writes, per gas
interleave, its TRUE time since the first acquired line (= the pneumotach record start) as NdArray
'pneumo_gas_true_s'. raw.load_pneumotach_from_arrays maps the pneumotach clock onto the synthetic ilvtime by
piecewise-linear interpolation t_syn = interp(t_pne + offset, gas_true_s, ilvtime) (linear extrapolation at the
ends) BEFORE cropping / integrating. offset = header double 'pneumo_offset_s' (0 outside the window; dyn_recon passes
-5.8 inside 2024-07-29..09-10 per XeCS F78). Gate: header long 'pneumoclock' (1 = map, 0 = first-packet as before;
default 0 -> bit-exact old behaviour). Check: PCA's minima-to-EE fit (vendored as pneumoclock.py) runs on the mapped
trace against the SIGNAL EE list and is EXPORTED (pneumo_clock_chk_a/b/n_matched/resid_sd/ok) + the aligned trace
(pneumo_time_aligned / pneumo_vol_aligned) + T_cal (pneumo_clock_T_cal_s) + the mean true spacing
(pneumo_clock_ilv_true_ms).
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PCA_HELPER = Path('/Users/hoomham/Hooman/Work/Codes/2026_PCA_registration/workspace/helpers/scripts/_pneumo_clock_fit.py')


def patch(path, pairs):
    p = ROOT / path
    s = p.read_text()
    for old, new in pairs:
        n = s.count(old)
        assert n == 1, f'{path}: expected 1 match, got {n} for:\n{old[:160]}'
        s = s.replace(old, new)
    p.write_text(s)
    print(f'patched {path} ({len(pairs)} edits)')


# ---------------------------------------------------------------- pneumoclock.py (vendored check)
src = PCA_HELPER.read_text()
head, body = src.split('import numpy as np', 1)
vend = ('"""Pneumotach clock CHECK (vendored 2026-10-27 from the PCA session\'s helper\n'
        '2026_PCA_registration/workspace/helpers/scripts/_pneumo_clock_fit.py, author: 2026-pca-registration-04 Claude session,\n'
        '2026-10-26). In the fork the pneumotach clock is MAPPED onto the gas interleaves from the MDH timestamps\n'
        '(convert_siemens_to_mrd + raw.load_pneumotach_from_arrays); this fit only verifies the mapping: on a correctly\n'
        'mapped trace a ~ 1, b ~ 0 and >= 60 % of the k0 EE match within 0.3 s. Original docstring follows.\n'
        + head.split('"""', 2)[1].strip() + '\n"""\nimport numpy as np' + body.split("\n\nif __name__ == '__main__':")[0] + '\n')
(ROOT / 'pneumoclock.py').write_text(vend)
print('wrote pneumoclock.py')

# ---------------------------------------------------------------- gtypes.py
patch('gtypes.py', [
    ("    binpou = 0       # 1: normalise each interleave's bin weights to sum 1 (partition of unity)\n",
     "    binpou = 0       # 1: normalise each interleave's bin weights to sum 1 (partition of unity)\n"
     "    pneumoclock = 0  # 1: map the pneumotach clock onto the gas interleaves from the MDH timestamps (fork 2026-10-27)\n"
     "    pneumo_offset_s = 0.0  # added to the pneumotach clock before mapping (record-start offset, e.g. -5.8 s in the 2024-07..09 window)\n"),
])

# ---------------------------------------------------------------- convert_siemens_to_mrd.py
patch('convert_siemens_to_mrd.py', [
    ("    raw = twix.image.unsorted().astype('complex64')\n",
     "    raw = twix.image.unsorted().astype('complex64')\n"
     "    # MDH timestamps (2.5 ms units) of every line, acquisition order = unsorted() order; used to map the\n"
     "    # pneumotach clock onto the gas interleaves (the log starts at the first acquired line)\n"
     "    try:\n"
     "        ts_ms = np.asarray(twix.image.timestamp, dtype='float64').ravel() * 2.5\n"
     "    except Exception:\n"
     "        ts_ms = None\n"),
    ("                numspec_raw=0, dtspec_us=0.0)\n",
     "                numspec_raw=0, dtspec_us=0.0, ts_ms=ts_ms)\n"),
    ("        if pneumotach_file is not None:\n            pn = _parse_pneumotach(pneumotach_file)\n            yield mrd.StreamItem.NdArrayDouble(\n                mrd.NdArray(data=pn, meta={'pneumotach': [mrd.ArrayMetaValue.String('1')]}))\n",
     "        if pneumotach_file is not None:\n            pn = _parse_pneumotach(pneumotach_file)\n            yield mrd.StreamItem.NdArrayDouble(\n                mrd.NdArray(data=pn, meta={'pneumotach': [mrd.ArrayMetaValue.String('1')]}))\n"
     "            gts = _gas_true_times(meta.get('ts_ms'))\n"
     "            if gts is not None:\n"
     "                yield mrd.StreamItem.NdArrayDouble(\n"
     "                    mrd.NdArray(data=gts, meta={'pneumo_gas_true_s': [mrd.ArrayMetaValue.String('1')]}))\n"),
    ("def _parse_pneumotach(pneumo_file):",
     '''def _gas_true_times(ts_ms):
    """True time (s) of every GAS interleave since the first acquired line, from the MDH timestamps.
    The spectral (cal) block is the leading run of lines with spacing > 20 ms; the dynamic block alternates gas /
    dissolved with spacings ~7.5 / ~10 ms, gas first. Returns None when timestamps are unavailable."""
    if ts_ms is None or len(ts_ms) < 100:
        return None
    d = np.diff(ts_ms)
    dyn = np.where(d < 20.0)[0]
    if dyn.size == 0:
        return None
    nspec = int(dyn[0])                                  # first line of the dynamic block
    gas = np.arange(nspec, len(ts_ms), 2)
    t = (ts_ms[gas] - ts_ms[0]) * 1e-3
    print(f'pneumo clock: cal block {nspec} lines = {t[0]:.3f} s before gas interleave 0; '
          f'{len(gas)} gas interleaves, mean true spacing {1e3 * (t[-1] - t[0]) / max(len(gas) - 1, 1):.4f} ms', flush=True)
    return t.astype('float64')


def _parse_pneumotach(pneumo_file):'''),
    ("        mrd.UserParameterDoubleType(name='binsigma', value=float(params.get('binsigma', g_defaults.binsigma))),  # soft-bin kernel width (bins)\n",
     "        mrd.UserParameterDoubleType(name='binsigma', value=float(params.get('binsigma', g_defaults.binsigma))),  # soft-bin kernel width (bins)\n"
     "        mrd.UserParameterDoubleType(name='pneumo_offset_s', value=float(params.get('pneumo_offset_s', g_defaults.pneumo_offset_s))),\n"),
    ("        mrd.UserParameterLongType(name='binpou',     value=int(params.get('binpou', g_defaults.binpou))),        # partition-of-unity bin weights\n",
     "        mrd.UserParameterLongType(name='binpou',     value=int(params.get('binpou', g_defaults.binpou))),        # partition-of-unity bin weights\n"
     "        mrd.UserParameterLongType(name='pneumoclock', value=int(params.get('pneumoclock', g_defaults.pneumoclock))),  # map pneumotach clock from MDH timestamps\n"),
    ("    parser.add_argument('--binpou',     type=int,   default=gvar.binpou,     help='1: partition-of-unity bin weights')\n",
     "    parser.add_argument('--binpou',     type=int,   default=gvar.binpou,     help='1: partition-of-unity bin weights')\n"
     "    parser.add_argument('--pneumoclock', type=int, default=gvar.pneumoclock, help='1: map the pneumotach clock onto the gas interleaves from the MDH timestamps')\n"
     "    parser.add_argument('--pneumo-offset-s', type=float, default=gvar.pneumo_offset_s, dest='pneumo_offset_s',\n"
     "                        help='seconds added to the pneumotach clock before mapping (record-start offset)')\n"),
    ("                  binsigma=args.binsigma, binpou=args.binpou)\n",
     "                  binsigma=args.binsigma, binpou=args.binpou,\n                  pneumoclock=args.pneumoclock, pneumo_offset_s=args.pneumo_offset_s)\n"),
])

# ---------------------------------------------------------------- raw.py
patch('raw.py', [
    ("    def load_pneumotach_from_arrays(self, t, P):\n",
     "    def load_pneumotach_from_arrays(self, t, P, tmap=None, offset_s=0.0):\n"),
    ("        bt = bintype.PNEUMOTACH\n        t = np.asarray(t, dtype = 'float64').copy()\n        P = np.asarray(P, dtype = 'float64').copy()\n        if(len(P) > 51):\n",
     "        bt = bintype.PNEUMOTACH\n        t = np.asarray(t, dtype = 'float64').copy()\n        P = np.asarray(P, dtype = 'float64').copy()\n"
     "        self.pneumo_clock = {'mapped': 0, 'offset_s': float(offset_s)}\n"
     "        if(tmap is not None and len(tmap) == len(self.ilvtime) and len(tmap) > 2):\n"
     "            # map the pneumotach clock (0 = first acquired line = record start, + offset) onto the synthetic\n"
     "            # ilvtime: gas interleave j is at true time tmap[j] and synthetic time ilvtime[j]; piecewise-linear\n"
     "            # between interleaves, linear extrapolation beyond the ends (fork 2026-10-27, draft §20/§23)\n"
     "            tt = t + offset_s\n"
     "            tm = np.asarray(tmap, dtype = 'float64')\n"
     "            ts = np.asarray(self.ilvtime, dtype = 'float64')\n"
     "            t = np.interp(tt, tm, ts)\n"
     "            lo = tt < tm[0]\n"
     "            hi = tt > tm[-1]\n"
     "            t[lo] = ts[0] + (tt[lo] - tm[0]) * (ts[1] - ts[0]) / (tm[1] - tm[0])\n"
     "            t[hi] = ts[-1] + (tt[hi] - tm[-1]) * (ts[-1] - ts[-2]) / (tm[-1] - tm[-2])\n"
     "            self.pneumo_clock.update({'mapped': 1, 'T_cal_s': float(tm[0]),\n"
     "                                      'ilv_true_ms': float(1e3 * (tm[-1] - tm[0]) / (len(tm) - 1))})\n"
     "            print(f'pneumotach: clock mapped onto the gas interleaves (T_cal {tm[0]:.3f} s, offset {offset_s:+.2f} s, '\n"
     "                  f'true spacing {self.pneumo_clock[\"ilv_true_ms\"]:.4f} ms vs synthetic {1e3 * self.TR:.4f} ms)', file=stderr)\n"
     "        elif(tmap is not None):\n"
     "            print(f'pneumotach: clock map length {len(tmap)} != interleaves {len(self.ilvtime)} -- first-packet clock kept', file=stderr)\n"
     "        if(len(P) > 51):\n"),
    ("        if pneumo_arr is not None:\n            self.load_pneumotach_from_arrays(pneumo_arr[0], pneumo_arr[1])\n",
     "        if pneumo_arr is not None:\n"
     "            pmap = meta.get('pneumo_tmap') if isinstance(meta, dict) else None\n"
     "            poff = float(meta.get('pneumo_offset_s', 0.0)) if isinstance(meta, dict) else 0.0\n"
     "            self.load_pneumotach_from_arrays(pneumo_arr[0], pneumo_arr[1], tmap=pmap, offset_s=poff)\n"),
    ("                if(bt == bintype.PNEUMOTACH):\n                    p = np.polyfit(self.volmeasEEtime[bt] - np.mean(self.volmeasEEtime[bt]), minv, 8)\n",
     "                if(bt == bintype.PNEUMOTACH and len(self.volmeasEEtime[bintype.SIGNAL]) >= 3):\n"
     "                    # clock CHECK (pneumoclock.py, PCA helper): on the mapped trace a ~ 1, b ~ 0; exported, never applied\n"
     "                    try:\n"
     "                        from pneumoclock import fit_pneumo_clock\n"
     "                        ca, cb, info = fit_pneumo_clock(self.volmeastime[bt], self.volmeasvol[bt], self.volmeasEEtime[bintype.SIGNAL])\n"
     "                        self.pneumo_clock.update({'chk_a': ca, 'chk_b': cb, 'chk_n_matched': info['n_matched'],\n"
     "                                                  'chk_n_ee': info['n_ee'], 'chk_resid_sd_s': info['resid_sd_s'], 'chk_ok': int(info['ok'])})\n"
     "                        print(f'pneumotach: clock check vs k0 EE: a {ca:.5f} b {cb:+.2f} s, {info[\"n_matched\"]}/{info[\"n_ee\"]} EE matched, '\n"
     "                              f'resid {info[\"resid_sd_s\"]:.3f} s, ok {info[\"ok\"]}' + ('' if info['ok'] else '  <-- CLOCK NOT ALIGNED'), file=stderr)\n"
     "                    except Exception as e:\n"
     "                        print(f'pneumotach: clock check failed: {e!r}', file=stderr)\n"
     "                if(bt == bintype.PNEUMOTACH):\n                    p = np.polyfit(self.volmeasEEtime[bt] - np.mean(self.volmeasEEtime[bt]), minv, 8)\n"),
])

# ---------------------------------------------------------------- tyger_recon.py
patch('tyger_recon.py', [
    ("                elif item.value.meta.get('pneumotach'):\n                    pneumo_arr = item.value.data\n",
     "                elif item.value.meta.get('pneumotach'):\n                    pneumo_arr = item.value.data\n"
     "                elif item.value.meta.get('pneumo_gas_true_s'):\n                    pneumo_tmap = item.value.data\n"),
    ("    pneumo_arr = None\n    gas_phase_traj_arr = None\n",
     "    pneumo_arr = None\n    pneumo_tmap = None\n    gas_phase_traj_arr = None\n"),
    ("            'sigminN': int(user_long.get('sigminN', 50))}\n",
     "            'sigminN': int(user_long.get('sigminN', 50))}\n"
     "    if int(user_long.get('pneumoclock', 0)) and pneumo_tmap is not None:\n"
     "        meta['pneumo_tmap'] = np.asarray(pneumo_tmap, dtype='float64')\n"
     "        meta['pneumo_offset_s'] = float(user_double.get('pneumo_offset_s', 0.0))\n"),
    ("    binexp.update(_bin_export(g, g_raw, bt))\n",
     "    binexp.update(_bin_export(g, g_raw, bt))\n"
     "    pc = getattr(g_raw, 'pneumo_clock', None)\n"
     "    if pc:\n"
     "        for k, v in pc.items():\n"
     "            binexp[f'pneumo_clock_{k}'] = np.array([float(v)], dtype='float32')\n"
     "        if len(g_raw.volmeastime[bintype.PNEUMOTACH]):\n"
     "            binexp['pneumo_time_aligned'] = np.asarray(g_raw.volmeastime[bintype.PNEUMOTACH], dtype='float32')\n"
     "            binexp['pneumo_vol_aligned'] = np.asarray(g_raw.volmeasvol[bintype.PNEUMOTACH], dtype='float32')\n"),
])

# ---------------------------------------------------------------- workspace pipeline
patch('workspace/pipeline/dyn_recon.py', [
    ("    if ds.get('binpou') is not None:\n        cmd += ['--binpou', str(ds['binpou'])]\n",
     "    if ds.get('binpou') is not None:\n        cmd += ['--binpou', str(ds['binpou'])]\n"
     "    if ds.get('pneumoclock'):\n"
     "        # XeCS F78: inside 2024-07-29..2024-09-10 the pneumotach record starts ~5.8 s before the first cal FID\n"
     "        off = -5.8 if '2024-07-29' <= str(ds.get('date', '')) <= '2024-09-10' else 0.0\n"
     "        cmd += ['--pneumoclock', '1', '--pneumo-offset-s', str(off)]\n"
     "        _log(f'pneumotach clock: MDH-timestamp mapping ON, record-start offset {off:+.1f} s')\n"),
    ("    ap.add_argument('--binpou', type=int, default=None, help='1: partition-of-unity bin weights (fork MIX patch)')\n",
     "    ap.add_argument('--binpou', type=int, default=None, help='1: partition-of-unity bin weights (fork MIX patch)')\n"
     "    ap.add_argument('--pneumoclock', action='store_true', help='map the pneumotach clock onto the gas interleaves from the MDH timestamps (fork 2026-10-27)')\n"),
    ("                binsigma=args.binsigma, binpou=args.binpou)",
     "                binsigma=args.binsigma, binpou=args.binpou, pneumoclock=args.pneumoclock)"),
])
PK = ("'pneumo_clock_mapped', 'pneumo_clock_offset_s', 'pneumo_clock_T_cal_s', 'pneumo_clock_ilv_true_ms', 'pneumo_clock_chk_a', "
      "'pneumo_clock_chk_b', 'pneumo_clock_chk_n_matched', 'pneumo_clock_chk_n_ee', 'pneumo_clock_chk_resid_sd_s', "
      "'pneumo_clock_chk_ok', 'pneumo_time_aligned', 'pneumo_vol_aligned'")
patch('workspace/pipeline/mrd_to_mat.py', [
    ("'mix_tau_I', 'mix_n_insp', 'mix_tB_med', 'mix_T_med')",
     f"'mix_tau_I', 'mix_n_insp', 'mix_tB_med', 'mix_T_med',\n            {PK})"),
])
patch('workspace/pipeline/post_process.py', [
    ("'mix_tau_I', 'mix_n_insp', 'mix_tB_med', 'mix_T_med')",
     f"'mix_tau_I', 'mix_n_insp', 'mix_tB_med', 'mix_T_med',\n            {PK})"),
])
print('ALL PATCHES APPLIED')
