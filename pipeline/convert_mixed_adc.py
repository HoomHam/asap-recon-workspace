"""Workspace wrapper of the root convert_siemens_to_mrd.py for twix files with MIXED ADC lengths.

2022 fancy_20220812 acqOrder 7 ('GP-DP-DP-Spectra') files interleave 1024-sample spectra (1 per 31 scans) with
512-sample imaging scans. mapVBVD (Python and MATLAB) sizes every scan to NCol = 1024 and returns garbage for
the 512-sample ones. Here the root converter is imported unchanged and only its reader `_read_twix` is
swapped: header metadata come from mapVBVD as before; raw data come from a direct MDH walk that keeps the
imaging scans (ushSamplesInScan == 512) in acquisition order and drops the spectra. Arm indexing downstream
(gas-scan counter mod nuniqueilvs) was verified on these files (2026-10-03, ASAP workspace
outputs/prev2_calib_2026-10-03): spectra do not advance the gas arm counter.

Same CLI as convert_siemens_to_mrd.py. dyn_recon.py uses it when ASAP_CONVERTER points here.
ASAP_DROP_CH=1[,..] or auto[:frac] (env; auto = drop coils with k0 SNR < frac x best, default 0.35) drops receive channels before writing (e.g. a near-dead coil that makes Steve's calcb
per-coil polynomial fit fail with < 10 mask voxels: 2022-10-26 MID00017, coil 1 k0/noise 7.7 vs 18-34).
"""
import argparse
import os
import pathlib
import sys

import numpy as np

PIPELINE = pathlib.Path(__file__).resolve().parent
REPO_ROOT = PIPELINE.parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(PIPELINE.parent / 'helpers'))
import convert_siemens_to_mrd as C          # noqa: E402  (root, read-only)
from gtypes import gvar                     # noqa: E402
from _prev2_static_recon import read_scans  # noqa: E402

_orig_read_twix = C._read_twix
NSAMP_IMG = 512


def _read_twix_mixed(dat_file):
    import mapvbvd
    twix = mapvbvd.mapVBVD(dat_file)
    if isinstance(twix, list):
        twix = twix[-1]
    meta = dict(TR_us=0.0, TE_us=0.0, DPoff=0.0, dtdyn_ns=0.0, numspec_raw=0, dtspec_us=0.0)
    y = twix.hdr.MeasYaps
    for key, k in [('TR_us', ('alTR', '0')), ('TE_us', ('alTE', '0')), ('DPoff', ('sWipMemBlock', 'adFree', '2')),
                   ('dtdyn_ns', ('sRXSPEC', 'alDwellTime', '0')), ('dtdyn_ns', ('sRXSPEC', 'alDwellTime', '1')),
                   ('dtspec_us', ('sWipMemBlock', 'alFree', '12'))]:
        try:
            meta[key] = float(y[k])
        except Exception:
            pass
    try:
        meta['numspec_raw'] = int(y[('sWipMemBlock', 'alFree', '10')] * y[('sWipMemBlock', 'alFree', '11')] + 0.1)
    except Exception:
        pass
    d, _, _, _ = read_scans(dat_file, nsamp_img=NSAMP_IMG)      # (nscan, nch, nsamp)
    spec = os.environ.get('ASAP_DROP_CH', '').strip()
    if spec.startswith('auto'):
        # per-coil SNR: mean |k0| (sample 1) over the strongest third of scans (gas-dominated) / per-coil
        # noise sd of the late readout (samples 400-511); drop coils below frac x the best coil.
        # Evidence (2022-10-26 HH fast): coil 1 at 0.23 x best lowered median gas SNR 9.1 -> 8.4 when kept.
        frac = float(spec.split(':')[1]) if ':' in spec else 0.35
        k0 = np.abs(d[:, :, 1])
        top = np.argsort(k0.sum(1))[-max(d.shape[0] // 3, 1):]
        snr = k0[top].mean(0) / np.abs(d[top][:, :, 400:]).std(axis=(0, 2))
        drop = [int(c) for c in np.nonzero(snr < frac * snr.max())[0]]
        print(f'convert_mixed_adc: per-coil k0 SNR {np.round(snr, 1).tolist()}; auto drop (< {frac} x max) -> {drop}',
              flush=True)
    else:
        drop = [int(c) for c in spec.split(',') if c.strip()]
    if drop and len(drop) < d.shape[1]:
        d = np.delete(d, drop, axis=1)
        print(f'convert_mixed_adc: ASAP_DROP_CH -> dropped channels {drop}', flush=True)
    raw = np.ascontiguousarray(d.transpose(2, 1, 0)).astype('complex64')   # (npts, nch, nilv) as mapvbvd unsorted
    print(f'convert_mixed_adc: {dat_file}: kept {raw.shape[2]} imaging scans of {NSAMP_IMG} samples, '
          f'{raw.shape[1]} ch (spectra dropped)', flush=True)
    return raw, meta


# ASAP_ROOT_READER=1 keeps the root reader (mapVBVD, sessions the root converter already handles, e.g. acq6) and only
# uses this wrapper for the header extras below, so the conversion stays identical to the root converter's
if os.environ.get('ASAP_ROOT_READER', '') != '1':
    C._read_twix = _read_twix_mixed

# ASAP_BSET=b0|b keeps one b-set of a 2022 'ADC' (diffusion) dynamic. Two interleaves exist (checked on k0
# 2026-10-21): 12-scan blocks b0 b b0 b b0 b | b b0 b b0 b b0, b0 <=> (s + s//6) even (F117; fancy_20220421/0422,
# 20221103 on 2023-11-09), and plain alternation (fancy_20220215/0323 ADC set-ups, 2022-02-24 000HH, 2022-04-20 017AK).
# The rule with the largest k0 b0/b contrast wins. The kept scans are consecutive in-set, so the downstream gas counter
# is the per-b-set arm counter and the pattern detector sees gas-only; TR (alTR = per-scan spacing on these, MDH clock)
# is doubled to the in-set spacing.
# ASAP_GASPHASE=3 keeps only the gas scans of a GP-DP-DP file (phase of 3 with the largest median k0) whose DP scans
# carry so much gas that the root detector calls it gas-only (2022-11-04 000HH: DP k0 at 0.84 x gas). alTR is the
# per-triplet spacing there (MDH: per-scan = alTR/3), so TR stays.
_BSET = os.environ.get('ASAP_BSET', '').strip()
_GASPHASE = int(os.environ.get('ASAP_GASPHASE', '0') or 0)
if _BSET or _GASPHASE:
    if _BSET not in ('', 'b0', 'b'):
        sys.exit(f'convert_mixed_adc: ASAP_BSET must be b0 or b, not {_BSET!r}')
    _inner_read_twix = C._read_twix

    def _read_twix_select(dat_file):
        raw, meta = _inner_read_twix(dat_file)
        if meta.get('numspec_raw', 0):
            sys.exit(f'convert_mixed_adc: scan selection with {meta["numspec_raw"]} initial spectra is not handled')
        s = np.arange(raw.shape[2])
        k0 = np.abs(raw[1]).sum(0)
        if _BSET:
            rules = {'blk12': (s + s // 6) % 2 == 0, 'alt0': s % 2 == 0, 'alt1': s % 2 == 1}
            ratios = {r: np.median(k0[m]) / np.median(k0[~m]) for r, m in rules.items()}
            rule = max(ratios, key=ratios.get)
            if ratios[rule] < 1.1:
                sys.exit(f'convert_mixed_adc: {dat_file}: no b0/b pattern (k0 b0/b {ratios})')
            keep = rules[rule] if _BSET == 'b0' else ~rules[rule]
            meta['TR_us'] *= 2
            what = f'ASAP_BSET={_BSET} rule {rule} (k0 b0/b {ratios[rule]:.2f})'
        else:
            med = [np.median(k0[p::_GASPHASE]) for p in range(_GASPHASE)]
            keep = s % _GASPHASE == int(np.argmax(med))
            what = f'ASAP_GASPHASE={_GASPHASE} phase {int(np.argmax(med))} (k0 by phase {np.round(np.array(med) / max(med), 2).tolist()})'
        raw = np.ascontiguousarray(raw[..., keep])
        print(f'convert_mixed_adc: {what}: kept {raw.shape[2]} of {len(s)} scans, TR {meta["TR_us"]:.0f} us', flush=True)
        return raw, meta

    C._read_twix = _read_twix_select

# extra MRD header user parameters read by tyger_recon.py from fork 262b16f on: 'excluderanges' (Steve's GUI
# exclude string, integer seconds on the gas-interleave clock iilv*TR) and 'sigminN' (SIGNAL EE-minimum half-window)
_EXTRA_LONG, _EXTRA_STR = {}, {}
_orig_write_header = C.mrd.BinaryMrdWriter.write_header


def _write_header_extra(self, header):
    up = header.user_parameters
    up.user_parameter_long = list(up.user_parameter_long) + [
        C.mrd.UserParameterLongType(name=k, value=int(v)) for k, v in _EXTRA_LONG.items()]
    up.user_parameter_string = list(up.user_parameter_string) + [
        C.mrd.UserParameterStringType(name=k, value=str(v)) for k, v in _EXTRA_STR.items()]
    if _EXTRA_LONG or _EXTRA_STR:
        print(f'convert_mixed_adc: header extras {_EXTRA_LONG} {_EXTRA_STR}', flush=True)
    return _orig_write_header(self, header)


C.mrd.BinaryMrdWriter.write_header = _write_header_extra

if __name__ == '__main__':
    p = argparse.ArgumentParser(description='convert_siemens_to_mrd with a mixed-ADC-length twix reader')
    p.add_argument('-i', '--input', required=True)
    p.add_argument('-o', '--output', required=True)
    p.add_argument('--gp-traj', default=None)
    p.add_argument('--dp-traj', default=None)
    p.add_argument('--pneumotach', default=None)
    p.add_argument('--seqname', default=None)
    p.add_argument('--ms', type=int, default=gvar.MS)
    p.add_argument('--is', type=int, default=gvar.IS, dest='IS')
    p.add_argument('--nbins', type=int, default=gvar.nbins)
    p.add_argument('--griddx', type=float, default=gvar.griddx)
    p.add_argument('--bindt', type=float, default=gvar.bindt)
    p.add_argument('--gplb', type=int, default=gvar.gplb)
    p.add_argument('--dplb', type=int, default=gvar.dplb)
    p.add_argument('--freqfilter', type=int, default=gvar.freqfilter)
    p.add_argument('--binning', default='SIGNAL', choices=['SIGNAL', 'PNEUMOTACH', 'DIAPHRAGM'])
    p.add_argument('--excluderanges', default=None, help="Steve's exclude ranges, e.g. '0-12, 66-93' (s, iilv*TR)")
    p.add_argument('--sigmin-n', type=int, default=None, help='SIGNAL EE-minimum half-window (Steve: 50)')
    a = p.parse_args()
    if a.excluderanges:
        _EXTRA_STR['excluderanges'] = a.excluderanges
    if a.sigmin_n is not None:
        _EXTRA_LONG['sigminN'] = a.sigmin_n
    params = dict(MS=a.ms, IS=a.IS, nbins=a.nbins, griddx=a.griddx, bindt=a.bindt, gplb=a.gplb, dplb=a.dplb,
                  freqfilter=a.freqfilter, binning=a.binning)
    C.convert_siemens_to_mrd(a.input, a.output, a.gp_traj, a.dp_traj, pneumotach_file=a.pneumotach, killpts=2,
                             params=params, seqname=a.seqname)
